library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(data.table)

# protein-level metadata 
protein_metadata <- read.table("data/derived/protein_metadata.tsv", header = TRUE) %>%
  mutate(gene_callers_id = as.character(gene_callers_id))

#------------------------------------------------------------
#  Generate taxonomy table for proteins (GTDB + NCBI nr)
#------------------------------------------------------------
rank_map <- c(
  d = "Domain",
  p = "Phylum",
  c = "Class",
  o = "Order",
  f = "Family",
  g = "Genus",
  s = "Species"
)

protein_taxonomy <- read_tsv("data/raw/CAT_output.ORF2LCA_GTDB.txt", guess_max = 1000000) |> 
  tidyr::separate(`# ORF`, into = c("c_prefix", "contig_num", "gene_callers_id"), sep = "_", remove = FALSE) |> 
  tidyr::unite("contig", c_prefix, contig_num, sep = "_", remove = TRUE) |> 
  filter(gene_callers_id %in% protein_metadata$gene_callers_id) |> 
  select(gene_callers_id, contig, lineage) |> 
  # handle NA lineage before splitting
  filter(!is.na(lineage)) |>                          # only process rows with lineage
  separate_longer_delim(lineage, delim = ";") |>
  filter(
    lineage != "",
    lineage != "root",
    !str_detect(lineage, "^\\d+$")
  ) |>
  extract(
    lineage, into = c("rank_letter", "name"),
    regex = "^([dpcofgs])__?(.*)$", remove = TRUE
  ) |>
  filter(!is.na(rank_letter)) |>
  mutate(rank = recode(rank_letter, !!!rank_map)) |>
  select(-rank_letter) |>
  pivot_wider(
    names_from  = rank,
    values_from = name,
    values_fn   = ~ .x[which(.x != "" & !is.na(.x))[1]],
    values_fill = NA_character_
  ) |>
  select(-any_of("NA")) |>
  relocate(any_of(c("Domain","Phylum","Class","Order","Family","Genus","Species")),
           .after = last_col()) |> 
  mutate(source="GDTB_LCA")


contig_taxonomy <- read_tsv("data/raw/CAT_output.contig2classification_GTDB.txt", guess_max = 1000000) |> 
  mutate(
    contig = as.character(`# contig`),
    source = "GDTB_contig",
    ORFs_on_contig = str_squish(str_remove(reason, "(?i)^\\s*based on\\s*"))
  ) %>%
  extract(
    reason,
    into = c("assigned_ORFs_on_contig", "total_ORFs_on_contig"),
    regex = "(?i)\\b(?:based on\\s*)?(\\d+)\\s*/\\s*(\\d+)\\s*ORFs\\b",
    remove = FALSE,
    convert = TRUE
  ) |> 
  mutate(Ratio=assigned_ORFs_on_contig/total_ORFs_on_contig) |> 
  filter(assigned_ORFs_on_contig>=2 & Ratio>=0.5) |> 
  select(contig, lineage) |> 
    # handle NA lineage before splitting
  filter(!is.na(lineage)) |>                          # only process rows with lineage
  separate_longer_delim(lineage, delim = ";") |>
  filter(
    lineage != "",
    lineage != "root",
    !str_detect(lineage, "^\\d+$")
  ) |>
  extract(
    lineage, into = c("rank_letter", "name"),
    regex = "^([dpcofgs])__?(.*)$", remove = TRUE
  ) |>
  filter(!is.na(rank_letter)) |>
  mutate(rank = recode(rank_letter, !!!rank_map)) |>
  select(-rank_letter) |>
  pivot_wider(
    names_from  = rank,
    values_from = name,
    values_fn   = ~ .x[which(.x != "" & !is.na(.x))[1]],
    values_fill = NA_character_
  ) |>
  select(-any_of("NA")) |>
  relocate(any_of(c("Domain","Phylum","Class","Order","Family","Genus","Species")),
           .after = last_col()) |> 
  mutate(source="GDTB_contig")


protein_taxonomy_GTDB <- protein_metadata |> select(gene_callers_id) |> left_join(protein_taxonomy |> select(gene_callers_id, contig), by = "gene_callers_id") |> left_join(contig_taxonomy, by = "contig") 
protein_taxonomy_GTDB<- protein_taxonomy_GTDB |> filter(is.na(Domain)) |> select(gene_callers_id, contig) |> left_join(protein_taxonomy, by = c("gene_callers_id","contig")) |> 
                        rbind(protein_taxonomy_GTDB |> filter(!is.na(Domain))) |> select(-contig) |> 
  mutate(across(where(is.character), ~ gsub("\\*", "", ., perl = TRUE)))

#complement missing taxonomy from blastp NCBI nr hits
protein_taxonomy_kaiju<-   read.table("data/raw/hits_proteins_kaiju_db_nr_2024-08-25_names.out", h=F, sep ="\t", fill = TRUE) |> 
  select(V2, V8) |> 
  tidyr::separate(V8, into = c("Domain","Phylum", "Order","Class", "Family","Genus", "Species"), sep = "; ",convert=TRUE) |> 
  rename(gene_callers_id=V2) |> mutate(gene_callers_id=as.character(gene_callers_id), source="kaiju_nr") |> 
  filter(gene_callers_id %in% protein_metadata$gene_callers_id) |> 
  select(gene_callers_id, Domain, Phylum, Class, Order, Family, Genus, Species, source) 

protein_taxonomy<- protein_taxonomy_GTDB |> filter(is.na(source)) |> select(gene_callers_id) |> left_join(protein_taxonomy_kaiju, by = "gene_callers_id") |> 
  rbind(protein_taxonomy_GTDB |> filter(!is.na(source))) |> 
  mutate(Order =case_when(
          Order =="Enterobacterales_A" ~ "Alteromonadales", 
          Order =="PCC-6307" ~ "Synechococcales", 
          Family=="Microbacteriaceae" ~ "Micrococcales",
          Family=="Maricaulaceae" ~ "Maricaulales",
          Family=="Acetobacteraceae" ~ "Rhodospirillales",
          Family=="Moraxellaceae" ~ "Moraxellales",
          Family=="Porticoccaceae" ~ "Cellvibrionales",
          Family=="Halieaceae" ~ "Cellvibrionales",
          Family=="Rhizobiaceae" ~ "Hyphomicrobiales",
          Family=="Streptomycetaceae" ~ "Kitasatosporales", TRUE ~ Order))

protein_taxonomy |> 
  write_tsv("data/derived/protein_taxonomy.tsv", col_names = TRUE)

#print session info and clean the workspace
sessionInfo()
rm(list = ls())

#unload all libraries loaded by this script
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
