# install.packages(c("readr","dplyr","tidyr","stringr","data.table"))

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(data.table)
library(ggplot2)

source("R/utils.R")

# protein-level metadata 
protein_metadata <- read.table("data/derived/protein_metadata.tsv", header = TRUE) %>%
  mutate(gene_callers_id = as.character(gene_callers_id))

#functional annotation using InterProScan
interpro_annotations<- lapply(c("InterPro", "Pfam","NCBIfam","Phobius","GOs"),
                                function(src){
                                annotation<- paste0(src,"_ann")
                                accession<- paste0(src,"_acc")
                                annotation_df <- read.table(paste0("data/raw/detected_proteins_",src,".tsv"),quote = "",  sep="\t", h=T) %>% 
                                  dplyr::rename(!!sym(annotation):= function.,
                                                !!sym(accession):= accession) %>% 
                                  select(-c("source")) %>% 
                                  group_by(gene_callers_id)%>%
                                  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
                                  summarise_all(funs(paste(unique(.), collapse='|')))
                                
                                return(annotation_df)
                              })  
#merge Interpro annotations into a table
interpro_annotations_df <- protein_metadata |> select(gene_callers_id) |> 
  left_join(interpro_annotations%>% purrr::reduce(full_join, by = "gene_callers_id") |>
  dplyr::rename("GO_terms"="GOs_acc") |> select(-GOs_ann) |> 
  tidyr::separate(gene_callers_id, into = c("c_prefix", "contig_num", "gene_callers_id"), sep = "_", remove = FALSE) %>%
  select(-c(c_prefix, contig_num)) %>%
  mutate(gene_callers_id = as.character(gene_callers_id)) , by = "gene_callers_id")

#import COG20 functions and categories
cog20_df <- read_tsv("data/raw/detected_proteins_COG_functions.txt", guess_max = 1000000, col_names = c("gene_callers_id", "source","COG20_acc", "COG20_ann", "evalue")) |> 
  mutate(COG20_acc=gsub("!!!.*", "", COG20_acc), COG20_ann=gsub("!!!.*", "", COG20_ann)) |> select(-c(source, evalue)) |> 
  left_join(read_tsv("data/raw/detected_proteins_COG_categories.txt", guess_max = 1000000, col_names = c("gene_callers_id", "source","COG20cat_acc", "COG20cat_ann", "evalue")) |> 
            mutate(COG20cat_acc=gsub("!!!.*", "", COG20cat_acc), COG20cat_ann=gsub("!!!.*", "", COG20cat_ann)) |> select(-c(source, evalue)), by = "gene_callers_id") |> 
  mutate(gene_callers_id = as.character(gene_callers_id))

#BLASTp against bacterial refseq
blastp_hits_df <- read.table("data/raw/detected_proteins_blastp_refseq.txt", h=T, sep ="\t") %>% 
  mutate(gene_callers_id =as.character(gene_callers_id)) %>%  select(gene_callers_id,blastp_acc,blastp_ann) |> 
  tidyr::separate(gene_callers_id, into = c("c_prefix", "contig_num", "gene_callers_id"), sep = "_", remove = FALSE) %>%
  select(-c(c_prefix, contig_num)) %>%
  mutate(gene_callers_id = as.character(gene_callers_id)) 


#PSIBLAST against TCDB
TCDB_id2fam <- readr::read_tsv(
  "data/raw/TCDB_id2fam.txt",
  skip = 1,
  col_names = c("TCID", "TC_subfamily", "TC_family", "Fam_abbreviation", "Superfamily"),
  guess_max = 1000000
) |> 
  left_join(readr::read_tsv("data/raw/TCDB_fam_defs.txt", guess_max = 1000000, col_names = c("TC_family", "TC_description")) , by = "TC_family") |> 
  left_join(readr::read_tsv("data/raw/TCDB_substrate.txt", guess_max = 1000000, col_names = c("TCID", "Substrate")) , by = "TCID") |> 
  mutate(TCDB_fam=paste(TC_description, Substrate, sep = "_"), TCDB_acc=TCID) 

TCDB_substrate <- readr::read_tsv("data/raw/TCDB_substrate.txt", guess_max = 1000000, col_names = c("TCID", "Substrate")) |>
  mutate(
    # strip CHEBI IDs from substrate, keep only the human-readable name
    substrate_clean = gsub("CHEBI:\\d+;", "", Substrate),
    substrate_clean = gsub("\\|", " | ", substrate_clean))

TCDB_hits_df <- read.table("data/raw/detected_proteins_TCDB_ref.txt", h=T, sep ="\t") |> 
                    group_by(gene_callers_id) %>% 
                    filter(pident>30  & evalue<0.001) %>% 
                    filter(bitscore >50) |> 
                    slice_max(pident, n = 1, with_ties=FALSE) |> 
                    tidyr::separate(sseqid, into = c("c_prefix", "source", "refseq","TCID"), sep = "\\|", remove = TRUE) %>%
                    mutate(TCDB_acc = str_extract(TCID, "^[^.]+\\.[^.]+\\.[^.]+")) |> 
                    left_join(readr::read_tsv("data/raw/TCDB_fam_defs.txt", guess_max = 1000000, col_names = c("TC_family", "TC_description")) , by = c("TCDB_acc"="TC_family")) |> 
                    select(gene_callers_id, TCID, TCDB_ann=function., TCDB_acc, TC_description) |> 
                    left_join(TCDB_substrate, by = "TCID") |> 
                    tidyr::separate(gene_callers_id, into = c("c_prefix", "contig_num", "gene_callers_id"), sep = "_", remove = FALSE) %>%
                    select(-c(c_prefix, contig_num)) %>%
                    mutate(gene_callers_id = as.character(gene_callers_id)) 

   
#deeploc localization estimates
source("R/deeplocpro_filtration.R")

#signal peptide predictions
source("R/signalP_filtration.R")

#merge all annotations
protein_annotations_df <- protein_metadata |> select(gene_callers_id) |>
  merge(interpro_annotations_df,  by = "gene_callers_id", all = TRUE) |>
  merge(blastp_hits_df,            by = "gene_callers_id", all = TRUE) |>
  merge(TCDB_hits_df,              by = "gene_callers_id", all = TRUE) |>
  merge(cog20_df , by = "gene_callers_id", all = TRUE) |> 
  merge(deeplocpro_filt|> select(gene_callers_id, deeplocpro_localization, deeplocpro_confidence), by = "gene_callers_id", all=TRUE) |> 
  merge(signalp_filt|> select(gene_callers_id, SignalP_Prediction, SignalP_confidence), by = "gene_callers_id", all=TRUE) 

#final sub-cellular localization prediction
source("R/protein_localization.R")

localization_df<- read_tsv("data/derived/protein_localization.tsv", guess_max = 1000000) |> mutate(gene_callers_id = as.character(gene_callers_id))

protein_annotations <- protein_metadata |> select(gene_callers_id) |>
  left_join(interpro_annotations_df,  by = "gene_callers_id") |> 
  left_join(blastp_hits_df,            by = "gene_callers_id") |>
  left_join(TCDB_hits_df,              by = "gene_callers_id") |>
  left_join(cog20_df , by = "gene_callers_id") |> 
  left_join(localization_df|> select(gene_callers_id, Localization_final, Localization_source, reliability_score) , by = "gene_callers_id") 

protein_annotations |> resolve_function() |>
  write_tsv("data/derived/protein_annotations.tsv", col_names = TRUE)

#print session info and clean the workspace
sessionInfo()
rm(list = ls())

#unload all libraries loaded by this script
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
