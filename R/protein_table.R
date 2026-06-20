# install.packages(c("readr","dplyr","tidyr","stringr","data.table"))

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(data.table)
library(ggplot2)

source("R/utils.R")

#------------------------------------------------------------
# 0) Inputs and basic settings
#------------------------------------------------------------
min_total_pep  <- 2  # ≥1 peptides total (Unique + Razor)
min_unique_pep <- 1  # ≥1 strictly unique peptide
include_razor_in_total <- TRUE

#------------------------------------------------------------
# 1) Read inputs + your sample mapping
#------------------------------------------------------------
pep  <- read_tsv("data/raw/SO289_nr100_PeptideGroups.txt",  guess_max = 100000)
prot <- read_tsv("data/raw/SO289_nr100_Proteins.txt", guess_max = 100000)

samples_df <- read.table("data/samples_metadata.txt", header = TRUE) |> 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 

stopifnot(all(c("File.ID","Sample_ID") %in% names(samples_df)))
if (!"Type" %in% names(samples_df)) {
  samples_df$Type <- ifelse(grepl("EVs$", samples_df$Sample_ID, ignore.case = TRUE), "EVs", "Cells")
}

#------------------------------------------------------------
# 2) Peptide-level filter and in-sample presence
#------------------------------------------------------------
# Peptide q-value column (prefer Percolator; fallback Qvality)
qcol_pep <- dplyr::case_when(
  "Percolator q-Value by Search Engine Sequest HT" %in% names(pep) ~
    "Percolator q-Value by Search Engine Sequest HT",
  "Qvality q-value" %in% names(pep) ~
    "Qvality q-value",
  TRUE ~ NA_character_
)
stopifnot(!is.na(qcol_pep))

# Protein keys
prot_key_in_pep  <- if ("Master Protein Accessions" %in% names(pep)) "Master Protein Accessions" else "Protein Accessions"
stopifnot(prot_key_in_pep %in% names(pep))
prot_key_in_prot <- "Accession"; stopifnot(prot_key_in_prot %in% names(prot))

# Found-in-sample columns in Peptide Groups
found_cols_pep <- grep("^Found in Sample in\\s", names(pep), value = TRUE)
stopifnot(length(found_cols_pep) > 0)

# Map "Found in Sample" strings to logical
found_map <- function(x) {
  # Treat "High" and "Peak Found" as TRUE, "Not Found" (or NA) as FALSE
  !(is.na(x) | x == "Not Found")
}

# Filter peptides at 1% FDR
pep_filt <- pep %>%
  mutate(q_pep = suppressWarnings(as.numeric(.data[[qcol_pep]]))) %>%
  filter(!is.na(q_pep) & q_pep <= 0.01) %>%
  mutate(
    # Use protein mapping to determine uniqueness
    n_groups = suppressWarnings(as.integer(`Number of Protein Groups`)),
    n_prots  = suppressWarnings(as.integer(`Number of Proteins`)),
    
    # Unique/shared flags
    is_unique = n_groups == 1L,
    is_shared = n_groups  > 1L,

     counts_for_total = if (include_razor_in_total) (is_unique | is_shared) else is_unique
  )

# Long format of peptide presence per file (parse File.ID from column name)
extract_file_id <- function(col_name) {
  m <- regexpr("\\bF\\d+\\b", col_name)
  ifelse(m > 0, regmatches(col_name, m), NA_character_)
}

pep_present_long <- pep_filt %>%
  select(pep_id = `Annotated Sequence`, protein = all_of(prot_key_in_pep), Modifications,
         all_of(found_cols_pep), is_unique, counts_for_total) %>%
  filter(!is.na(protein) & protein != "") %>%
  pivot_longer(cols = all_of(found_cols_pep), names_to = "found_col", values_to = "found_str") %>%
  mutate(found = found_map(found_str),
         File.ID = extract_file_id(found_col)) %>%
  filter(found, !is.na(File.ID)) %>%
  left_join(samples_df[, c("File.ID","Sample_ID","Type")], by = "File.ID")

# Clean peptide base sequence to avoid double counting
strip_flanks <- function(s) {
  s <- sub("^\\[[A-Z]\\]\\.", "", s)
  s <- sub("\\.\\[[A-Z]\\]$", "", s)
  s
}
pep_present_long <- pep_present_long %>%
  mutate(base_seq = strip_flanks(pep_id))

#------------------------------------------------------------
# 3) Count per protein × Sample_ID: unique and total peptides
#------------------------------------------------------------
pep_counts_ps <- pep_present_long %>%
  group_by(Sample_ID, protein, base_seq) %>%
  summarise(any_unique = any(is_unique),
            any_total  = any(counts_for_total),
            .groups = "drop") %>% 
  tidyr::separate_rows(protein, sep= "; ") %>%
  group_by(Sample_ID, protein) %>%
  summarise(n_unique = sum(any_unique),
            n_total  = sum(any_total),
            .groups = "drop")

# Apply per-sample rule: ≥1 unique AND ≥2 total
prot_pass_ps <- pep_counts_ps %>%
  mutate(pass = (n_unique >= min_unique_pep) & (n_total >= min_total_pep))

#------------------------------------------------------------
# 4) Filter Proteins.txt to 1% protein FDR and align keys
#------------------------------------------------------------
prot_fdr <- prot
if ("Exp q-value Combined" %in% names(prot_fdr)) {
  prot_fdr <- prot_fdr %>%
    mutate(q_prot = suppressWarnings(as.numeric(`Exp q-value Combined`))) %>%
    filter(!is.na(q_prot) & q_prot <= 0.01)
} else if ("Protein FDR Confidence Combined" %in% names(prot_fdr)) {
  prot_fdr <- prot_fdr %>% filter(`Protein FDR Confidence Combined` == "High")
} else {
  warning("No protein-level FDR column (using unfiltered Proteins.txt).")
}

# Coerce keys to character for join
prot_ids_fdr <- prot_fdr %>%
  filter(Master == "IsMasterProtein") %>%
  transmute(Accession = as.character(!!sym(prot_key_in_prot))) %>%
  distinct()

prot_pass_ps_fdr <- prot_pass_ps %>%
  mutate(Accession = as.character(protein)) %>%
  inner_join(prot_ids_fdr, by = "Accession") %>%
  filter(pass)

# Per-sample counts of qualifying proteins
proteins_identified_per_sample <- prot_pass_ps_fdr %>%
  group_by(Sample_ID) %>%
  summarise(n_proteins = n_distinct(Accession), .groups = "drop") %>%
  arrange(Sample_ID)

#------------------------------------------------------------
# 5) Build pass/fail matrix (protein × sample) for auditing
#------------------------------------------------------------
pass_matrix <- prot_pass_ps %>%
  mutate(Accession = as.character(protein)) %>%
  inner_join(prot_ids_fdr, by = "Accession") %>%
  select(Accession, Sample_ID, n_unique, n_total, pass)

pass_wide <- pass_matrix %>%
  select(Accession, Sample_ID, pass) %>%
  pivot_wider(names_from = Sample_ID, values_from = pass, values_fill = FALSE)

#------------------------------------------------------------
# 6) Mask Proteins abundances by per-sample evidence rule
#     (rename Abundance F# Sample -> Sample_ID using samples_df)
#------------------------------------------------------------
ab_renamed <- rename_abundance_cols(prot, samples_df) 
prot_abund_named <- ab_renamed$df
abundance_cols   <- ab_renamed$sample_cols

pass_lookup <- pass_matrix %>%
    select(Accession, Sample_ID, pass)

# Long abundance table, mask by pass
prot_abund_long <- prot_abund_named %>% 
    mutate(Accession = as.character(!!sym(prot_key_in_prot))) %>%
    select(Accession, all_of(abundance_cols)) %>%
    pivot_longer(cols = all_of(abundance_cols), names_to = "Sample_ID", values_to = "abundance") %>%
    left_join(pass_lookup, by = c("Accession","Sample_ID")) %>%
    mutate(abundance_masked = ifelse(pass, abundance, NA_real_))

prot_abund_masked <- prot_abund_long |>
  select(Accession, Sample_ID, abundance_masked) |>
  pivot_wider(names_from  = Sample_ID,
              values_from = abundance_masked) |>
  # remove proteins where ALL sample columns are NA
  filter(if_any(all_of(abundance_cols), ~ !is.na(.x)))

cat(sprintf(
  "Proteins retained after removing all-NA rows: %d / %d\n",
  nrow(prot_abund_masked),
  n_distinct(prot_abund_long$Accession)
))

  # Rebuild a full Proteins table with masked abundances
prot_masked_full <- prot_abund_named %>%
    mutate(Accession = as.character(!!sym(prot_key_in_prot))) %>%
    select(-all_of(abundance_cols)) %>%
    left_join(prot_abund_masked, by = "Accession") |> 
      # remove proteins where ALL sample columns are NA
    filter(if_any(all_of(abundance_cols), ~ !is.na(.x)))


#------------------------------------------------------------
# 7) Write protein table
#------------------------------------------------------------
prot_masked_full |> 
  select(Accession) |> 
  write_tsv("data/derived/detected_protein_gcids.tsv", col_names = FALSE)

prot_masked_full |> 
  mutate(
    gene_callers_id = as.character(Accession)
  ) |> 
  select(gene_callers_id, all_of(abundance_cols)) |> 
  write_tsv("data/derived/protein_abundance.tsv", col_names = TRUE)

#------------------------------------------------------------
# 8) Generate protein metadata table
#------------------------------------------------------------

protein_sequences<- read_tsv("data/raw/protein_seqs.txt",  guess_max = 100000, col_names = c("gene_callers_id","AA_sequence")) |> 
  filter(gene_callers_id %in% prot_masked_full$Accession) |> 
  mutate(gene_callers_id = as.character(gene_callers_id))


protein_metadata<- prot_masked_full |> 
  mutate(
    gene_callers_id = as.character(Accession),
    Length = as.integer(`Number of AAs`),
    MW_kDa = as.numeric(`MW in kDa`),
    Unique_pep     = as.integer(`Number of Unique Peptides`),
    Total_pep    = as.integer(`Number of Peptides`),
    PSMs    = as.integer(`Number of PSMs`)
  ) |> 
  left_join(protein_sequences, by = "gene_callers_id") |>
  select(gene_callers_id, Length, MW_kDa, Unique_pep, Total_pep, PSMs, AA_sequence)

protein_metadata |>
  write_tsv("data/derived/protein_metadata.tsv", col_names = TRUE)


#print session info and clean the workspace
sessionInfo()
rm(list = ls())

#unload all libraries loaded by this script
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
