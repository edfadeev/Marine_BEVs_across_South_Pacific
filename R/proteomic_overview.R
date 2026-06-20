require(dplyr)
require(DEqMS)
require(vegan)
require(tibble)
require(readr)

# load custom ggplot theme and colour palettes
source("R/utils.R")

############################
#import samples metadata
############################

samples_df<- read.table("data/samples_metadata.txt", header = TRUE) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 

############################
#import derived data tables
############################
# protein abundance table
run_source_if_missing("data/derived/protein_abundance.tsv", "R/protein_table.R")
protein_abund <- read.table("data/derived/protein_abundance.tsv", header = TRUE) |> 
  tibble::column_to_rownames("gene_callers_id")

# taxonomic classification per protein 
run_source_if_missing("data/derived/protein_taxonomy.tsv", "R/protein_taxonomy.R")
protein_taxonomy <- read.table("data/derived/protein_taxonomy.tsv", header = TRUE, sep = "\t") %>%
  mutate(gene_callers_id = as.character(gene_callers_id))

# functional annotations per protein 
run_source_if_missing("data/derived/protein_annotations.tsv", "R/protein_annotation.R")
protein_annotations <- read.csv("data/derived/protein_annotations.tsv", sep = "\t") %>%
  mutate(gene_callers_id = as.character(gene_callers_id))

# protein-level metadata
protein_metadata <- read.table("data/derived/protein_metadata.tsv", header = TRUE) %>%
  mutate(gene_callers_id = as.character(gene_callers_id))

############################
#Number of identified proteins
############################
#bacterial 
bac_gcids_taxa<- protein_taxonomy %>% 
  filter(Domain=="Bacteria")%>% pull(gene_callers_id)

# keep only bacterial IDs that actually exist as rownames and warn if any are missing
bact_ids <- as.character(bac_gcids_taxa)
present_ids <- intersect(rownames(protein_abund), bact_ids)

# safe subset (preserve matrix shape)
protein_abund_bac <- protein_abund[present_ids, , drop = FALSE]

# convert to data frame and make gene_callers_id an explicit column
protein_abund_bac_df <- protein_abund_bac %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_callers_id")

# then continue your pipeline using protein_abund_bac_df
prot_per_sample <- protein_abund_bac_df %>%
  reshape2::melt(id.vars = "gene_callers_id",
                 variable.name = "Sample_ID",
                 value.name = "Abund") %>%
  filter(Abund > 0) %>%
  left_join(samples_df, by = "Sample_ID") %>%
  group_by(Region, Station_ID, Fraction) %>%
  dplyr::summarize(N_prot = n(), .groups = "drop") %>%
  tidyr::spread(Fraction, N_prot) %>%
  dplyr::rename("Cells_Bac" = "Cells", "BEVs_Bac" = "BEVs")


#identify viral proteins
vir_gcids_taxa<- protein_taxonomy %>% 
  filter(grepl("Viruses", Domain))%>% pull(gene_callers_id)

vir_gcids_ann <-protein_annotations %>% 
  filter(grepl("phage|virus|capsid|Tail sheath",Pfam_ann, ignore.case = TRUE)|
           grepl("phage|virus|capsid|Tail sheath",InterPro_ann, ignore.case = TRUE)|
           NCBIfam_acc %in% c("TIGR01554", "TIGR02126", "TIGR01543")|Function=="virion structural protein") %>% pull(gene_callers_id)
vir_gcids<- c(vir_gcids_taxa, vir_gcids_ann)

#summarize how many viral proteins per sample
# keep only viral IDs that actually exist as rownames and warn if any are missing
vir_ids <- as.character(vir_gcids)
present_ids <- intersect(rownames(protein_abund), vir_ids)

# safe subset (preserve matrix shape)
protein_abund_vir <- protein_abund[present_ids, , drop = FALSE]

# convert to data frame and make gene_callers_id an explicit column
protein_abund_vir_df <- protein_abund_vir %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_callers_id")

viral_prot_per_sample <- protein_abund_vir_df %>%
  reshape2::melt(id.vars = "gene_callers_id",
                 variable.name = "Sample_ID",
                 value.name = "Abund") %>%
  filter(Abund > 0) %>%
  left_join(samples_df, by = "Sample_ID") %>%
  group_by(Region, Station_ID, Fraction) %>%
  dplyr::summarize(N_prot = n(), .groups = "drop") %>%
  tidyr::spread(Fraction, N_prot) %>%
  dplyr::rename("Cells_Vir" = "Cells", "BEVs_Vir" = "BEVs")


#unclassified 
uncl_gcids_taxa<- protein_taxonomy %>% 
  filter(is.na(Phylum) | Phylum=="")%>% pull(gene_callers_id)

uncl_ids <- as.character(uncl_gcids_taxa)
present_ids <- intersect(rownames(protein_abund), uncl_ids)

# safe subset (preserve matrix shape)
protein_abund_uncl <- protein_abund[present_ids, , drop = FALSE]

# convert to data frame and make gene_callers_id an explicit column
protein_abund_uncl_df <- protein_abund_uncl %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_callers_id")

uncl_prot_per_sample <- protein_abund_uncl_df %>%
  reshape2::melt(id.vars = "gene_callers_id",
                 variable.name = "Sample_ID",
                 value.name = "Abund") %>%
  filter(Abund > 0) %>%
  left_join(samples_df, by = "Sample_ID") %>%
  group_by(Region, Station_ID, Fraction) %>%
  dplyr::summarize(N_prot = n(), .groups = "drop") %>%
  tidyr::spread(Fraction, N_prot) %>%
  dplyr::rename("Cells_uncl" = "Cells", "BEVs_uncl" = "BEVs")

#merged table and write into file
protein_abund %>%
   as.data.frame() %>%
  tibble::rownames_to_column("gene_callers_id") |> 
  reshape2::melt(id.vars = "gene_callers_id",
                 variable.name = "Sample_ID",
                 value.name = "Abund") %>%
  filter(Abund > 0) %>%
  left_join(samples_df, by = "Sample_ID") %>%
  group_by(Region, Station_ID, Fraction) %>%
  dplyr::summarize(N_prot = n(), .groups = "drop") %>%
  tidyr::spread(Fraction, N_prot) %>%
  dplyr::rename("Cells_Total" = "Cells", "BEVs_Total" = "BEVs") |> 
  left_join(prot_per_sample, by =c("Station_ID", "Region")) %>% 
  left_join(viral_prot_per_sample, by =c("Station_ID", "Region")) %>% 
  left_join(uncl_prot_per_sample, by =c("Station_ID", "Region")) %>% 
  arrange(rev(Station_ID)) %>% 
  write.table(file="Tables/Table_S1-Proteomics_overview.txt", col.names =TRUE, row.names = FALSE, quote = FALSE)


#proportions of viral proteins
protein_abund %>%
   as.data.frame() %>%
  tibble::rownames_to_column("gene_callers_id") |> 
  reshape2::melt(id.vars = "gene_callers_id",
                 variable.name = "Sample_ID",
                 value.name = "Abund") %>%
  filter(Abund > 0) %>%
  left_join(samples_df, by = "Sample_ID") %>%
  group_by(Region, Station_ID, Fraction) %>%
  dplyr::summarize(N_prot = n(), .groups = "drop") %>%
  tidyr::spread(Fraction, N_prot) %>%
  dplyr::rename("Cells_Total" = "Cells", "BEVs_Total" = "BEVs") |> 
  left_join(prot_per_sample, by =c("Station_ID", "Region")) %>% 
  left_join(viral_prot_per_sample, by =c("Station_ID", "Region")) %>% 
  mutate(Cells_Bac_prop = round(Cells_Bac/Cells_Total,2),
         BEVs_Bac_prop = round(BEVs_Bac/BEVs_Total,2),
         Cells_Vir_prop = round(Cells_Vir/Cells_Total,2),
         BEVs_Vir_prop = round(BEVs_Vir/BEVs_Total,2)) 

############################
#Calculate NSAF value and plot cpomposition per sample
############################
# counts of proteins present (Abundance>0) per Station_ID × Fraction
counts_df <- protein_abund |>
  tibble::rownames_to_column(var = "gene_callers_id") |>
  reshape2::melt(variable = "Sample_ID", value.name = "Abundance") |>
  left_join(samples_df, by = "Sample_ID") |>
  left_join(protein_taxonomy, by = "gene_callers_id") |>
  mutate(Domain = case_when(gene_callers_id %in% vir_gcids ~ "Viruses",
                            TRUE ~ Domain),
         present = !is.na(Abundance) & Abundance > 0) |>
  group_by(Station_ID, Fraction) |>
  summarize(
    total_proteins = sum(present, na.rm = TRUE),
    virus_proteins = sum(present & Domain == "Viruses", na.rm = TRUE),
    .groups = "drop"
  )

# total NSAF per Station_ID × Fraction (to position labels)
nsaf_mat <- add_nsaf_protein_abund(protein_abund, protein_metadata)

nsaf_totals <- nsaf_mat |>
  tibble::rownames_to_column(var = "gene_callers_id") |>
  reshape2::melt(variable = "Sample_ID", value.name = "NSAF") |>
  left_join(samples_df, by = "Sample_ID") |>
  group_by(Station_ID, Fraction) |>
  summarize(total_NSAF = sum(NSAF, na.rm = TRUE), .groups = "drop")

# combine and create label and y position (small offset above column)
labels <- nsaf_totals |>
  left_join(counts_df, by = c("Station_ID", "Fraction")) |>
  mutate(
    ratio_label = paste0(virus_proteins, "/", total_proteins),
    y = ifelse(is.na(total_NSAF), 0, total_NSAF * 1.03)
  )

# original plotting pipeline + labels
nsaf_mat |>
  tibble::rownames_to_column(var = "gene_callers_id") |> reshape2::melt(variable = "Sample_ID") |>
  left_join(protein_taxonomy, by = "gene_callers_id") |>
  left_join(samples_df, by = "Sample_ID") |>
  mutate(Domain = case_when(gene_callers_id %in% vir_gcids ~ "Viruses",
                            TRUE ~ Domain)) |>
  group_by(Station_ID, Fraction, Domain) |>
  summarize(Total = sum(value, na.rm = TRUE), .groups = "drop") |>
  mutate(Domain = factor(Domain, levels = c("Archaea","Bacteria", "Viruses", "Unassigned"))) |>
  ggplot(aes(x = Station_ID, y = Total, fill = Domain)) +
  geom_col() +
  facet_grid(~Fraction, scales = "free_y") + theme_bw() +
  labs(y = "Total Normalized Spectral Abundance Factor (NSAF)", x = "Station ID") +
  scale_fill_manual(values = c("Bacteria" = "#1f78b4",  "Archaea" = "#e31a1c", "Viruses" = "#ff7f00", "Unassigned" = "gray50")) +
  theme(legend.title = element_blank()) +
  geom_text(
    data = labels,
    mapping = aes(x = Station_ID, y = y, label = ratio_label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3
  )+
  theme_EF+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#save the plot
ggsave("./Figures/Fig_S4_metaP_tax_comp.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

############################
#Exclude non-bacterial proteins and export filtered abundance tables for downstream analyses
############################
protein_abund_filt<- protein_abund |>  as.data.frame() |> rownames_to_column("gene_callers_id") |>
  filter(!gene_callers_id %in% uncl_ids) 

protein_abund_filt |> 
  write_tsv("data/derived/protein_abundance_filt.tsv", col_names = TRUE)

# log2 transformation — preserve protein IDs as rownames
protein_abund.log2 <- protein_abund_filt |>
  tibble::column_to_rownames("gene_callers_id") |>
  as.matrix() |>
  log2()

# verify
cat(sprintf(
  "log2 matrix: %d proteins x %d samples\n  NAs: %d / %d (%.1f%%)\n  Range: %.2f - %.2f\n",
  nrow(protein_abund.log2), ncol(protein_abund.log2),
  sum(is.na(protein_abund.log2)),
  prod(dim(protein_abund.log2)),
  100 * mean(is.na(protein_abund.log2)),
  min(protein_abund.log2, na.rm = TRUE),
  max(protein_abund.log2, na.rm = TRUE)
))

protein_abund.log2 |> as.data.frame() |> rownames_to_column("gene_callers_id") |>
  write_tsv("data/derived/protein_abundance_filt_log2.tsv", col_names = TRUE)

#median normalization of the data
prot.dat.log2_norm <- DEqMS::equalMedianNormalization(protein_abund.log2)
row.names(prot.dat.log2_norm) <- rownames(protein_abund.log2)

prot.dat.log2_norm |> as.data.frame() |> rownames_to_column("gene_callers_id") |>
  write_tsv("data/derived/protein_abundance_filt_log2_norm.tsv", col_names = TRUE)


 


############################
#Export raw abudnance table with annotations
############################
export_prot_abund <- protein_abund_filt |>
  left_join(protein_metadata,    by = "gene_callers_id") |>
  left_join(
    protein_annotations |>
      select(gene_callers_id, InterPro_acc, InterPro_ann,
             NCBIfam_acc, NCBIfam_ann,
             Pfam_acc, Pfam_ann,
             Localization_final),
    by = "gene_callers_id"
  ) |>
  left_join(protein_taxonomy, by = "gene_callers_id") |>
  # sort: Domain → Phylum → Order → descending mean abundance
  mutate(
    mean_abund = rowMeans(
      across(all_of(samples_df$Sample_ID)),
      na.rm = TRUE
    )
  ) |>
  arrange(Domain, Phylum, Order, desc(mean_abund)) |>
  select(-mean_abund)

# ── build styled workbook ──────────────────────────────────────────────────
wb2 <- openxlsx::createWorkbook()

# styles — same as Fig_S2 script
header_style <- openxlsx::createStyle(
  fontColour     = "#FFFFFF",
  fgFill         = "#2E75B6",
  halign         = "CENTER",
  textDecoration = "bold",
  border         = "Bottom"
)
row_style <- openxlsx::createStyle(fgFill = "#F2F7FB")

# ── sheet 1: full abundance table ─────────────────────────────────────────
openxlsx::addWorksheet(wb2, sheetName = "Protein_abundance")
openxlsx::writeData(wb2,
                    sheet       = "Protein_abundance",
                    x           = export_prot_abund,
                    startRow    = 1,
                    startCol    = 1,
                    headerStyle = header_style)

# alternating row fill
even_rows <- seq(3, nrow(export_prot_abund) + 1, by = 2)
openxlsx::addStyle(wb2,
                   sheet      = "Protein_abundance",
                   style      = row_style,
                   rows       = even_rows,
                   cols       = seq_len(ncol(export_prot_abund)),
                   gridExpand = TRUE)

openxlsx::freezePane(wb2,  sheet = "Protein_abundance", firstRow = TRUE)
openxlsx::setColWidths(wb2, sheet = "Protein_abundance",
                       cols   = seq_len(ncol(export_prot_abund)),
                       widths = "auto")


# ── save ───────────────────────────────────────────────────────────────────
openxlsx::saveWorkbook(wb2,
                       file      = "Tables/Table_S2-Protein_abundance.xlsx",
                       overwrite = TRUE)

cat(sprintf(
  "Saved Table_S2:\n  Sheet 1 'Protein_abundance' : %d proteins x %d columns\n  ",
  nrow(export_prot_abund),
  ncol(export_prot_abund)
))



############################
#print session info and clean the workspace
############################  
sessionInfo()
rm(list = ls())
gc()
  
  
  
  
  
  
  
