require(tidyr)
require(dplyr)
require(vegan)
require(stringr)
require(broom)
require(DEqMS)

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
run_source_if_missing("data/derived/protein_abundance_filt_log2_norm.tsv", "R/proteomic_overview.R")
protein_abund_log2_norm <- read.table("data/derived/protein_abundance_filt_log2_norm.tsv", header = TRUE)|> 
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

# -----------------------------------------------------------------------------
# 2. Dissimilarity between cellular proteomes
# -----------------------------------------------------------------------------
# Euclidean distance on log2-normalised abundances approximates the L2 norm
# in protein expression space. We test whether Region explains a significant
# proportion of the total variance (PERMANOVA) and whether within-group
# dispersions are homogeneous (betadisper — a key assumption of adonis2).

prot.dat.log2_norm_Cells<- protein_abund_log2_norm |> select(all_of(samples_df %>% filter(Fraction == "Cells") %>% pull(Sample_ID))) %>%
                            filter(rowSums(!is.na(.)) > 0) # keep proteins observed in at least one sample

# pairwise Euclidean distance matrix (samples x samples)
protein_dist_matrix <- vegdist(t(prot.dat.log2_norm_Cells), method = "euclidean", na.rm = TRUE)


# PERMANOVA: does Region explain proteome composition?
# R2 = fraction of total variance explained by Region; p-value from permutations
adonis2(protein_dist_matrix ~ Region, samples_df %>% filter(Fraction == "Cells"), permutations = 999)

# post-hoc pairwise PERMANOVA: which pairs of regions are significantly different?
# Bonferroni correction applied across the 6 pairwise comparisons
pairwiseAdonis::pairwise.adonis(protein_dist_matrix,
                                factors = samples_df %>% filter(Fraction == "Cells") %>% pull(Region),
                                p.adjust.m = 'bonferroni', perm = 999)

# betadisper: test homogeneity of within-group dispersion
# Significant result means groups differ in spread (not just centroid position)
# — this should be interpreted alongside adonis2 results
bd <- betadisper(protein_dist_matrix, samples_df %>% filter(Fraction == "Cells") %>% pull(Region))
permutest(bd, permutations = 999)

# visualise within-region spread
plot(bd, main = "Distance to centroid by Region")     # PCoA with centroids
boxplot(bd, xlab = "Region", ylab = "Distance to centroid")  # dispersion boxplot

# -----------------------------------------------------------------------------
# 4. NMDS ordination — Fig S2
# -----------------------------------------------------------------------------

prot_NMDS <- metaMDS(protein_dist_matrix,
                     k         = 4,          # number of ordination axes
                     maxit     = 999,
                     trymax    = 999,         # try many random starts to find global minimum
                     wascores  = FALSE,       # do not compute weighted average species scores
                     autotransform = FALSE,   # data already log2-normalised; skip auto-transform
                     tidy      = "sites")

# combine NMDS scores with sample metadata and cluster assignments for plotting
prot_NMDS.scores <- as.data.frame(scores(prot_NMDS)) %>%
  tibble::rownames_to_column(var = "Sample_ID") %>%
  left_join(samples_df, by = "Sample_ID", copy = TRUE) 

# plot NMDS coloured by Region
prot_NMDS.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2,
             colour = Region, label = Sample_ID)) +
  geom_point(size = 7, colour = "black") +   # black outline
  geom_point(size = 5) +                     # coloured fill
  scale_color_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  theme_EF +
  #geom_text(size = 3, vjust = -1.5) +
  guides(colour = guide_legend(title = "Region", ncol = 4)) +
  coord_fixed() +
  labs(title = sprintf("NMDS (euclidean), stress = %.3f", prot_NMDS$stress))+
  theme(legend.position = "bottom")

ggsave("./Figures/Fig_S2-Cell_prot_NMDS.pdf",
       plot = last_plot(),
       units = "mm", width = 90, height = 90, scale = 2, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Differential abundance between adjacent regions (DEqMS)
# -----------------------------------------------------------------------------
# Comparisons run pairwise between adjacent regions along the transect:
#   UP vs TRAN | TRAN vs GYRE | GYRE vs WEST
# Proteins are called enriched if |logFC| > 1 and sca.adj.pval < 0.1.

enrichment_tests_list <- lapply(c("UP_TRAN","TRAN_GYRE","GYRE_WEST"), function(x) {

  # subset metadata to the two regions being compared
  samples_meta_sub <- samples_df %>% filter(Fraction == "Cells") %>%
    mutate(Enr.group = factor(Region, levels = c("WEST","GYRE","TRAN","UP"))) %>%
    filter(Enr.group %in% c(str_split_1(x, "_")))

  # design matrix: one column per group, no intercept
  # each column indicates group membership (0/1) for each sample
  cond   <- factor(samples_meta_sub$Enr.group, levels = c(stringr::str_split_1(x, "_")))
  design <- model.matrix(~0 + cond)
  colnames(design) <- gsub("cond", "", colnames(design))

  # contrast: difference between the two regions (e.g. UP - TRAN)
  contrast <- makeContrasts(
    contrasts = paste(stringr::str_split_1(x, "_")[1], stringr::str_split_1(x, "_")[2], sep = "-"),
    levels = design
  )

  regA_sample_IDs <- samples_meta_sub %>% filter(Enr.group == str_split_1(x, "_")[1]) %>% pull(Sample_ID)
  regB_sample_IDs <- samples_meta_sub %>% filter(Enr.group == str_split_1(x, "_")[2]) %>% pull(Sample_ID)

  # apply same prevalence filter as step 2, then additionally require each
  # protein to be observed in > 1 sample per region (avoids single-point fits)
  prot.dat.log2_norm.filter <- protein_abund_log2_norm[, c(regA_sample_IDs, regB_sample_IDs), drop = FALSE]

  prot.dat.log2_norm.filter <- prot.dat.log2_norm.filter[
    rowSums(!is.na(prot.dat.log2_norm.filter[, regA_sample_IDs])) > 1 &
    rowSums(!is.na(prot.dat.log2_norm.filter[, regB_sample_IDs])) > 1, ]

  # step 1: fit linear model per protein (limma)
  fit1 <- lmFit(prot.dat.log2_norm.filter, design)
  # step 2: apply contrast (computes logFC for the comparison)
  fit2 <- contrasts.fit(fit1, contrasts = contrast)
  # step 3: empirical Bayes moderation of variance (shrinks per-protein variance
  #         towards a global prior — stabilises estimates for low-n comparisons)
  fit3 <- eBayes(fit2)

 # step 4: DEqMS correction — attach peptide counts and re-shrink variance

  # Reorder rows to match coefficients’ row order
  total_pep <- protein_metadata %>%
            filter(gene_callers_id %in% rownames(fit3$coefficients)) %>%
            mutate(gene_callers_id = factor(gene_callers_id, levels = rownames(fit3$coefficients))) %>%
            arrange(gene_callers_id) |> pull(Total_pep)
  

  fit3$count <- total_pep
  fit4 <- spectraCounteBayes(fit3)

  # extract results and classify enrichment direction
  DEqMS.results <- outputResult(fit4, coef_col = 1) %>%
    dplyr::rename("gene_callers_id" = "gene") %>%
    mutate(
      # Enr.reg: which region is enriched? |logFC| > 1 and sca.adj.pval < 0.1
      Enr.reg = case_when(
        logFC > 1  & sca.adj.pval < 0.1 ~ str_split_1(x, "_")[1],
        -1 > logFC & sca.adj.pval < 0.1 ~ str_split_1(x, "_")[2],
        TRUE ~ "Not.enr"),
      log.sca.pval = -log10(sca.P.Value),  # for volcano plotting
      Comp = x                              # comparison label e.g. "UP_TRAN"
    )
  return(DEqMS.results)
})

# combine all three pairwise comparisons into one dataframe
DEqMS_results_regions <- bind_rows(enrichment_tests_list) |> filter(Enr.reg != "Not.enr")


# summary: unique enriched proteins per region x Class x Order
# with total per region and proportion (one row per unique protein via distinct())
DEqMS_results_taxa_summary <- DEqMS_results_regions %>%
  left_join(protein_taxonomy, by = "gene_callers_id") %>%
  distinct(Enr.reg, Class, Order, gene_callers_id) %>%  unique() |> # one row per unique protein
  dplyr::count(Enr.reg, Class, Order, name = "n_proteins") %>%
  group_by(Enr.reg) %>%
  mutate(
    n_total_region = sum(n_proteins),             # total enriched proteins in that region
    prop = round(n_proteins / n_total_region, 3)  # proportion of region total
  ) %>%
  ungroup() %>%
  mutate(Enr.reg = factor(Enr.reg, levels = c("UP", "TRAN", "GYRE", "WEST"))) %>%
  arrange(Enr.reg, desc(n_proteins))


#plot functional differences between regions
DEqMS_results_regions_functions <- DEqMS_results_regions %>% 
  left_join(protein_annotations , by ="gene_callers_id") |> 
    left_join(protein_taxonomy %>% unique(), by ="gene_callers_id") |> 
    mutate(
    Category = categorize_single_loc(Function, Localization_final,
                                     patterns_bacteria_loc))


DEqMS_results_regions_totals <- DEqMS_results_regions_functions %>%
  group_by(Comp, Enr.reg) %>%
  summarize(Total_p = n(), .groups = "drop")

DEqMS_results_regions_CAT <- DEqMS_results_regions_functions |>
  group_by(Comp, Enr.reg, COG20cat_ann) |>
  summarise(
    count         = n(),
    log2fold_mean = mean(logFC),
    log2fold_se   = se(logFC),
    .groups       = "drop"
  ) |>
  left_join(DEqMS_results_regions_totals, by = c("Enr.reg", "Comp")) |>
  mutate(
    Prop    = count / Total_p,
    Enr.reg = factor(Enr.reg, levels = c("WEST","GYRE","TRAN", "UP"))
  )

DEqMS_results_regions_CAT |> 
  filter(COG20cat_ann %in% #sel_cats
    c(    DEqMS_results_regions_CAT |> 
      filter( count >= 3) |>
      pull(COG20cat_ann) |> unique()
  )) |>
  mutate(COG20cat_ann=gsub("_"," ",COG20cat_ann)) |> 
  mutate(COG20cat_ann=case_when(is.na(COG20cat_ann) ~ "Unclassified", TRUE ~ COG20cat_ann)) |> 
  ggplot(aes(y     = COG20cat_ann,
             x     = log2fold_mean,
             fill  = Enr.reg,
             label = count)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = log2fold_mean - log2fold_se,
                     xmax = log2fold_mean + log2fold_se),
                 height = 0.25, linewidth = 0.4, colour = "grey40") +
  geom_point(aes(size = Prop), shape = 21, alpha = 0.9) +
  scale_size_continuous(range  = c(4, 14),
                        labels = scales::percent_format(accuracy = 0.1)) +
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  geom_text(size = 3.5, #nudge_y = -0.45, nudge_x = -1, 
    colour = "black") +
  xlim(c(-6,6)) +
  facet_grid(. ~Comp, space = "free", switch = "y")+
  labs(
    x        = expression("Mean log"[2]*" fold change (BEVs / Cells)"),
    y        = NULL,
    fill     = "Region",
    #size     = "% of sar11bacterial\nenriched proteins",
    #title    = "sar11bacterial protein enrichment in BEVs vs Cells",
    #subtitle = "Bubble size = proportion of total enriched sar11bacterial proteins | n = protein count"
  ) +
  theme_EF +
  theme(
    legend.position    = "bottom",
    axis.text.x        = element_text(angle = 90, vjust = 0.5),
    axis.text.y        = element_text(size = 14),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
    strip.text         = element_text(face = "bold")
  )

ggsave("./Figures/Fig_3-Cell_prot_enr_region.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

# -----------------------------------------------------------------------------
# 7. Export results to Excel — one sheet per comparison
# -----------------------------------------------------------------------------
# Columns selected: comparison metadata | DEqMS statistics | taxonomy |
#                   functional annotations | peptide/PSM counts | sequence
# Sorted: Enr.reg (factor UP > TRAN > GYRE > WEST) then logFC descending

export_df <- DEqMS_results_regions |>
  filter(Enr.reg != "Not.enr") |>
  select(any_of(c(
    "Comp", "Enr.reg",
    names(DEqMS_results_regions),
    "Phylum", "Class", "Order", "Family", "Genus",
    names(protein_annotations)[str_starts(names(protein_annotations), "InterPro_")],
    names(protein_annotations)[str_starts(names(protein_annotations), "NCBIfam_")],
    names(protein_annotations)[str_starts(names(protein_annotations), "Pfam_")],
    "PSMs", "Unique_pep", "Total_pep", "Length", "AA_sequence"
  ))) |>
  mutate(
    Enr.reg = factor(Enr.reg, levels = c("UP","TRAN","GYRE","WEST"))
  ) |>
  arrange(Comp, Enr.reg, desc(logFC))

# split into one data frame per comparison
sheets <- split(export_df, export_df$Comp)

# order sheets along the transect
sheet_order <- intersect(c("UP_TRAN","TRAN_GYRE","GYRE_WEST"), names(sheets))
sheets <- sheets[sheet_order]

# build workbook — one styled sheet per comparison
wb <- openxlsx::createWorkbook()

# header style: bold + light blue fill
header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF",
  fgFill     = "#2E75B6",
  halign     = "CENTER",
  textDecoration = "bold",
  border     = "Bottom"
)

# alternating row style for readability
row_style <- openxlsx::createStyle(fgFill = "#F2F7FB")

purrr::walk(sheet_order, function(comp) {
  df <- sheets[[comp]]
  openxlsx::addWorksheet(wb, sheetName = comp)
  openxlsx::writeData(wb, sheet = comp, x = df,
                      startRow = 1, startCol = 1,
                      headerStyle = header_style)
  # alternate row fill for every other data row
  even_rows <- seq(3, nrow(df) + 1, by = 2)  # +1 for header offset
  openxlsx::addStyle(wb, sheet = comp,
                     style = row_style,
                     rows  = even_rows,
                     cols  = seq_len(ncol(df)),
                     gridExpand = TRUE)
  # freeze top row
  openxlsx::freezePane(wb, sheet = comp, firstRow = TRUE)
  # auto-fit column widths (capped at 40 characters)
  openxlsx::setColWidths(wb, sheet = comp,
                         cols   = seq_len(ncol(df)),
                         widths = "auto")
})

openxlsx::saveWorkbook(wb,
                       file      = "Tables/Table_S3_DEqMS_results_regions.xlsx",
                       overwrite = TRUE)

cat(sprintf(
  "Saved Table_S3 — %d sheets: %s\n  Rows per sheet: %s\n",
  length(sheets),
  paste(sheet_order, collapse = " | "),
  paste(sapply(sheets, nrow), collapse = " | ")
))

# print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
