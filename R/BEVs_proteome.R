require(tidyr)
require(tibble)
require(dplyr)
require(vegan)
require(stringr)
require(broom)
require(DEqMS)
require(cowplot)

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
run_source_if_missing("data/derived/protein_abundance_filt.tsv", "R/proteomic_overview.R")
protein_abund_filt<- read.table("data/derived/protein_abundance_filt.tsv", header = TRUE)|> 
  tibble::column_to_rownames("gene_callers_id")

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


############################
#Protein overlap between fractions
############################
prot_overlap_ls <- lapply(paste0(samples_df$Station_ID,"_"), function(s){
    st<- protein_abund_filt %>% rownames_to_column("gene_callers_id") |> 
      dplyr::select(gene_callers_id, starts_with(match=s)) %>% 
      reshape2::melt(value.name = "Abund") %>% 
      mutate(Fraction = case_when(grepl("Cells", variable)~ "Cells", grepl("EVs", variable)~ "EVs")) %>% 
      filter(Abund>0)
    
    overlap <- VennDiagram::calculate.overlap(list("Cells"= st %>% filter(Fraction =="Cells") %>% dplyr::select(gene_callers_id) %>%  pull(), 
                                                   "EVs"=st %>% filter(Fraction =="EVs") %>% dplyr::select(gene_callers_id) %>%  pull()))
    
    return(data.frame(Station_ID = s, 
                      Cells_prot = length(overlap$a1),
                      EVs_prot = length(overlap$a2),
                      Shared_prot = length(overlap$a3)))
  })
  
prot_overlap<- bind_rows(prot_overlap_ls) %>% 
    mutate(EVs_shared_prop = round(Shared_prot/EVs_prot,2)) %>% unique()
  
  
prot_overlap %>% 
    filter(Shared_prot>0) %>% 
    dplyr::summarize(Min=min(EVs_shared_prop),
              Max=max(EVs_shared_prop),
              Mean = mean(EVs_shared_prop),
              SE = se(EVs_shared_prop))


############################
#dissimilarity of Cells vs. BEV proteomes
############################
#count how many detections were for each protein
EV_sample_IDs<- samples_df %>% filter(Fraction =="BEVs") %>% pull(Sample_ID)
Cell_sample_IDs<- samples_df %>% filter(Fraction =="Cells") %>% pull(Sample_ID)

BEVs_gcids<- protein_abund_filt[,EV_sample_IDs] |> as.data.frame() |> rownames_to_column("gene_callers_id") |>
                reshape2::melt(id.vars = "gene_callers_id",
                        variable.name = "Sample_ID",
                        value.name = "Abund") %>% filter(Abund > 0) %>% pull(gene_callers_id) %>% unique()

Cells_gcids<- protein_abund_filt[,Cell_sample_IDs] |> as.data.frame() |> rownames_to_column("gene_callers_id") |>
                reshape2::melt(id.vars = "gene_callers_id",
                        variable.name = "Sample_ID",
                        value.name = "Abund") %>% filter(Abund > 0) %>% pull(gene_callers_id) %>% unique()

protein_abund_filt_no_na<- protein_abund_filt[intersect(BEVs_gcids, Cells_gcids),]
protein_abund_filt_no_na[is.na(protein_abund_filt_no_na)] <- 0
protein_abund_filt_hell <- vegan::decostand(protein_abund_filt_no_na, method = "hellinger")

#test differences between fraction
protein_dist_matrix <- vegdist(t(protein_abund_filt_hell), method = "euclidean",na.rm = TRUE)

adonis2(protein_dist_matrix ~ Fraction, samples_df,permutations=999)  

# post-hoc pairwise PERMANOVA: which pairs of regions are significantly different?
# Bonferroni correction applied across the 6 pairwise comparisons
pairwiseAdonis::pairwise.adonis(protein_dist_matrix,
                                factors = samples_df %>% mutate(factor=paste(Region, Fraction, sep="_")) %>% pull(factor),
                                p.adjust.m = 'bonferroni', perm = 999)

# betadisper: test homogeneity of within-group dispersion
# Significant result means groups differ in spread (not just centroid position)
# — this should be interpreted alongside adonis2 results
bd <- betadisper(protein_dist_matrix, samples_df %>% mutate(factor=paste(Region, Fraction, sep="_")) %>% pull(factor))
permutest(bd, permutations = 999)

# visualise spread
plot(bd, main = "Distance to centroid by Region")     # PCoA with centroids
boxplot(bd, xlab = "Region", ylab = "Distance to centroid")  # dispersion boxplot

############################
#Protein enrichment analysis between fractions in each region
############################
#carry out enrichment tests based on regions
enrichment_tests_list <- lapply(c("WEST","GYRE", "TRAN"), function(x) {
  
  samples_meta_sub<- samples_df %>% 
    mutate(Enr.group = factor(Region, levels =c("WEST","GYRE", "TRAN"))) %>% 
    filter(Enr.group ==x)
  
  #prepare experiment design matrix
  cond <-factor(samples_meta_sub$Fraction,
                levels = c("BEVs","Cells"))
  design <- model.matrix(~0+cond) # 0 means no intercept for the linear model
  colnames(design) <- gsub("cond","",colnames(design))
  contrast <-  makeContrasts(contrasts="BEVs-Cells",levels=design)
  
  #count how many detections were for each protein
  EV_sample_IDs<- samples_meta_sub %>% filter(Fraction =="BEVs") %>% pull(Sample_ID)
  Cell_sample_IDs<- samples_meta_sub %>% filter(Fraction =="Cells") %>% pull(Sample_ID)
  
  prot.dat.log2_norm.filter <- protein_abund_log2_norm[,c(EV_sample_IDs,Cell_sample_IDs)]
  
  # Filter proteins that were observed in at least two samples in each fraction
  prot.dat.log2_norm.filter <- prot.dat.log2_norm.filter[rowSums(!is.na(prot.dat.log2_norm.filter[, c(EV_sample_IDs)]))>=2 &
                                                           rowSums(!is.na(prot.dat.log2_norm.filter[, c(Cell_sample_IDs)]))>=2,]
  
  #run linear model
  fit1<- lmFit(prot.dat.log2_norm.filter, design)
  fit2 <- contrasts.fit(fit1,contrasts = contrast)
  fit3<- eBayes(fit2)
  


  #correct bias of variance estimate based on number of PSMs per protein
  # Reorder rows to match coefficients’ row order
  total_PSMs <- protein_metadata %>%
            filter(gene_callers_id %in% rownames(fit3$coefficients)) %>%
            mutate(gene_callers_id = factor(gene_callers_id, levels = rownames(fit3$coefficients))) %>%
            arrange(gene_callers_id) |> pull(PSMs)
  
  fit3$count <- total_PSMs

  fit4 = spectraCounteBayes(fit3)
  
  #results
  DEqMS.results <- outputResult(fit4,coef_col = 1) %>% 
                 dplyr::rename("gene_callers_id"="gene") %>% 
                 mutate(Enr.frac = case_when(logFC>1 & sca.adj.pval<0.1  
                                             ~ "BEVs",
                                             -1> logFC & sca.adj.pval<0.1  
                                             ~ "Cells", TRUE ~ "Not.enr"),
                        log.sca.pval = -log10(sca.P.Value),
                        Enr.reg =x) 
  
  return(DEqMS.results)
})

############################
#enrichment results 
############################
DEqMS_results_frac<- bind_rows(enrichment_tests_list) |>   filter(Enr.frac!="Not.enr")

DEqMS_results_not.enr<- bind_rows(enrichment_tests_list) |>   filter(Enr.frac=="Not.enr")


#total of different proteins
DEqMS_results_frac %>%  
  select(Enr.frac, gene_callers_id) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())

#total of proteins epr region
DEqMS_results_totals<- DEqMS_results_frac %>% 
  select(Enr.reg, Enr.frac, gene_callers_id) %>% 
  group_by(Enr.reg, Enr.frac) %>% 
  summarize(Total_p = n())


write.table(DEqMS_results_frac, "data/derived/DEqMS_results_fractions.txt", col.names =T, row.names = F, quote = F)

############################
# localization of the enriched proteins
############################
#plot localization of the enriched proteins
DEqMS_results_loc<- DEqMS_results_frac %>% 
  left_join(protein_annotations,  by="gene_callers_id") %>% 
  mutate(Localization_final=case_when(Localization_final=="Cell wall & surface"~"Extracellular",  TRUE~ Localization_final)) %>% 
  select(Enr.frac, gene_callers_id, Localization_final) %>% unique() |>
  group_by(Enr.frac, Localization_final) %>% 
  dplyr::summarize(N_p = n()) 
  
#total of proteins per fraction region
DEqMS_results_frac_totals<- DEqMS_results_frac %>% 
  select(Enr.frac, gene_callers_id) %>% unique() |>
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())

DEqMS_results_loc_prop.p<-DEqMS_results_loc %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(DEqMS_results_frac_totals, by = c("Enr.frac")) %>% 
  mutate(Prop=100*N_p/Total_p,
         Localization_final=factor(Localization_final, levels=c("Cytoplasmic","Cytoplasmic Membrane",
                                          "Periplasmic", "Outer Membrane",
                                          "Extracellular","Unassigned"))) %>%  
  ggplot(aes(x=Enr.frac, y= Prop, fill = Localization_final))+
  geom_col()+
  #geom_text(aes(label=Total_p), y=100, vjust=0.1, size = 5)+
  scale_fill_manual(values = c("#e5c185","#c7522a", "#74a892", "#008585", "#fbf2c4","red4"))+
  #facet_grid(cols=vars(Enr.reg),scales="free_x",space="free_x",switch="x") +
  ylab("Proportion of enriched proteins (%)")+
  theme_EF+
  guides(fill = guide_legend(title="Sub-cellular localization", ncol = 2))+
  theme(legend.position = "none",
        axis.title.x = element_blank())

############################
# Statistical test: differences in localization proportions between fractions
############################
loc_levels <- c("Cytoplasmic","Cytoplasmic Membrane",
                "Periplasmic", "Outer Membrane",
                "Extracellular","Unassigned")


# ── 1. enriched protein lists ─────────────────────────────────────────────────
# BEVs-enriched: logFC > 1, sca.adj.pval < 0.1
# Cells-enriched: logFC < -1, sca.adj.pval < 0.1

enr_BEVs  <- DEqMS_results_frac |> 
  filter(Enr.frac == "BEVs") |>
  select(Enr.reg, gene_callers_id) |> unique() |> 
  mutate(pass=TRUE, Enr.frac="BEVs")

enr_Cells <- DEqMS_results_frac |> 
  filter(Enr.frac == "Cells") |>
  select(Enr.reg, gene_callers_id) |> unique() |> 
  mutate(pass=TRUE, Enr.frac="Cells")

# ── 2. per-station detection matrix ──────────────────────────────────────────
# long format: gene_callers_id x Sample_ID → detected (T/F)
# then join sample metadata to get Station_ID and Fraction

abund_BEVs <- protein_abund_log2_norm |> as.data.frame() |>
  tibble::rownames_to_column("gene_callers_id") |>
  pivot_longer(
    cols      = -gene_callers_id,
    names_to  = "Sample_ID",
    values_to = "log2_abund"
  ) |>
    filter(!is.na(log2_abund)) |>            # detected = non-NA
  left_join(
    samples_df |> select(Sample_ID, Station_ID, Fraction, Region),
    by = "Sample_ID"
  ) |>
  left_join(enr_BEVs, by = c("gene_callers_id", "Fraction" = "Enr.frac",  "Region" = "Enr.reg")) |> 
  filter(pass==TRUE) |> 
   select(-pass)
  

abund_cells <- protein_abund_log2_norm |> as.data.frame() |>
  tibble::rownames_to_column("gene_callers_id") |>
  pivot_longer(
    cols      = -gene_callers_id,
    names_to  = "Sample_ID",
    values_to = "log2_abund"
  ) |>
    filter(!is.na(log2_abund)) |>            # detected = non-NA
  left_join(
    samples_df |> select(Sample_ID, Station_ID, Fraction, Region),
    by = "Sample_ID"
  ) |>
  left_join(enr_Cells, by = c("gene_callers_id", "Fraction" = "Enr.frac",  "Region" = "Enr.reg")) |> 
  filter(pass==TRUE) |> 
   select(-pass)
  

abund_long <- bind_rows(abund_BEVs, abund_cells)
 

# ── 3. join localization ──────────────────────────────────────────────────────
abund_long_loc <- abund_long |>
  left_join(
    protein_annotations |>
      select(gene_callers_id, Localization_final) |>
      mutate(
        Localization_final = case_when(
          Localization_final == "Cell wall & surface" ~ "Extracellular",
          Localization_final == "" | is.na(Localization_final) ~ "Unassigned",
          TRUE ~ Localization_final
        ),
        Localization_final = factor(
          Localization_final,
          levels = c("Cytoplasmic", "Cytoplasmic Membrane",
                     "Periplasmic", "Outer Membrane",
                     "Extracellular", "Unassigned")
        )
      ),
    by = "gene_callers_id"
  )

# ── 4. count proteins per station × fraction × localization ──────────────────
loc_per_station <- abund_long_loc |>
  # one row per unique protein per station per fraction
  # (protein may be detected in multiple samples at same station)
  distinct(gene_callers_id, Station_ID, Fraction,
           Region, Localization_final) |>
  dplyr::count(Station_ID, Fraction, Region,
               Localization_final, name = "N_prot") |>
  group_by(Station_ID, Fraction) |>
  mutate(
    Total    = sum(N_prot),
    Prop     = N_prot / Total,
    Prop_pct = 100 * Prop
  ) |>
  ungroup() |>
  mutate(
    Station_ID = factor(Station_ID,
                        levels = c("SO289_44","SO289_43","SO289_41","SO289_39",
                                   "SO289_37","SO289_34","SO289_33","SO289_32",
                                   "SO289_30","SO289_27","SO289_23","SO289_20",
                                   "SO289_17","SO289_16","SO289_13","SO289_12",
                                   "SO289_9","SO289_6","SO289_3","SO289_1")),
    Region     = factor(Region, levels = c("WEST","GYRE","TRAN","UP"))
  )

# ── 5. statistical test: differences between fractions ───────────────────────
paired_results <- purrr::map_dfr(loc_levels, function(loc) {

  # aggregate across physical fractions BEFORE pivoting
  # (avoids duplicate Station_ID × Fraction → list-column error)
  df_wide <- loc_per_station |>
    filter(Localization_final == loc) |>
    group_by(Station_ID, Fraction, Region) |>
    dplyr::summarise(
      N_prot = sum(N_prot),
      Total  = sum(Total),
      .groups = "drop"
    ) |>
    mutate(Prop = N_prot / Total) |>
    pivot_wider(
      names_from  = Fraction,
      values_from = c(Prop, N_prot, Total),
      values_fill = 0
    )

  # need both fractions and at least 3 stations
  if (!all(c("Prop_BEVs","Prop_Cells") %in% names(df_wide)) |
      nrow(df_wide) < 3) {
    return(tibble(
      Localization = loc, n_stations = nrow(df_wide),
      mean_BEVs = NA, mean_Cells = NA, mean_diff = NA,
      median_diff = NA, CI_lo = NA, CI_hi = NA,
      W = NA, p_wilcox = NA,
      t = NA, p_ttest = NA,
      p.adj = NA, sig = NA, direction = NA
    ))
  }

  # Wilcoxon signed-rank (non-parametric, paired)
  w_test <- wilcox.test(
    df_wide$Prop_BEVs,
    df_wide$Prop_Cells,
    paired     = TRUE,
    exact      = FALSE,
    conf.int   = TRUE,
    conf.level = 0.95
  )

  # paired t-test (parametric — reported alongside for comparison)
  t_test <- t.test(
    df_wide$Prop_BEVs,
    df_wide$Prop_Cells,
    paired = TRUE
  )

  tibble(
    Localization = loc,
    n_stations   = nrow(df_wide),
    mean_BEVs    = round(mean(df_wide$Prop_BEVs,  na.rm = TRUE), 4),
    mean_Cells   = round(mean(df_wide$Prop_Cells, na.rm = TRUE), 4),
    mean_diff    = round(mean(df_wide$Prop_BEVs - df_wide$Prop_Cells,
                              na.rm = TRUE), 4),
    median_diff  = round(w_test$estimate,      4),   # Hodges-Lehmann estimator
    CI_lo        = round(w_test$conf.int[1],   4),
    CI_hi        = round(w_test$conf.int[2],   4),
    W            = w_test$statistic,
    p_wilcox     = w_test$p.value,
    t            = round(t_test$statistic,     3),
    p_ttest      = t_test$p.value,
    p.adj        = NA_real_
  )
})

# FDR correction across 6 localization categories
paired_results <- paired_results |>
  mutate(
    p.adj     = p.adjust(p_wilcox, method = "BH"),
    sig       = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.1   ~ ".",
      TRUE          ~ "ns"
    ),
    direction = case_when(
      mean_diff > 0 & p.adj < 0.05 ~ "BEVs enriched",
      mean_diff < 0 & p.adj < 0.05 ~ "Cells enriched",
      TRUE                          ~ "ns"
    ),
    Localization_final = factor(Localization, levels = loc_levels)
  )

cat("\n── Paired Wilcoxon results ──────────────────────────────────────────\n")
paired_results |>
  select(Localization, n_stations, mean_BEVs, mean_Cells,
         mean_diff, CI_lo, CI_hi, W, p_wilcox, p.adj, sig) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  print()

# ── 6b. boxplot — replace glmm_results with paired_results ───────────────────
loc_boxplot <- loc_per_station |>
  left_join(
    paired_results |> select(Localization_final, sig),
    by = "Localization_final"
  ) |>
  ggplot(aes(x = Fraction, y = Prop_pct, fill =Region)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  #geom_jitter(aes(shape = Region), width = 0.15, size = 2, alpha = 0.8) +
  geom_text(
    data        = paired_results,
    aes(x = 1.5, y = Inf, label = sig),
    inherit.aes = FALSE,
    vjust = 1.5, size = 5, colour = "black"
  ) +
  #scale_fill_manual(values   = c("BEVs"  = "#ff0000",
  #                               "Cells" = "#0000ff")) +
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00")) +
  facet_wrap(~ Localization_final, scales = "free_y", nrow = 2) +
  labs(x        = NULL,
       y        = "Enriched proteins (%)",
       fill     = "Fraction", colour = "Region", shape = "Region",
       #title    = "Localization proportion: BEVs vs Cells enriched proteins",
       #subtitle = "Paired Wilcoxon signed-rank | * p.adj<0.05 | ** p.adj<0.01 | *** p.adj<0.001"
      ) +
  theme_EF +
  theme(legend.position = "none",
        strip.text      = element_text(face = "bold"))


cowplot::plot_grid(DEqMS_results_loc_prop.p, loc_boxplot, ncol = 2, align = "hv")

#save the plot
ggsave("./Figures/cobined_cells_BEVs_enr_prot_local.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

############################
# taxonomy of the enriched proteins
############################
#whoich are the most enriched taxa in the enriched proteins?
DEqMS_results_order_prop<- DEqMS_results_frac %>% filter(Enr.frac!="Not.enr") |> 
  left_join(protein_taxonomy, by = "gene_callers_id") |> 
  select(Enr.reg, Enr.frac, gene_callers_id, Phylum, Class, Order) %>% unique() |>
  group_by(Enr.reg, Enr.frac, Phylum, Class, Order) %>% 
  dplyr::summarize(N_p = n()) |> 
  left_join(DEqMS_results_totals, by = c("Enr.reg","Enr.frac") )%>% 
  mutate(Prop = 100 * N_p / Total_p) |> 
  group_by(Enr.frac, Order) |> 
  summarise(min_prop = min(Prop), max_prop=max(Prop)) 



#plot taxa of the enriched proteins
DEqMS_results_order_tax<- DEqMS_results_frac %>% filter(Enr.frac!="Not.enr") |> 
  left_join(protein_taxonomy, by = "gene_callers_id") |> 
  select(Enr.reg, Enr.frac, gene_callers_id, Phylum, Class, Order) %>% unique() |>
  group_by(Enr.reg, Enr.frac, Phylum, Class, Order) %>% 
  dplyr::summarize(N_p = n())




# ── build sorted order label map: Class A→Z, then Order A→Z within class ─────
# Label format: "Order (X)" where X = first letter of Class
# e.g. Flavobacteriales in Bacteroidia  → "Flavobacteriales (B)"
#      Alteromonadales in Gammaproteobacteria → "Alteromonadales (G)"
# Flavobacteriales (B) sorts before Alteromonadales (G) because B < G
order_label_map <- DEqMS_results_order_tax |> 
  filter(!is.na(Order), !is.na(Class), Enr.frac!="Not.enr") |>
  mutate(
    Order        = gsub("^ ",         "", Order),
    Class        = gsub("^ ",         "", Class),
    Order_clean  = gsub("Candidatus ", "", Order),
    Class_clean  = gsub("Candidatus ", "", Class),
    Class_letter = toupper(substr(Class_clean, 1, 1)),
    Label        = paste0(Order_clean, " (", Class_letter, ")")
  ) |>
  distinct(Order_clean, Class_clean, Class_letter, Label) |>
  arrange(Class_clean, Order_clean)

# factor levels: sorted named orders first, fallback categories last
ordered_labels <- c(unique(order_label_map$Label), "Other taxa < 1%", "unclassified")

DEqMS_results_order_tax_prop<- DEqMS_results_order_tax %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(DEqMS_results_totals, by = c("Enr.reg","Enr.frac") )%>% 
  mutate(
    Prop         = 100 * N_p / Total_p,
    Order        = gsub("^ ", "", Order),
    Class        = gsub("^ ", "", Class),
    Order_clean  = gsub("Candidatus ", "", Order),
    Class_clean  = gsub("Candidatus ", "", Class),
    Class_letter = toupper(substr(Class_clean, 1, 1)),
    Taxa = case_when(
      Prop < 1                    ~ "Other taxa < 1%",
      is.na(Order) & is.na(Class) ~ "unclassified",
      is.na(Order)                ~ "unclassified",
      TRUE ~ paste0(Order_clean, " (", Class_letter, ")")
    ),
    Taxa     = factor(Taxa, levels = ordered_labels),
    Enr.reg     = factor(Enr.reg, levels = c("WEST","GYRE","TRAN","UP"))
  ) |> droplevels()

taxa_cols<- setNames(tol21rainbow, unique(DEqMS_results_order_tax_prop$Taxa)[seq_along(tol21rainbow)])


DEqMS_results_order_tax_prop.p<- DEqMS_results_order_tax_prop%>%
  ggplot(aes(x = Enr.frac, y = Prop, fill = Taxa)) +
  geom_col() +
  scale_fill_manual(
    values   = taxa_cols, na.value = "grey70"
  ) +
  ylab("Enriched proteins (%)") +
  facet_grid(cols=vars(Enr.reg),scales="free_x",space="free_x",switch="x") +
  theme_EF +
  guides(fill = guide_legend(nrow = 4)) +
  theme(legend.position = "none",
        axis.title.x    = element_blank())


DEqMS_results_order_tax_no_reg<- DEqMS_results_frac |> filter(Enr.frac!="Not.enr") |> 
  left_join(protein_taxonomy, by = "gene_callers_id") |> 
  select(Enr.frac, gene_callers_id, Phylum, Class, Order) |> unique() |>
  group_by(Enr.frac, Phylum, Class, Order) |> 
  dplyr::summarize(N_p = n()) 

DEqMS_results_order_totals<- bind_rows(enrichment_tests_list) %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") |> 
  select(gene_callers_id, Order) |> unique() |>
  group_by(Order) %>% dplyr::summarize(Total_p = n()) 

#plot ratio of enriched proteins
DEqMS_results_order_tax_ratio.p<- DEqMS_results_order_tax_no_reg %>% 
    left_join(DEqMS_results_order_totals, by = c("Order")) %>% 
  filter(Total_p>=10, !is.na(Order)) |> 
  mutate(Prop         = 100 * N_p / Total_p) |> 
  select(-c(N_p, Total_p)) |>
  tidyr::spread(Enr.frac, Prop) |> filter(!is.na(BEVs)&!is.na(Cells)&!is.na(Phylum)) %>%  
  mutate(
    Ratio         = BEVs/Cells,
    Order        = gsub("^ ", "", Order),
    Class        = gsub("^ ", "", Class),
    Order_clean  = gsub("Candidatus ", "", Order),
    Class_clean  = gsub("Candidatus ", "", Class),
    Class_letter = toupper(substr(Class_clean, 1, 1)),
    Taxa = case_when(
      BEVs <= 1 &        Cells <= 1  ~ "Other taxa < 1%",
      is.na(Order) & is.na(Class) ~ "unclassified",
      is.na(Order)                ~ "unclassified",
      TRUE ~ paste0(Order_clean, " (", Class_letter, ")")
    ),
    Taxa     = factor(Taxa, levels = ordered_labels)
  ) %>% 
  ggplot(aes(x=Cells, y= BEVs, colour= Taxa, label = Taxa))+
  geom_point(size=5)+
  ggrepel::geom_label_repel(color = "red", nudge_x = .5, nudge_y = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.8)+
  scale_colour_manual(values = taxa_cols)+
  labs(x="Enriched proteins in cells (%)", y = "Enriched proteins in BEVs (%)")+
  theme_EF+
  theme(legend.position = "none")

cowplot::plot_grid(DEqMS_results_order_tax_prop.p, DEqMS_results_order_tax_ratio.p, ncol = 2, align = "hv")

#save the plot
ggsave("./Figures/cobined_cells_BEVs_enr_prot_taxa.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

cowplot::plot_grid(DEqMS_results_order_tax_prop.p, DEqMS_results_order_tax_ratio.p, ncol = 2, align = "hv")

#save the plot
ggsave("./Figures/cobined_cells_BEVs_enr_prot_taxa_legend.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 180, 
       scale = 2,
       dpi = 300)





############################
# Export results
############################
# one sheet per fraction (BEVs / Cells), sorted by logFC descending
# styled: bold blue header | alternating row fill | frozen top row | auto widths

export_frac_df <- DEqMS_results_frac |>
  filter(Enr.frac != "Not.enr") |>
  left_join(protein_annotations,
            by = "gene_callers_id",
            relationship = "many-to-many") |>
  unique() |>
  left_join(protein_taxonomy, by = "gene_callers_id") |>
  left_join(protein_metadata,  by = "gene_callers_id") |>
  select(any_of(c(
    "Enr.frac", "gene_callers_id",
    "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
    "sca.P.Value", "sca.adj.pval", "log.sca.pval",
    names(protein_annotations)[str_starts(names(protein_annotations), "InterPro_")],
    names(protein_annotations)[str_starts(names(protein_annotations), "NCBIfam_")],
    names(protein_annotations)[str_starts(names(protein_annotations), "Pfam_")],
    "Domain", "Phylum", "Class", "Order", "Family", "Genus",
    "PSMs", "Unique_pep", "Total_pep", "Length"
  ))) |>
  arrange(Enr.frac, desc(logFC))

# split into one data frame per fraction
frac_sheets     <- split(export_frac_df, export_frac_df$Enr.frac)
frac_sheet_order <- intersect(c("BEVs", "Cells"), names(frac_sheets))
frac_sheets      <- frac_sheets[frac_sheet_order]

# ── styled workbook ───────────────────────────────────────────────────────────
wb_frac <- openxlsx::createWorkbook()

header_style <- openxlsx::createStyle(
  fontColour     = "#FFFFFF",
  fgFill         = "#2E75B6",
  halign         = "CENTER",
  textDecoration = "bold",
  border         = "Bottom"
)
row_style <- openxlsx::createStyle(fgFill = "#F2F7FB")

purrr::walk(frac_sheet_order, function(frac) {
  df <- frac_sheets[[frac]]
  openxlsx::addWorksheet(wb_frac, sheetName = frac)
  openxlsx::writeData(wb_frac,
                      sheet       = frac,
                      x           = df,
                      startRow    = 1,
                      startCol    = 1,
                      headerStyle = header_style)
  # alternating row fill
  even_rows <- seq(3, nrow(df) + 1, by = 2)
  openxlsx::addStyle(wb_frac,
                     sheet      = frac,
                     style      = row_style,
                     rows       = even_rows,
                     cols       = seq_len(ncol(df)),
                     gridExpand = TRUE)
  openxlsx::freezePane(wb_frac,   sheet = frac, firstRow = TRUE)
  openxlsx::setColWidths(wb_frac, sheet = frac,
                         cols   = seq_len(ncol(df)),
                         widths = "auto")
})

openxlsx::saveWorkbook(wb_frac,
                       file      = "Tables/Table_S3_DEqMS_results_Fractions.xlsx",
                       overwrite = TRUE)

cat(sprintf(
  "Saved Table_S4 — %d sheets: %s\n  Rows per sheet: %s\n",
  length(frac_sheets),
  paste(frac_sheet_order, collapse = " | "),
  paste(sapply(frac_sheets, nrow), collapse = " | ")
))

  


#print session info and clean the workspace
sessionInfo()
rm(list = ls())
#gc()
