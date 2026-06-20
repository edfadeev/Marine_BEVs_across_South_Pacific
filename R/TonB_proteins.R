#Explore operons with enriched tonB proteins in BEVs
require(tidyr)
require(tibble)
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

tonB_context<- read.table("data/derived/tonB_context.txt", header = TRUE, sep = "\t") |> 
  mutate(gene_callers_id = as.character(gene_callers_id))

############################
#Select only outer membrane and extracellular proteins for enrichment analysis
############################
nonCyto_proteins <- protein_annotations |> filter(Localization_final %in% c("Periplasmic", "Outer Membrane", "Extracellular") & reliability_score>=2 ) |> pull(gene_callers_id) 

protein_abund_log2_norm_nonCyto <- protein_abund_log2_norm[rownames(protein_abund_log2_norm) %in% nonCyto_proteins,]

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
  
  prot.dat.log2_norm.filter <- protein_abund_log2_norm_nonCyto[,c(EV_sample_IDs,Cell_sample_IDs)]
  
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
                 mutate(log.sca.pval = -log10(sca.P.Value),
                        Enr.reg =x) 
  
  return(DEqMS.results)
})

############################
#enrichment results 
############################
DEqMS_results_nonCyto<- bind_rows(enrichment_tests_list) |>  
                          mutate(Enr.frac = case_when(logFC>1& sca.adj.pval<0.1  
                                             ~ "BEVs",
                                             -1> logFC & sca.adj.pval<0.1  
                                             ~ "Cells", TRUE ~ "Not.enr"))
  
#total of different proteins
DEqMS_results_nonCyto %>%  filter(Enr.frac!="Not.enr") |> 
  select(Enr.frac, gene_callers_id) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())

#total of proteins epr region
DEqMS_results_nonCyto_totals<- DEqMS_results_nonCyto %>% 
  select(Enr.reg, Enr.frac, gene_callers_id) %>% 
  group_by(Enr.reg, Enr.frac) %>% 
  summarize(Total_p = n())

############################
# functions of enriched proteins
############################
DEqMS_results_nonCyto_annotations<- DEqMS_results_nonCyto %>% filter(Enr.frac!="Not.enr") |> 
  left_join(protein_annotations,  by="gene_callers_id") %>% 
   mutate(Category = categorize_single_loc(Function, Localization_final,
                                     patterns_bacteria_loc)) |> 
  left_join(protein_taxonomy %>% unique(), by ="gene_callers_id") |> 
  mutate(Category = case_when(grepl("solute", Pfam_ann, ignore.case=TRUE)|grepl("substrate", blastp_ann, ignore.case=TRUE)  ~ "Solute binding proteins (SBP)", 
                              grepl("peptidoglycan-associated lipoprotein", blastp_ann, ignore.case=TRUE)  ~ "Peptidoglycan / cell wall", 
                              grepl("Dcap", Pfam_ann, ignore.case=TRUE) ~ "Outer membrane porins (other)", 
                              NCBIfam_acc=="TIGR04057" ~"SusC (polysaccharide utilization loci)",
                              TRUE ~ Category))


DEqMS_results_nonCyto_CAT_df <- DEqMS_results_nonCyto_annotations |>
  group_by(Enr.reg, Enr.frac, Category) |>
  summarise(
    count         = n(),
    log2fold_mean = mean(logFC),
    log2fold_se   = se(logFC),
    .groups       = "drop"
  ) |>
  left_join(DEqMS_results_nonCyto_totals, by = c("Enr.reg","Enr.frac")) |>
  mutate(
    Prop    = count / Total_p,
    Enr.reg = factor(Enr.reg, levels = c("WEST","GYRE","TRAN"))
  )


FUN_categories <- enframe(bacteria_categories_list, name = "group", value = "category") %>%
  unnest_longer(category) %>%
  rename(category_label = category)

category_order_present <- FUN_categories[
  FUN_categories %in% unique(DEqMS_results_nonCyto_CAT_df$Category)
]
category_order_present <- c(
  setdiff(unique(DEqMS_results_nonCyto_CAT_df$Category), category_order_present),
  category_order_present
)

present_in_results <- unique(DEqMS_results_nonCyto_CAT_df$Category)
category_order_present <- FUN_categories %>%
                            dplyr::filter(category_label %in% present_in_results) %>%
                            dplyr::distinct(category_label) %>%          # keep list order, drop dupes
                            dplyr::pull(category_label)

 category_order <- c( setdiff(present_in_results, category_order_present),category_order_present)


plot_CAT_nonCyto<- DEqMS_results_nonCyto_annotations |> select(Category, gene_callers_id)|> unique() |> 
  group_by(Category) |> 
  summarise(count = n(), .groups = "drop") |> filter(count>=2) |>pull(Category)

  
#plot
DEqMS_results_nonCyto_CAT_df |> filter(Category %in% plot_CAT_nonCyto) |>
  left_join(FUN_categories, by = c("Category" = "category_label")) |> 
  mutate(Category = factor(Category, levels = category_order),
         group = factor(group, levels = names(bacteria_categories_list))) |>
  ggplot(aes(y     = Category,
             x     = log2fold_mean,
             fill  = Enr.reg,
             label = count)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = log2fold_mean - log2fold_se,
                     xmax = log2fold_mean + log2fold_se),
                 height = 0.25, linewidth = 0.4, colour = "grey40") +
  geom_point(aes(size = Prop), shape = 21, alpha = 0.9) +
  geom_text(size = 3.5, #nudge_y = -0.25, 
    colour = "black") +
  scale_size_continuous(range  = c(4, 14),
                        labels = scales::percent_format(accuracy = 0.1)) +
  scale_fill_manual(values = c("WEST" = "#009E73",
                               "GYRE" = "#F0E442",
                               "TRAN" = "#0072B2")) +
  #facet_wrap(~ Enr.reg, nrow = 1) +
  facet_grid(group ~ Enr.reg,scales = "free", space = "free", switch = "y")+
  xlim(c(-8,8))+
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



#save the plot
ggsave("./Figures/nonCyto_enr_proteins_regions.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90,
       scale  = 2,
       dpi = 300)



############################
# TBDRs
############################
DEqMS.results_taxa_df <- DEqMS_results_nonCyto_annotations %>% 
  group_by(Enr.reg, Enr.frac, Category, Class, Order) |>
  summarise(
    count         = n(),
    log2fold_mean = mean(logFC),
    log2fold_se   = se(logFC),
    .groups       = "drop"
  ) |>
  left_join(DEqMS_results_nonCyto_totals, by = c("Enr.reg","Enr.frac")) |>
  mutate(
    Prop    = count / Total_p,
    Enr.reg = factor(Enr.reg, levels = c("WEST","GYRE","TRAN"))
  )



DEqMS_results_nonCyto_annotations |> filter(grepl("TonB-", Category)) |> left_join(tonB_context, by = "gene_callers_id", full) |> 
  select(gene_callers_id, NCBIfam_ann, Function, context=Context, Enr.frac, TCDB_ann, Substrate, logFC, Enr.reg) |> 
  mutate(context=case_when(grepl("siderophore", NCBIfam_ann) ~  "Siderophore",
                           grepl("heme", NCBIfam_ann) ~  "Heme",
                           grepl("N-acetyl-beta-D-glucosamine", Substrate) ~"GlcNAc",
                           grepl("siderophore|Yersiniabactin|ferrienterobactin", Substrate) ~"Siderophore",
                           grepl("cobalamin|cyanocob(III)alamin", Substrate) ~"Cobalamine",
                           grepl("maltooligosaccharide", Substrate) ~"Saccharides",
                           grepl("SusC", Function) & is.na(context)~"Saccharides",
                           grepl("Carbohydrates", context)~"Saccharides",
                           grepl("SusC", context)~"Saccharides",
                           grepl("Zinc", context) ~"Zinc",
                           grepl("Ammonium|Phosphate", context)| is.na(context)~"Undefined",
                           TRUE~context)) |>
  left_join(protein_taxonomy, by = "gene_callers_id") |> 
  #filter(!is.na(Order)) |>
  mutate(Order = case_when(Order %in% c("Cellvibrionales", "SAR86", "Alteromonadales", "Flavobacterales") ~ Order,
                           TRUE ~ "Other taxa")) |>
  group_by(Enr.reg, Enr.frac, Order, context) |>
  summarise(
    count         = n(),
    log2fold_mean = mean(logFC),
    log2fold_se   = se(logFC),
    .groups       = "drop"
  ) |>
    mutate(
    Enr.reg = factor(Enr.reg, levels = c("WEST","GYRE","TRAN"))
  ) |> 
  mutate(Order =factor(Order, levels = c("SAR86", "Cellvibrionales", "Alteromonadales", "Flavobacterales", "Other taxa"))) |> 
  ggplot(aes(y     = context,
             x     = log2fold_mean,
             fill  = Enr.reg,
             label = count)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = log2fold_mean - log2fold_se,
                     xmax = log2fold_mean + log2fold_se),
                 height = 0.25, linewidth = 0.4, colour = "grey40") +
  geom_point(shape = 21, alpha = 0.9, size =8) +
  geom_text(size = 6, #nudge_y = -0.45, 
    colour = "black") +
  scale_size_continuous(range  = c(4, 14),
                        labels = scales::percent_format(accuracy = 0.1)) +
  scale_fill_manual(values = c("WEST" = "#009E73",
                               "GYRE" = "#F0E442",
                               "TRAN" = "#0072B2")) +
  #facet_wrap(~ Enr.reg, nrow = 1) +
  xlim(c(-7,7))+
  facet_grid(Order ~ Enr.reg,scales = "free_y", space = "free", switch = "y")+
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


#save the plot
ggsave("./Figures/TBDRs_enr_proteins_regions.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 180,
       scale  = 2,
       dpi = 300)


DEqMS_results_nonCyto_annotations |> filter(grepl("TonB-|SusC", Category)) |> 
                              group_by(Category, Enr.frac) |> unique() |>
                              summarise(Total_p = n(), .groups = "drop")



DEqMS_results_nonCyto_annotations |> filter(grepl("TonB-", Category)) |> left_join(tonB_context, by = "gene_callers_id", full) |>
  select(gene_callers_id, NCBIfam_ann, Function, context=Context, Enr.frac, TCDB_ann, Substrate, logFC, Enr.reg) |> 
  mutate(context=case_when(grepl("siderophore", NCBIfam_ann) ~  "Siderophore",
                           grepl("heme", NCBIfam_ann) ~  "Heme",
                           grepl("N-acetyl-beta-D-glucosamine", Substrate) ~"GlcNAc",
                           grepl("siderophore|Yersiniabactin|ferrienterobactin", Substrate) ~"Siderophore",
                           grepl("cobalamin|cyanocob(III)alamin", Substrate) ~"Cobalamine",
                           grepl("maltooligosaccharide", Substrate) ~"Oligosacharides",
                           grepl("SusC", Function) & is.na(context)~"SusC",
                           grepl("Zinc", context) ~"Zinc",
                           grepl("Ammonium|Phosphate", context)| is.na(context)~"Undefined",
                           TRUE~context)) |>
  select(gene_callers_id, Enr.frac,  context) |> unique() |> 
  group_by(Enr.frac,  context) |> 
summarise(count = n(), .groups = "drop") 



DEqMS_results_nonCyto_annotations |> filter(grepl("Sus", Category)) |>
  select(Enr.frac,Category, gene_callers_id) |> unique() |>
  group_by(Category, Enr.frac) |> summarise(Total_p = n(), .groups = "drop") 



DEqMS_results_nonCyto_annotations |> filter(grepl("SBP", Category)) |> 
  select(Enr.frac,Category, gene_callers_id) |> unique() |>
  group_by(Category, Enr.frac) |> summarise(Total_p = n(), .groups = "drop")




DEqMS_results_nonCyto_annotations |> filter(grepl("porin", Category)) 


#functional annotation using InterProScan
interpro_annotations<- lapply(c("InterPro", "Pfam","NCBIfam"
                                ),
                                function(src){
                                annotation<- paste0(src,"_ann")
                                accession<- paste0(src,"_acc")
                                annotation_df <- read.table(paste0("data/raw/proteins_",src,".tsv"),quote = "",  sep="\t", h=T) %>% 
                                  dplyr::rename(!!sym(annotation):= function.,
                                                !!sym(accession):= accession) %>% 
                                  select(-c("source")) %>% 
                                  group_by(gene_callers_id)%>%
                                  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
                                  summarise_all(funs(paste(unique(.), collapse='|')))
                                
                                return(annotation_df)
                              })  
#merge Interpro annotations into a table
interpro_annotations_ORFs <- interpro_annotations%>% purrr::reduce(full_join, by = "gene_callers_id") |>
  tidyr::separate(gene_callers_id, into = c("c_prefix", "contig_num", "gene_callers_id"), sep = "_", remove = FALSE) %>%
  tidyr::unite("contig", c_prefix, contig_num, sep = "_", remove = TRUE) |>
  mutate(gene_callers_id = as.character(gene_callers_id)) 


DEqMS_results_nonCyto_annotations |> filter(grepl("TonB-", Category)) |> 
  pull(gene_callers_id) |> unique() -> enriched_tonB_gcids


enriched_contigs<- interpro_annotations_ORFs |> 
                          filter(contig %in% c(interpro_annotations_ORFs |> 
                              group_by(contig) |> summarize(count=n()) |> filter(count>=2) |> pull(contig)))|> 
                              mutate(Enriched=case_when(gene_callers_id %in% enriched_tonB_gcids ~ "yes",
                                                             TRUE~"")) 

enriched_contigs_filt<- enriched_contigs |> 
  filter(contig %in% c(enriched_contigs |> group_by(contig) |> summarize(Enriched_sum=sum(Enriched!="")) |> 
                          filter(Enriched_sum>=1) |> pull(contig)))|>
  select(contig, gene_callers_id, Enriched, names(interpro_annotations_ORFs)) |>
  left_join(DEqMS_results_nonCyto |> select(gene_callers_id, Enr.frac) |> unique(), by = "gene_callers_id") |> 
  left_join(protein_taxonomy , by ="gene_callers_id") |> 
  left_join(DEqMS_results_nonCyto_annotations |> select(gene_callers_id, Category, Function) |> unique(), by = "gene_callers_id") |> 
  select(contig, gene_callers_id, Enriched, Enr.frac,  Category, Function, Class, Order, Family, names(interpro_annotations_ORFs)) 






DEqMS_results_tonB_totals <- DEqMS_results_nonCyto_annotations |> filter(grepl("TonB-", Category)) |> 
                              group_by(Category, Enr.reg, Enr.frac) |> unique() |>
                              summarise(Total_p = n(), .groups = "drop")


DEqMS.results_tonB_taxa_df <- DEqMS_results_nonCyto_annotations |> filter(grepl("TonB-", Category)) |> left_join(tonB_context, by = "gene_callers_id", full) |>
  select(gene_callers_id, Function, context=Context, Enr.frac, TCDB_ann, Substrate, logFC, Enr.reg) |> 
  mutate(context=case_when(grepl("N-acetyl-beta-D-glucosamine", Substrate) ~"GlcNAc",
                           grepl("siderophore|Yersiniabactin|ferrienterobactin", Substrate) ~"Siderophore",
                           grepl("cobalamin|cyanocob(III)alamin", Substrate) ~"Cobalamine",
                           grepl("maltooligosaccharide", Substrate) ~"Oligosacharides",
                           grepl("SusC", Function) & is.na(context)~"SusC",
                           grepl("Zinc", context) ~"Zinc",
                           grepl("Ammonium|Phosphate", context)| is.na(context)~"Undefined",
                           TRUE~context)) |>
  group_by(Enr.reg, Enr.frac,  context) |>
  summarise(
    count         = n(),
    log2fold_mean = mean(logFC),
    log2fold_se   = se(logFC),
    .groups       = "drop"
  ) |>
  left_join(DEqMS_results_tonB_totals, by = c("Enr.reg","Enr.frac")) |>
  mutate(
    Prop    = count / Total_p,
    Enr.reg = factor(Enr.reg, levels = c("WEST","GYRE","TRAN"))
  )

DEqMS_results_tonB_totals <- DEqMS_results_nonCyto_annotations |> filter(grepl("TonB-", Category)) |>
  select(Enr.frac, gene_callers_id) |> unique() |>
  group_by(Enr.frac) |> summarise(Total_p = n(), .groups = "drop") 






############################
# Export results
############################
# one sheet per fraction (BEVs / Cells), sorted by logFC descending
# styled: bold blue header | alternating row fill | frozen top row | auto widths

export_df <- DEqMS_results_nonCyto |>
  filter(Enr.frac != "Not.enr") |>
  left_join(protein_annotations,
            by = "gene_callers_id",
            relationship = "many-to-many") |>
  unique() |>
  left_join(protein_taxonomy, by = "gene_callers_id") |>
  left_join(protein_metadata,  by = "gene_callers_id") |>
  select(any_of(c(
    "Enr.reg","Enr.frac", "gene_callers_id",
    "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
    "sca.P.Value", "sca.adj.pval", "log.sca.pval",
    names(protein_annotations)[str_starts(names(protein_annotations), "InterPro_")],
    names(protein_annotations)[str_starts(names(protein_annotations), "NCBIfam_")],
    names(protein_annotations)[str_starts(names(protein_annotations), "Pfam_")],
    "Domain", "Phylum", "Class", "Order", "Family", "Genus",
    "PSMs", "Unique_pep", "Total_pep", "Length"
  ))) |>
  mutate(
    Enr.reg = factor(Enr.reg, levels = c("TRAN","GYRE","WEST"))
  ) |>
  arrange(Enr.reg, Enr.frac, desc(logFC))

# split into one data frame per comparison
sheets <- split(export_df, export_df$Enr.reg)

# order sheets along the transect
sheet_order <- intersect(c("UP","TRAN","GYRE","WEST"), names(sheets))
sheets <- sheets[sheet_order]

# ── styled workbook ───────────────────────────────────────────────────────────
wb <- openxlsx::createWorkbook()

header_style <- openxlsx::createStyle(
  fontColour     = "#FFFFFF",
  fgFill         = "#2E75B6",
  halign         = "CENTER",
  textDecoration = "bold",
  border         = "Bottom"
)
row_style <- openxlsx::createStyle(fgFill = "#F2F7FB")

purrr::walk(sheet_order, function(Enr.reg) {
  df <- sheets[[Enr.reg]]
  openxlsx::addWorksheet(wb, sheetName = Enr.reg)
  openxlsx::writeData(wb, sheet = Enr.reg, x = df,
                      startRow = 1, startCol = 1,
                      headerStyle = header_style)
  # alternate row fill for every other data row
  even_rows <- seq(3, nrow(df) + 1, by = 2)  # +1 for header offset
  openxlsx::addStyle(wb, sheet = Enr.reg,
                     style = row_style,
                     rows  = even_rows,
                     cols  = seq_len(ncol(df)),
                     gridExpand = TRUE)
  # freeze top row
  openxlsx::freezePane(wb, sheet = Enr.reg, firstRow = TRUE)
  # auto-fit column widths (capped at 40 characters)
  openxlsx::setColWidths(wb, sheet = Enr.reg,
                         cols   = seq_len(ncol(df)),
                         widths = "auto")
})

openxlsx::saveWorkbook(wb,
                       file      = "Tables/Table_S4_DEqMS_results_Fractions_nonCyto.xlsx",
                       overwrite = TRUE)

cat(sprintf(
  "Saved Table_S4 — %d sheets: %s\n  Rows per sheet: %s\n",
  length(sheets),
  paste(sheet_order, collapse = " | "),
  paste(sapply(sheets, nrow), collapse = " | ")
))

  


#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
