require(dplyr)
require(tidyr)
require(vegan)
require(stringr)
require(DEqMS)

#load ggplot theme and colours
source("source/ggplot_parameters.R")

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#import normalized protein abundance table and metadata
samples_df<- read.table("data/samples_meta.txt", header = TRUE) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 

prot.dat.log2_norm<- read.table("data/protein_abund_log2_norm.txt",  header = TRUE)

protein_metadata<- read.table("data/protein_metadata.txt",  header = TRUE)%>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

protein_taxonomy<- read.table("data/protein_taxonomy.txt",  header = TRUE, sep ="\t") %>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

protein_annotations<- read.csv("data/protein_annotations.txt", sep ="\t") %>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

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
                levels = c("EVs","Cells"))
  design <- model.matrix(~0+cond) # 0 means no intercept for the linear model
  colnames(design) <- gsub("cond","",colnames(design))
  contrast <-  makeContrasts(contrasts="EVs-Cells",levels=design)
  
  #count how many detections were for each protein
  EV_sample_IDs<- samples_meta_sub %>% filter(Fraction =="EVs") %>% pull(Sample_ID)
  Cell_sample_IDs<- samples_meta_sub %>% filter(Fraction =="Cells") %>% pull(Sample_ID)
  
  prot.dat.log2_norm.filter <- prot.dat.log2_norm[,c(EV_sample_IDs,Cell_sample_IDs)]
  
  # Filter proteins that were observed in at least two samples in each fraction
  prot.dat.log2_norm.filter <- prot.dat.log2_norm.filter[rowSums(!is.na(prot.dat.log2_norm.filter[, c(EV_sample_IDs)]))>1 &
                                                           rowSums(!is.na(prot.dat.log2_norm.filter[, c(Cell_sample_IDs)]))>1,]
  
  #run linear model
  fit1<- lmFit(prot.dat.log2_norm.filter, design)
  fit2 <- contrasts.fit(fit1,contrasts = contrast)
  fit3<- eBayes(fit2)
  
  #correct bias of variance estimate based on number of PSMs per protein
  fit3$count <- protein_metadata[rownames(fit3$coefficients), "Number.of.PSMs"]

  fit4 = spectraCounteBayes(fit3)
  
  #results
  DEqMS.results <- outputResult(fit4,coef_col = 1) %>% 
                 dplyr::rename("gene_callers_id"="gene") %>% 
                 mutate(Enr.frac = case_when(logFC>1 & sca.adj.pval<0.1  
                                             ~ "BEVs",
                                             -1> logFC & sca.adj.pval<0.1  
                                             ~ "Cells", TRUE ~ "Not.enr"),
                        log.sca.pval = -log10(sca.P.Value),
                        Enr.group =x) 
  
  return(DEqMS.results)
})

############################
#enrichment results 
############################
DEqMS_results<- bind_rows(enrichment_tests_list) %>% 
  filter(Enr.frac!="Not.enr")

write.table(DEqMS_results, "data/DEqMS_results_regions.txt", col.names =T, row.names = F, quote = F)

#total of different proteins
DEqMS_results %>%  
  select(Enr.frac, gene_callers_id) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())

#total of different orders
DEqMS_results %>% select(Enr.frac, gene_callers_id) %>% unique() %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  select(Enr.frac, Class, Order) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_orders = n())

#proteins per order
DEqMS_results %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  select(Enr.frac, Phylum, Class, Order, gene_callers_id) %>% unique() %>% 
  group_by(Enr.frac, Phylum, Class, Order) %>% 
  summarize(N_p = n()) %>% 
  tidyr::spread(Enr.frac, N_p) %>% View()

#proteins per order per region
DEqMS_results %>% 
  left_join(protein_taxonomy %>% unique(), by = "gene_callers_id") %>% 
  select(Enr.group, Enr.frac, Phylum, Class, Order, gene_callers_id) %>% unique() %>% 
  group_by(Enr.group, Enr.frac, Phylum, Class, Order) %>% 
  summarize(N_p = n()) %>% 
  tidyr::spread(Enr.frac, N_p) %>% View()

#total of proteins epr region
DEqMS_results_totals<- DEqMS_results %>% 
  select(Enr.group, Enr.frac, gene_callers_id) %>% 
  group_by(Enr.group, Enr.frac) %>% 
  summarize(Total_p = n())

############################
#taxonomy overview of enriched proteins
############################
#plot taxa of the enriched proteins
DEqMS_results_Class_tax<- DEqMS_results %>% 
  left_join(protein_taxonomy , by = "gene_callers_id") %>% 
  group_by(Enr.group, Enr.frac, Class, Order) %>% 
  summarize(N_p = n()) 

DEqMS_results_Class_totals<- DEqMS_results_Class_tax %>% 
  group_by(Enr.group, Enr.frac) %>% summarize(Total_p = sum(N_p)) 

DEqMS_results_Class_tax %>% 
  left_join(DEqMS_results_Class_totals, by = c("Enr.group", "Enr.frac")) %>% 
  mutate(Prop=100*N_p/Total_p, Order=gsub("^ ","",Order),Class=gsub("^ ","",Class) ) %>% 
  mutate(Taxa=case_when(Prop<1~"Other taxa < 1%", 
                        is.na(Order) & is.na(Class)~"unclassified",
                        is.na(Order) & !is.na(Class) ~ paste0(Class, "_uncl"),
                        grepl("Pelagibacterales", Order)~"Pelagibacterales",
                        Order=="Candidatus Puniceispirillales"~"Puniceispirillales",
                        TRUE~Order)) %>%  
  mutate(Taxa = factor(Taxa, levels=c( "Caulobacterales", 
                                       "Hyphomicrobiales" , 
                                       "Hyphomonadales","Kordiimonadales", "Minwuiales","Parvularculales","Pelagibacterales" ,"Puniceispirillales",
                                       "Rhodobacterales" , "Rhodospirillales" , "Sphingomonadales" , "Alphaproteobacteria_uncl", "Burkholderiales", 
                                       "Chitinophagales", "Cytophagales" , "Synechococcales", "Flavobacteriales" ,
                                       "Alteromonadales" , "Cellvibrionales" , "Chromatiales", "Kangiellales", "Lysobacterales" ,"Moraxellales",
                                       "Oceanospirillales" , "Pseudomonadales" , "Sphingobacteriales", "Gammaproteobacteria_uncl", 
                                       "Other taxa < 1%", "Viral", "unclassified")),
         Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN","UP")),
         Enr.frac=ifelse(Enr.frac=="EVs","BEVs",Enr.frac)) %>%  
  ggplot(aes(x=Enr.frac, y= Prop, fill = Taxa))+
  geom_col()+
  scale_fill_manual(values = tol21rainbow)+
  facet_grid(cols=vars(Enr.group),scales="free_x",space="free_x",switch="x") +
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/prot_enr_Tax_regions.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

############################
#export results for excel 
############################
DEqMS_results %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(protein_annotations %>% unique(), by = "gene_callers_id", relationship = "many-to-many") %>% unique() %>% 
  left_join(protein_taxonomy%>% unique(), by = "gene_callers_id") %>% 
  left_join(protein_metadata, by="gene_callers_id") %>% 
  mutate(deeploc=case_when(deeploc=="Cell wall & surface"~"Extracellular", TRUE~ deeploc)) %>% 
  select(names(DEqMS_results), Phylum, Class,Order,Family,Genus, deeploc, starts_with("InterPro_"), starts_with("NCBIfam_"),
         starts_with("Pfam_"),Number.of.PSMs,Number.of.Unique.Peptide,Number.of.Razor.Peptide,Length, aa_sequence) %>% 
  openxlsx::write.xlsx(., colNames = TRUE, 
                       file= 'Tables/Table_S4-DEqMS_regions.xlsx')



#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
