require(dplyr)
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
############################
<<<<<<< HEAD
#import normalized protein abundance table and metadata
############################
samples_df<- read.table("data/samples_meta.txt", header = TRUE) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 

prot.dat.log2_norm<- read.table("data/protein_abund_log2_norm.txt",  header = TRUE)

protein_metadata<- read.table("data/protein_metadata.txt",  header = TRUE)%>% 
  mutate(gene_callers_id=as.character(gene_callers_id))
=======
#Protein overlap between fractions
############################
prot_overlap_ls <- lapply(paste0(samples_df$Station_ID,"_"), function(s){
  st<- protein_abund_no_vir %>% rownames_to_column("gene_callers_id") %>% 
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
  summarize(Min=min(EVs_shared_prop),
            Max=max(EVs_shared_prop),
            Mean = mean(EVs_shared_prop),
            SE = se(EVs_shared_prop))

############################
#dissimilarity of Cells vs. BEV proteomes
############################
#replace NAs with 0 for dissimilarity calculation
prot.dat.MDS<- prot.dat.log2_norm

prot.dat.MDS[is.na(prot.dat.MDS)]<- 0
prot.dat.MDS<- prot.dat.MDS[!rowSums(prot.dat.MDS)==0,] # remove proteins that were not observed in any cellular sample


#test whether the differences between the runs are significant
protein_dist_matrix <- vegdist(t(prot.dat.MDS), method = "euclidean",na.rm = TRUE)

adonis2(protein_dist_matrix ~ Fraction, samples_df,permutations=999)
>>>>>>> e48d2ba4b277485019716045c8d9a68ab236f549

protein_taxonomy<- read.table("data/protein_taxonomy.txt",  header = TRUE, sep ="\t") %>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

protein_annotations<- read.csv("data/protein_annotations.txt", sep ="\t") %>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

DEqMS_results<- read.table("data/DEqMS_results_regions.txt", h=T)%>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

############################
#Protein enrichment analysis between fractions
############################
#prepare experiment design matrix
cond <-factor(samples_df$Fraction,
              levels = c("EVs","Cells"))
design <- model.matrix(~0+cond) # 0 means no intercept for the linear model
colnames(design) <- gsub("cond","",colnames(design))
contrast <-  makeContrasts(contrasts="EVs-Cells",levels=design)

#count how many detections were for each protein
EV_sample_IDs<- samples_df %>% filter(Fraction =="EVs") %>% pull(Sample_ID)
Cell_sample_IDs<- samples_df %>% filter(Fraction =="Cells") %>% pull(Sample_ID)

<<<<<<< HEAD
=======
prot.dat.log2_norm.filter<- prot.dat.log2_norm[, c(EV_sample_IDs, Cell_sample_IDs)]

>>>>>>> e48d2ba4b277485019716045c8d9a68ab236f549
# Filter proteins that were observed in at least two samples in each fraction
prot.dat.log2_norm.filter <- prot.dat.log2_norm[rowSums(!is.na(prot.dat.log2_norm[, c(EV_sample_IDs)]))>1 &
                                                         rowSums(!is.na(prot.dat.log2_norm[, c(Cell_sample_IDs)]))>1,]

#run linear model
fit1<- lmFit(prot.dat.log2_norm.filter, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3<- eBayes(fit2)

#correct bias of variance estimate based on number of PSMs per protein
fit3$count <- protein_metadata[rownames(fit3$coefficients), "Number.of.PSMs"]

fit4 = spectraCounteBayes(fit3)

DEqMS.results <- outputResult(fit4,coef_col = 1) %>% 
  dplyr::rename("gene_callers_id"="gene") %>% 
  mutate(Enr.frac = case_when(logFC>1 & sca.adj.pval<0.1  
                              ~ "BEVs",
                              -1> logFC & sca.adj.pval<0.1  
                              ~ "Cells", TRUE ~ "Not.enr")) %>% 
  left_join(protein_metadata , by = "gene_callers_id")

#total of different proteins
DEqMS.results %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())

############################
# Export results
############################
#export results for excel 
DEqMS.results %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(protein_annotations %>% unique(), by = "gene_callers_id", relationship = "many-to-many") %>% unique() %>% 
  left_join(protein_taxonomy%>% unique(), by = "gene_callers_id") %>% 
  mutate(deeploc=case_when(deeploc=="Cell wall & surface"~"Extracellular", TRUE~ deeploc)) %>% 
  select(names(DEqMS.results), Phylum, Class,Order,Family,Genus, deeploc, starts_with("InterPro_"), starts_with("NCBIfam_"),
         starts_with("Pfam_")) %>% 
  openxlsx::write.xlsx(., colNames = TRUE, 
                       file= 'Tables/Table_S3-DEqMS_results_Fractions.xlsx')

############################
# localization of the enriched proteins
############################
#plot localization of the enriched proteins
DEqMS_results_loc<- DEqMS.results %>% 
  left_join(protein_annotations,  by="gene_callers_id") %>% 
  mutate(deeploc=case_when(deeploc=="Cell wall & surface"~"Extracellular", TRUE~ deeploc)) %>% 
  group_by(Enr.frac, deeploc) %>% 
  summarize(N_p = n()) 
  
DEqMS_results_loc_totals<- DEqMS_results_loc %>% 
  group_by(Enr.frac) %>% summarize(Total_p = sum(N_p)) 

DEqMS_results_loc %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(DEqMS_results_loc_totals, by = c("Enr.frac")) %>% 
  mutate(Prop=100*N_p/Total_p,
         Enr.frac=ifelse(Enr.frac=="EVs","BEVs",Enr.frac),
         deeploc=factor(deeploc, levels=c("Cytoplasmic","Cytoplasmic Membrane",
                                          "Periplasmic", "Outer Membrane",
                                          "Extracellular"))) %>%  
  ggplot(aes(x=Enr.frac, y= Prop, fill = deeploc))+
  geom_col()+
  #geom_text(aes(label=Total_p), y=100, vjust=0.1, size = 5)+
  scale_fill_manual(values = c("#e5c185","#c7522a", "#74a892", "#008585", "#fbf2c4","red4"))+
  #facet_grid(cols=vars(Enr.group),scales="free_x",space="free_x",switch="x") +
  ylab("Proportion of enriched proteins (%)")+
  theme_EF+
  guides(fill = guide_legend(title="Sub-cellular localization", ncol = 2))+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

#save the plot
ggsave("./Figures/prot_enr_local.pdf",
       plot = last_plot(),
       units = "mm",
       width = 90,
       #height = 90, 
       scale = 2,
       dpi = 300)


############################
# taxonomy of the enriched proteins
############################
#summarize enriched proteins by genera
DEqMS_results_Genus_tax<- DEqMS.results %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  group_by(Enr.frac, Phylum, Class, Order, Family, Genus) %>% 
  summarize(N_p = n()) 

#plot taxa of the enriched proteins
DEqMS_results_Class_tax<- DEqMS.results %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  group_by(Enr.frac, Phylum, Class, Order) %>% 
  summarize(N_p = n()) 

DEqMS_results_Class_tax %>% filter(Enr.frac!="Not.enr") %>% tidyr::spread(Enr.frac, N_p) %>% View()

DEqMS_results_Class_totals<- DEqMS_results_Class_tax %>% 
  group_by(Enr.frac) %>% summarize(Total_p = sum(N_p)) 

DEqMS_results_Class_tax %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(DEqMS_results_Class_totals, by = "Enr.frac") %>% 
  mutate(Prop=100*N_p/Total_p, Order=gsub("^ ","",Order),Class=gsub("^ ","",Class) ) %>% 
  mutate(Taxa=case_when(Prop<1~"Other taxa < 1%", 
                        is.na(Order) & is.na(Class)~"unclassified",
                        is.na(Order) & !is.na(Class) ~ paste0(Class, "_uncl"),
                        grepl("Pelagibacterales", Order)~"Pelagibacterales",
                        Order=="Candidatus Puniceispirillales"~"Puniceispirillales",
                        TRUE~Order)) %>% 
  mutate(Taxa = factor(Taxa, levels=c( "Hyphomicrobiales" , "Minwuiales","Parvularculales","Pelagibacterales" ,"Puniceispirillales",
                                       "Rhodobacterales" , "Rhodospirillales" , "Sphingomonadales" , "Alphaproteobacteria_uncl", "Burkholderiales", 
                                       "Cytophagales" , "Synechococcales", "Flavobacteriales" ,
                                       "Alteromonadales" , "Cellvibrionales" , "Lysobacterales" , "Moraxellales",
                                       "Oceanospirillales" , "Pseudomonadales" , "Xanthomonadales" , "Gammaproteobacteria_uncl", "Pseudomonadota", "Unknown_Uroviricota",
                                       "Other taxa < 1%", "unclassified")),
         Enr.frac=ifelse(Enr.frac=="EVs","BEVs",Enr.frac)) %>%   
  ggplot(aes(x=Enr.frac, y= Prop, fill = Taxa))+
  geom_col()+
  #geom_text(aes(label=Total_p), y=101, vjust=0.1, size = 5)+
  scale_fill_manual(values = c(tol21rainbow))+
  ylab("Proportion of enriched proteins (%)")+
  theme_EF+
  guides(fill = guide_legend(nrow = 6))+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

#save the plot
ggsave("./Figures/prot_enr_taxa.pdf",
       plot = last_plot(),
       units = "mm",
       width = 90,
       #height = 90, 
       scale = 2,
       dpi = 300)


#summarize enriched proteins by Order
DEqMS_results_Order_tax<- DEqMS.results %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  group_by(Enr.frac, Phylum, Class, Order) %>% 
  summarize(N_p = n()) %>% 
  tidyr::spread(Enr.frac, N_p)

#plot ratio of enriched proteins
DEqMS_results_Order_tax %>% 
  filter(!is.na(BEVs)&!is.na(Cells)&!is.na(Phylum)) %>%  
  mutate(Order=case_when(is.na(Order) & !is.na(Class)~ paste0(Class,"_uncl"),
                         is.na(Order) & is.na(Class)~ paste0(Phylum,"_uncl"),TRUE~Order),
         Class=case_when(is.na(Class)~ paste0(Phylum,"_uncl"),TRUE~Class)) %>% 
  mutate(Ratio=BEVs/Cells,
         Taxa=case_when(BEVs>2&Cells>2 ~ Order, TRUE~"Other taxa"),
         Class=case_when(BEVs>2&Cells>2 ~ Class, TRUE~"Other taxa")) %>% 
  mutate(Taxa=factor(Taxa, levels=c(DEqMS_results_Order_tax %>% filter(!is.na(Order)) %>% pull(Order) %>% unique(), 
                                    "Alphaproteobacteria_uncl", "Gammaproteobacteria_uncl", "Bacteroidota_uncl", "Other taxa" )),
         Class=factor(Class, levels=c(DEqMS_results_Order_tax %>% filter(!is.na(Class)) %>% pull(Class) %>% unique() %>% sort(), 
                                      "Bacteroidota_uncl", "Other taxa" ))) %>% 
  mutate(Taxa=gsub("Candidatus ","", Taxa)) %>% 
  ggplot(aes(x=Cells, y= BEVs, colour= Class, label = Taxa))+
  geom_point(size=5)+
  geom_label_repel(color = "red", nudge_x = .5, nudge_y = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.8)+
  scale_colour_manual(values = c(rev(tol21rainbow)))+
  labs(x="Enriched proteins in cells", y = "Enriched proteins in BEVs")+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/prot_cells_BEVs_ratio.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 180, 
       scale = 2,
       dpi = 300)

############################
# Export results
############################
#export results for excel 
DEqMS.results %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(protein_annotations, by = "gene_callers_id", relationship = "many-to-many") %>% unique() %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  select(Enr.frac, logFC, sca.adj.pval, starts_with("InterPro_"), starts_with("NCBIfam_"),
         starts_with("Pfam_"), Domain, Phylum, Class,Order,Family,Genus) %>% 
  openxlsx::write.xlsx(., colNames = TRUE, 
                       file= 'Tables/Table_S3-DEqMS_results_Fractions.xlsx')

  
  
#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
