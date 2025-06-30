require(dplyr)
require(vegan)
require(stringr)
require(DEqMS)

#load ggplot theme and colours
source("source/ggplot_parameters.R")
source("source/PD_results2R.R")

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

##########################################
#Dissimilarity between cellular proteomes#
##########################################
#replace NAs with 0 for dissimilarity calculation
prot.dat.MDS<- prot.dat.log2_norm[, c(samples_df %>% 
                                        filter(Fraction =="Cells") %>% 
                                        pull(Sample_ID))]
prot.dat.MDS[is.na(prot.dat.MDS)]<- 0
prot.dat.MDS<- prot.dat.MDS[!rowSums(prot.dat.MDS)==0,] # remove proteins that were not observed in any cellular sample


#test whether the differences between the runs are significant
protein_dist_matrix <- vegdist(t(prot.dat.MDS), method = "euclidean",na.rm = TRUE)

adonis2(protein_dist_matrix ~ Region, samples_df %>% filter(Fraction =="Cells"),permutations=999)

#posthoc to check which regions are different
pairwiseAdonis::pairwise.adonis(protein_dist_matrix,factors= samples_df %>% filter(Fraction =="Cells") %>% pull(Region), 
                                p.adjust.m='bonferroni', perm = 999)

##########################################
#Figure 3 - NMDS                         #
##########################################
# Running NMDS in vegan (metaMDS)
prot_NMDS <-  metaMDS(t(prot.dat.MDS),
                      distance = "euclidean",
                      k = 4,
                      maxit = 999, 
                      trymax = 999,
                      wascores = FALSE,
                      autotransform = FALSE,
                      tidy= "sites",
                      na.rm = TRUE)

#plot
prot_NMDS.scores <- as.data.frame(scores(prot_NMDS)) %>% 
  tibble::rownames_to_column(var ="Sample_ID") %>% 
  left_join(samples_df, 
            by = "Sample_ID", copy = TRUE)

prot_NMDS.scores %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, 
             colour = Region, label = Sample_ID))+
  geom_point(size =7, colour = "black")+
  geom_point(size =5)+
  #geom_text(nudge_y = -2, size =4)+
  scale_color_manual(values = c("#009E73", "#F0E442", "#0072B2", 
                                "#D55E00"))+
  theme_EF+
  guides(fill = guide_legend(title="Region", ncol = 4))+
  tune::coord_obs_pred()+
  #coord_fixed()+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Fig_S2-Cell_prot_NMDS.pdf",
       plot = last_plot(),
       units = "mm",
       width = 90,
       height = 90, 
       scale = 2,
       dpi = 300)

##########################################
# Functional differences between region  #
##########################################
enrichment_tests_list <- lapply(c("UP_TRAN","TRAN_GYRE","GYRE_WEST"), function(x) {
  
  samples_meta_sub<- samples_df %>% filter(Fraction =="Cells") %>% 
    mutate(Enr.group = factor(Region, levels =c("WEST","GYRE", "TRAN","UP"))) %>% 
    filter(Enr.group  %in% c(str_split_1(x, "_")))
  
  #prepare experiment design matrix
  cond <-factor(samples_meta_sub$Enr.group, levels = c(stringr::str_split_1(x, "_")))
  design <- model.matrix(~0+cond) # 0 means no intercept for the linear model
  colnames(design) <- gsub("cond","",colnames(design))
  contrast <-  makeContrasts(contrasts=paste(stringr::str_split_1(x, "_")[1], stringr::str_split_1(x, "_")[2], sep ="-"),levels=design)
  
  #replace zero with NA
  prot.dat<- prot.dat.log2_norm[,samples_meta_sub$Sample_ID]
  prot.dat[prot.dat==0] <- NA
  
  #count how many detections were for each protein
  regA_sample_IDs<- samples_meta_sub %>% filter(Enr.group ==str_split_1(x, "_")[1]) %>% pull(Sample_ID)
  regB_sample_IDs<- samples_meta_sub %>% filter(Enr.group ==str_split_1(x, "_")[2]) %>% pull(Sample_ID)
  
  # Filter proteins that were observed in at least two samples in each fraction
  prot.dat.log2_norm.filter<- prot.dat.log2_norm[, c(regA_sample_IDs, regB_sample_IDs)]
  
 
  prot.dat.log2_norm.filter <- prot.dat.log2_norm.filter[rowSums(!is.na(prot.dat.log2_norm.filter[, c(regA_sample_IDs)]))>1 &
                                                           rowSums(!is.na(prot.dat.log2_norm.filter[, c(regB_sample_IDs)]))>1,]

  #run linear model
  fit1<- lmFit(prot.dat.log2_norm.filter, design)
  fit2 <- contrasts.fit(fit1,contrasts = contrast)
  fit3<- eBayes(fit2)
  
  #correct bias of variance estimate based on minimum number of PSMs per protein
  fit3$count <- protein_metadata[rownames(fit3$coefficients), "Number.of.PSMs"]
  
  fit4 = spectraCounteBayes(fit3)
  
  #results
  DEqMS.results <- outputResult(fit4,coef_col = 1) %>% 
    dplyr::rename("gene_callers_id"="gene") %>% 
    mutate(Enr.reg = case_when(logFC>1 & sca.adj.pval<0.1  
                               ~ str_split_1(x, "_")[1],
                               -1> logFC & sca.adj.pval<0.1  
                               ~ str_split_1(x, "_")[2], TRUE ~ "Not.enr"),
           log.sca.pval = -log10(sca.P.Value),
           Comp= x) 
  return(DEqMS.results)
})

#generate results dataframe 
DEqMS_results<- bind_rows(enrichment_tests_list) %>% 
  left_join(protein_annotations %>% select(gene_callers_id,InterPro_acc,InterPro_ann,NCBIfam_acc,NCBIfam_ann,Pfam_acc,Pfam_ann), by = "gene_callers_id", relationship = "many-to-many") %>% unique() %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") 


#overview of total enriched proteins
DEqMS_results %>% filter(Enr.reg!="Not.enr") %>% 
  group_by(Enr.reg) %>% unique() %>% 
  summarize(count=n())

############################
#export results for excel 
############################
DEqMS_results %>% filter(Enr.reg!="Not.enr") %>% 
  select(Comp, Enr.reg, logFC, sca.adj.pval, starts_with("InterPro_"), starts_with("NCBIfam_"),
                         starts_with("Pfam_"), Domain, Phylum, Class,Order,Family,Genus) %>% 
  openxlsx::write.xlsx(., colNames = TRUE, 
                     file= 'Tables/Table_S2-DEqMS_results_regions.xlsx')


#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()