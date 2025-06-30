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
#Explore viral proteins in BEVs fraction
############################
vir_gcids_taxa<- protein_taxonomy %>% 
  filter(grepl(".*viri|Unknown_Viruses", Phylum_refseq))%>% pull(gene_callers_id)

vir_gcids_ann <-protein_metadata %>% filter(grepl("phage|virus|capsid|Tail sheath",Pfam_ann, ignore.case = TRUE)|
                                               grepl("phage|virus|capsid|Tail sheath",InterPro_ann, ignore.case = TRUE)|
                                               NCBIfam_acc %in% c("TIGR01554", "TIGR02126", "TIGR01543")) %>% pull(gene_callers_id)

#summarize how many viral proteins per sample
protein_abund[c(vir_gcids_taxa, vir_gcids_ann),] %>%
  reshape2::melt(variable.name = "Sample_ID", value.name = "Abund") %>% 
  filter(Abund>0) %>% 
  left_join(samples_df, by = "Sample_ID") %>% 
  group_by(Region,Station_ID,Fraction) %>% 
  summarize(N_prot=n()) %>% 
  tidyr::spread(Fraction, N_prot)

############################
#Protein overlap between fractions
############################
prot_overlap_ls <- lapply(paste0(samples_df$Station_ID,"_"), function(s){
  st<- protein_abund %>% rownames_to_column("gene_callers_id") %>% 
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
  summarize(Mean = mean(EVs_shared_prop),
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

############################
#Protein enrichment analysis between fractions in each region
############################
enrichment_tests_list <- lapply(c("UP","TRAN","GYRE","WEST"), function(x) {
  
  samples_meta_sub<- samples_df %>% 
    mutate(Enr.group = factor(Region, levels =c("WEST","GYRE", "TRAN","UP"))) %>% 
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
  
  prot.dat.log2_norm.filter<- prot.dat.log2_norm[, c(EV_sample_IDs, Cell_sample_IDs)]
  
  # Filter proteins that were observed in at least two samples in each fraction
  prot.dat.log2_norm.filter <- prot.dat.log2_norm.filter[rowSums(!is.na(prot.dat.log2_norm.filter[, c(EV_sample_IDs)]))>1 &
                                rowSums(!is.na(prot.dat.log2_norm.filter[, c(Cell_sample_IDs)]))>1,]

    #run linear model
  fit1<- lmFit(prot.dat.log2_norm.filter, design)
  fit2 <- contrasts.fit(fit1,contrasts = contrast)
  fit3<- eBayes(fit2)
  
  #correct bias of variance estimate based on number of PSMs per protein
  fit3$count <- protein_metadata %>% filter(gene_callers_id %in% rownames(fit3$coefficients)) %>% pull(Number.of.PSMs)
  
  fit4 = spectraCounteBayes(fit3)
  
  #results
  if(x=="UP"){
    DEqMS.results <- outputResult(fit4,coef_col = 1) %>% 
      dplyr::rename("gene_callers_id"="gene") %>% 
      mutate(Enr.frac = case_when(logFC>1 & sca.adj.pval<0.2  
                                  ~ "EVs",
                                  -1> logFC & sca.adj.pval<0.2  
                                  ~ "Cells", TRUE ~ "Not.enr"),
             log.sca.pval = -log10(sca.P.Value),
             Enr.group =x) } else {
               DEqMS.results <- outputResult(fit4,coef_col = 1) %>% 
                 dplyr::rename("gene_callers_id"="gene") %>% 
                 mutate(Enr.frac = case_when(logFC>1 & sca.adj.pval<0.1  
                                             ~ "EVs",
                                             -1> logFC & sca.adj.pval<0.1  
                                             ~ "Cells", TRUE ~ "Not.enr"),
                        log.sca.pval = -log10(sca.P.Value),
                        Enr.group =x) 
             }
  
  
  return(DEqMS.results)
})

############################
#enrichment results 
############################
DEqMS_results<- bind_rows(enrichment_tests_list) %>% 
  filter(Enr.frac!="Not.enr") %>% 
  left_join(protein_metadata , by = "gene_callers_id")

#total of different proteins
DEqMS_results %>% 
  select(Enr.frac, gene_callers_id) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())

#summarize number of enriched proteins
DEqMS_results_totals<- DEqMS_results %>% 
  group_by(Enr.group, Enr.frac) %>% 
  summarize(Total_p = n())

#check enrichment of viral proteins
DEqMS_results_vir<- DEqMS_results %>% filter(gene_callers_id %in% c(vir_gcids_taxa, vir_gcids_ann)) 

DEqMS_results_vir%>% 
  group_by(Enr.group, Enr.frac) %>% 
  summarize(n())

#export results for excel 
openxlsx::write.xlsx(DEqMS_results, 
                     file= 'data/DEqMS_results_1704.xlsx')

#generate faa file
enr_prot_AAs<- DEqMS_results %>% 
  select(gene_callers_id, aa_sequence) %>% unique() %>% 
  mutate(names=paste(">",gene_callers_id, sep="")) %>% 
  select(names, aa_sequence)

enr_prot_AAs_fasta <- c(rbind(enr_prot_AAs$names, enr_prot_AAs$aa_sequence))

write(x = enr_prot_AAs_fasta, file = "data/enriched_prot_AAs.faa")

############################
# localization of the enriched proteins
############################
DEqMS_results %>% 
  mutate(PSORTb_ann = case_when(grepl("TRANSMEMBRANE", Phobius_acc) ~ "Transmembranal",
                                is.na(PSORTb_ann) & grepl("NON_CYTOPLASMIC_DOMAIN", Phobius_acc) & !grepl("in the cytoplasm", Phobius_ann) ~ "OuterMembrane/Extracellular",
                                is.na(PSORTb_ann) & grepl("in the cytoplasm", Phobius_ann) & !grepl("NON_CYTOPLASMIC_DOMAIN", Phobius_acc) ~ "CytoplasmicMembrane",
                                TRUE ~ PSORTb_ann)) %>% 
  mutate(PSORTb_ann = case_when(grepl("Extracellular|OuterMembrane|Periplasmic|Transmembranal|CytoplasmicMembrane", PSORTb_ann) ~ "Non-cytoplasmic", 
                                grepl("Cytoplasmic", PSORTb_ann) ~ "Cytoplasmic", 
                                TRUE ~  PSORTb_ann)) %>% 
  group_by(Enr.frac, PSORTb_ann) %>% 
  summarize(Total_p = n()) %>% 
  tidyr::spread(Enr.frac, Total_p)%>% 
  mutate(Total = sum(Cells, EVs, Not.enr))


#plot localization of the enriched proteins
DEqMS_results_loc<- DEqMS_results %>% 
  mutate(PSORTb_ann = case_when(grepl("TRANSMEMBRANE", Phobius_acc) ~ "Transmembranal",
                                is.na(PSORTb_ann) & grepl("NON_CYTOPLASMIC_DOMAIN", Phobius_acc) & !grepl("in the cytoplasm", Phobius_ann) ~ "OuterMembrane/Extracellular",
                                is.na(PSORTb_ann) & grepl("in the cytoplasm", Phobius_ann) & !grepl("NON_CYTOPLASMIC_DOMAIN", Phobius_acc) ~ "CytoplasmicMembrane",
                                TRUE ~ PSORTb_ann)) %>% 
  mutate(PSORTb_ann = case_when(grepl("Extracellular|OuterMembrane", PSORTb_ann) ~ "OuterMembrane/Extracellular", 
                                grepl("Cytoplasmic", PSORTb_ann) ~ "Cytoplasmic", 
                                TRUE ~  PSORTb_ann),
         Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN","UP"))) %>% 
  group_by(Enr.group, Enr.frac, PSORTb_ann) %>% 
  summarize(N_p = n()) 

DEqMS_results_loc_totals<- DEqMS_results_loc %>% 
  group_by(Enr.group, Enr.frac) %>% summarize(Total_p = sum(N_p)) 

DEqMS_results_loc %>% 
  mutate(PSORTb_ann =ifelse(is.na(PSORTb_ann),"Undefined", PSORTb_ann)) %>% 
  left_join(DEqMS_results_loc_totals, by = c("Enr.group", "Enr.frac")) %>% 
  mutate(Prop=100*N_p/Total_p,
         Enr.frac=ifelse(Enr.frac=="EVs","BEVs",Enr.frac)) %>%  
  ggplot(aes(x=Enr.frac, y= Prop, fill = PSORTb_ann))+
  geom_col()+
  geom_text(aes(label=Total_p), y=100, vjust=0.1, size = 5)+
  scale_fill_manual(values = c("#e5c185","#c7522a", "#74a892", "#008585", "#fbf2c4"))+
  facet_grid(cols=vars(Enr.group),scales="free_x",space="free_x",switch="x") +
  ylab("Proportion of enriched proteins (%)")+
  theme_EF+
  guides(fill = guide_legend(title="Sub-cellular localization", nrow = 2))+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

#save the plot
ggsave("./Figures/prot_enr_local.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

############################
#taxonomy overview of enriched proteins
############################
#plot taxa of the enriched proteins
DEqMS_results_Class_tax<- DEqMS_results %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  group_by(Enr.frac, Class, Order) %>% unique() %>% 
  summarize(N_p = n()) 



#plot taxa of the enriched proteins
DEqMS_results_Class_tax<- DEqMS_results %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  group_by(Enr.group, Enr.frac, Class, Order) %>% 
  summarize(N_p = n()) 

DEqMS_results_Class_totals<- DEqMS_results_Class_tax %>% 
  group_by(Enr.group, Enr.frac) %>% summarize(Total_p = sum(N_p)) 

DEqMS_results_Class_tax %>% 
  left_join(DEqMS_results_Class_totals, by = c("Enr.group", "Enr.frac")) %>% 
  mutate(Prop=N_p/Total_p, Order=gsub("^ ","",Order),Class=gsub("^ ","",Class) ) %>% 
  mutate(Taxa=case_when(Prop<0.01~"Other taxa < 1%", 
                        is.na(Order)~"unclassified",
                        Order=="NA" ~"unclassified",
                        Order=="Candidatus Pelagibacterales"~"Pelagibacterales",
                        Order=="Candidatus Puniceispirillales"~"Puniceispirillales",
                        TRUE~Order)) %>% 
  mutate(Taxa = factor(Taxa, levels=c( "Hyphomicrobiales" , "Minwuiales","Parvularculales","Pelagibacterales" ,"Puniceispirillales",
                                       "Rhodobacterales" , "Rhodospirillales" , "Sphingomonadales" , "Burkholderiales", 
                                       "Cytophagales" , "Synechococcales", "Flavobacteriales" ,
                                       "Alteromonadales" , "Cellvibrionales" , "Lysobacterales" ,
                                       "Oceanospirillales" , "Pseudomonadales" , "Xanthomonadales" , "Unknown_Pseudomonadota", "Unknown_Uroviricota",
                                       "Other taxa < 1%", "unclassified")),
         Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN","UP")),
         Enr.frac=ifelse(Enr.frac=="EVs","BEVs",Enr.frac)) %>%   
  ggplot(aes(x=Enr.frac, y= Prop, fill = Taxa))+
  geom_col()+
  #geom_text(aes(label=Total_p), y=1, vjust=0.1, size = 5)+
  scale_fill_manual(values = c(tol21rainbow,"gray50"))+
  facet_grid(cols=vars(Enr.group),scales="free_x",space="free_x",switch="x") +
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/prot_enr_Tax.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

############################
#iron related porins in cyanobacteria
############################
DEqMS_results %>%left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  filter(grepl("Cyanophyceae", Class)) %>% 
  mutate(Iron=ifelse(grepl("porin|ferr", InterPro_ann, ignore.case = TRUE),"Yes","No")) %>% 
  group_by(Enr.group, Enr.frac, Iron) %>% 
  summarize(Total_p = n()) %>% 
  tidyr::spread(Enr.frac, Total_p)
