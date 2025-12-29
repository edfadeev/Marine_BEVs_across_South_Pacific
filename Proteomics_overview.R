require(dplyr)
require(DEqMS)
require(vegan)
require(tibble)


#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

############################
#import normalized protein abundance table and metadata
############################
samples_df<- read.table("data/samples_meta.txt", header = TRUE) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 

protein_abund<- read.table("data/protein_abund.txt",  header = TRUE)

protein_metadata<- read.table("data/protein_metadata.txt",  header = TRUE)%>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

protein_taxonomy<- read.table("data/protein_taxonomy.txt",  header = TRUE, sep ="\t") %>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

protein_annotations<- read.csv("data/protein_annotations.txt", sep ="\t") %>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

DEqMS_results<- read.table("data/DEqMS_results_regions.txt", h=T)%>% 
  mutate(gene_callers_id=as.character(gene_callers_id))


############################
#Export raw abudnance table with annotations
############################
protein_abund %>% tibble::rownames_to_column("gene_callers_id") %>% 
  left_join(protein_metadata, by ="gene_callers_id") %>% 
  left_join(protein_annotations %>% select(c("gene_callers_id","InterPro_acc","InterPro_ann","NCBIfam_acc","NCBIfam_ann", "Pfam_acc","Pfam_ann", "deeploc")), by ="gene_callers_id") %>% 
  left_join(protein_taxonomy, by ="gene_callers_id") %>% 
  openxlsx::write.xlsx(., colNames = TRUE, 
                       file= 'Tables/Table_S2-Protein_abundance.xlsx')




############################
#Number of identified proteins
############################
#bacterial 
bac_gcids_taxa<- protein_taxonomy %>% 
  filter(Domain=="Bacteria")%>% pull(gene_callers_id)

prot_per_sample <- protein_abund[bac_gcids_taxa,] %>%
  reshape2::melt(variable.name = "Sample_ID", value.name = "Abund") %>% 
  filter(Abund>0) %>% 
  left_join(samples_df, by = "Sample_ID") %>% 
  group_by(Region,Station_ID,Fraction) %>% 
  summarize(N_prot=n()) %>% 
  tidyr::spread(Fraction, N_prot)%>% 
  dplyr::rename("Cells_Bac"="Cells", "BEVs_Bac"="EVs")


#identify viral proteins
vir_gcids_taxa<- protein_taxonomy %>% 
  filter(grepl("Viruses", Domain))%>% pull(gene_callers_id)

vir_gcids_ann <-protein_annotations %>% 
  filter(grepl("phage|virus|capsid|Tail sheath",Pfam_ann, ignore.case = TRUE)|
           grepl("phage|virus|capsid|Tail sheath",InterPro_ann, ignore.case = TRUE)|
           NCBIfam_acc %in% c("TIGR01554", "TIGR02126", "TIGR01543")) %>% pull(gene_callers_id)
vir_gcids<- c(vir_gcids_taxa, vir_gcids_ann)

#summarize how many viral proteins per sample
viral_prot_per_sample <- protein_abund[vir_gcids,] %>%
  reshape2::melt(variable.name = "Sample_ID", value.name = "Abund") %>% 
  filter(Abund>0) %>% 
  left_join(samples_df, by = "Sample_ID") %>% 
  group_by(Region,Station_ID,Fraction) %>% 
  summarize(N_prot=n()) %>% 
  tidyr::spread(Fraction, N_prot) %>% 
  dplyr::rename("Cells_vir"="Cells", "BEVs_vir"="EVs")

#unclassified 
uncl_gcids_taxa<- protein_taxonomy %>% 
  filter(is.na(Domain))%>% pull(gene_callers_id)

uncl_prot_per_sample <- protein_abund[uncl_gcids_taxa,] %>%
  reshape2::melt(variable.name = "Sample_ID", value.name = "Abund") %>% 
  filter(Abund>0) %>% 
  left_join(samples_df, by = "Sample_ID") %>% 
  group_by(Region,Station_ID,Fraction) %>% 
  summarize(N_prot=n()) %>% 
  tidyr::spread(Fraction, N_prot)%>% 
  dplyr::rename("Cells_uncl"="Cells", "BEVs_uncl"="EVs")

#merged table and write into file
protein_abund %>%
  reshape2::melt(variable.name = "Sample_ID", value.name = "Abund") %>% 
  filter(Abund>0) %>% 
  left_join(samples_df, by = "Sample_ID") %>% 
  group_by(Region,Station_ID,Fraction) %>% 
  summarize(N_prot=n()) %>% 
  tidyr::spread(Fraction, N_prot)%>% 
  dplyr::rename("Cells"="Cells", "BEVs"="EVs") %>% 
  left_join(prot_per_sample, by =c("Station_ID", "Region")) %>% 
  left_join(viral_prot_per_sample, by =c("Station_ID", "Region")) %>% 
  left_join(uncl_prot_per_sample, by =c("Station_ID", "Region")) %>% 
  arrange(rev(Station_ID)) %>% 
  write.table(file="Tables/Table_S1-Proteomics_overview.txt", col.names =TRUE, row.names = FALSE, quote = FALSE)

###########################################
#remove viral and unclassified proteins from downstream analysis
###########################################
protein_abund_clean<- protein_abund[!row.names(protein_abund) %in% c(vir_gcids,uncl_gcids_taxa),]

write.table(protein_abund_clean, "data/protein_abund_clean.txt", row.names = TRUE, col.names = TRUE, quote = FALSE)

#log2 transformation of protein abundances
protein_abund.log2 <- log2(as.matrix(protein_abund_clean))

#median normalization of the data
prot.dat.log2_norm = equalMedianNormalization(protein_abund.log2)

write.table(prot.dat.log2_norm, "data/protein_abund_log2_norm.txt", row.names = TRUE, col.names = TRUE, quote = FALSE)

############################
#Protein overlap between fractions
############################
prot_overlap_ls <- lapply(paste0(samples_df$Station_ID,"_"), function(s){
    st<- protein_abund_clean %>% rownames_to_column("gene_callers_id") %>% 
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
prot.dat.log2_norm_no_NA<- prot.dat.log2_norm

prot.dat.log2_norm_no_NA[is.na(prot.dat.log2_norm_no_NA)]<- 0

#test differences between fraction
protein_dist_matrix <- vegdist(t(prot.dat.log2_norm_no_NA), method = "euclidean",na.rm = TRUE)

adonis2(protein_dist_matrix ~ Fraction, samples_df,permutations=999)  
  
  
############################
#print session info and clean the workspace
############################  
sessionInfo()
rm(list = ls())
gc()
  
  
  
  
  
  
  
