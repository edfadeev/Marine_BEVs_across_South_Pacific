require(dplyr)
require(tibble)

############################
#import metadata
############################
samples_df<- read.table("./11_PROTEIN/R_data/samples_meta.txt", sep ="\t", header = TRUE, row.names = 1) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 

############################
#Import PD output
############################
protein_metadata<- read.table("./11_PROTEIN/SO289_clust99_RC_PeptideGroups.txt", sep ="\t", header = TRUE) %>% 
  dplyr::rename(gene_callers_id=Master.Protein.Accessions) %>% 
  tidyr::separate_rows(gene_callers_id, sep= "; ") %>% 
  #calculate total of unique peptide abundance by protein
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(all_of(samples_df$File.ID), ~ samples_df$Sample_ID)%>% 
  filter(Number.of.PSMs>1, gene_callers_id!="") %>% #include only peptides with at least two PSMs
  mutate(gene_callers_id=as.character(gene_callers_id),
         Unique.Peptides = case_when(grepl(";",Positions.in.Master.Proteins)~ 0, TRUE ~ 1)) %>% 
  group_by(gene_callers_id) %>% 
  summarize(Number.of.Peptides=n(), Number.of.PSMs=sum(Number.of.PSMs), Number.of.Unique.Peptide=sum(Unique.Peptides)) %>% 
  filter(Number.of.Unique.Peptide> 0) # remove proteins with less than 2 peptides and without unique peptides 

############################
#save gene callers IDs of the identified proteins for further annotation
############################
write(paste(protein_metadata$gene_callers_id, collapse = ','), file ="./11_PROTEIN/R_data/prot_GCIDs.txt")
