require(dplyr)
require(tibble)
require(ggplot2)
require(DEqMS)

#load ggplot theme and colours
source("source/ggplot_parameters.R")

############################
#import metadata
############################
samples_df<- read.table("data/samples_meta.txt", sep ="\t", header = TRUE, row.names = 1) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 

############################
#Import PD output
############################
protein_metadata<- read.table("data/PD_peptide_df.txt", sep ="\t", header = TRUE) %>% 
  filter(Number.of.PSMs>1, gene_callers_id!="") %>% #include only peptides with at least two PSMs
  mutate(gene_callers_id=as.character(gene_callers_id),
         Unique.Peptides = case_when(grepl(";",Positions.in.Master.Proteins)~ 0, TRUE ~ 1)) %>% 
  group_by(gene_callers_id) %>% 
  summarize(Number.of.Peptides=n(), Number.of.PSMs=sum(Number.of.PSMs), Number.of.Unique.Peptide=sum(Unique.Peptides)) %>% 
  filter(Number.of.Peptides>1, Number.of.Unique.Peptide> 0) %>% # remove proteins with less than 2 peptides and without unique peptides 
  left_join(read.table("data/protein_ann.txt", sep ="\t", header = TRUE, colClasses = "character"), by = "gene_callers_id")

#summarize abundance of peptides by protein
protein_abund<- read.table("data/PD_peptide_df.txt", sep ="\t", header = TRUE) %>% 
  filter(gene_callers_id %in% protein_metadata$gene_callers_id) %>% #include only proteins that match the filtering criteria
  group_by(gene_callers_id) %>% 
  summarise_at(vars(samples_df$Sample_ID),sum, na.rm = TRUE) %>% 
  tibble::column_to_rownames("gene_callers_id")%>% 
  mutate(across(where(is.numeric), ~na_if(., 0)))


#log2 transformation of protein abundances
protein_abund.log2 <- log2(as.matrix(protein_abund))

#median normalization of the data
prot.dat.log2_norm = equalMedianNormalization(protein_abund.log2)

############################
#Import protein taxonomy
############################
#based on NCBI Refseq 
protein_taxonomy_refseq<- protein_metadata %>% select(gene_callers_id, aa_sequence) %>%
  left_join(read.csv("data/genes_taxonomy_refseq.txt", sep="\t", h= T, colClasses = "character"), by ="gene_callers_id")  %>% 
  left_join(read.csv("data/taxon_names_refseq.txt", sep="\t", h= T, colClasses = "character") , by = "taxon_id") %>% 
  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
  dplyr::rename(Class = t_order ,Order = t_class, #switch columns due to some taxonomy parsing bug in anvio
                Phylum = t_phylum, Family = t_family,
                Genus = t_genus, Species=t_species) %>% 
  mutate(Order= case_when(grepl("Pelagibacterales", Order)~ "Pelagibacterales", TRUE ~ Order)) #correcting some space character in the name


protein_taxonomy_nr<- protein_metadata %>% select(gene_callers_id, aa_sequence) %>%
  left_join(read.csv("data/genes_taxonomy_nr.txt", sep="\t", h= T, colClasses = "character"), by ="gene_callers_id")  %>% 
  left_join(read.csv("data/taxon_names_nr.txt", sep="\t", h= T, colClasses = "character") , by = "taxon_id") %>% 
  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
  dplyr::rename(Class = t_order ,Order = t_class, #switch columns due to some taxonomy parsing bug in anvio
                Phylum = t_phylum, Family = t_family,
                Genus = t_genus, Species=t_species) %>% 
  mutate(Order= case_when(grepl("Pelagibacterales", Order)~ "Pelagibacterales", TRUE ~ Order)) #correcting some space character in the name


protein_taxonomy<- protein_taxonomy_refseq %>% filter(is.na(Class)) %>% select(gene_callers_id) %>% 
                    left_join(protein_taxonomy_nr, by = c("gene_callers_id")) %>% 
                    rbind(protein_taxonomy_refseq %>% filter(!is.na(Class)))

############################
#save gene callers IDs of the identified proteins for further annotation
############################
write(paste(row.names(protein_abund), collapse = ','), file ="./data/protein_GCIDs.txt")

############################
#summarize number of proteins
############################
protein_abund_per_sample <- protein_abund %>%
  reshape2::melt(variable.name = "Sample_ID", value.name = "Abund") %>% 
  filter(Abund>0) %>% 
  left_join(samples_df, by = "Sample_ID") %>% 
  group_by(Region,Station_ID,Fraction) %>% 
  summarize(N_prot=n()) 

#plot
protein_abund_per_sample %>% 
  ggplot(aes(x=Station_ID, y=N_prot))+
  geom_col()+
  facet_wrap(Fraction~., scales = "free")+
  xlab("Station ID")+ ylab("Number of proteins")+
  theme_EF+
  theme(axis.text.x=element_text(angle=90))
