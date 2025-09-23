require(dplyr)
require(tibble)
require(stringr)
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
protein_metadata<- read.table("data/SO289_clust99_RC_PeptideGroups.txt", sep ="\t", header = TRUE) %>% 
  dplyr::rename(gene_callers_id=Master.Protein.Accessions) %>%
  tidyr::separate_rows(gene_callers_id, sep= "; ") %>%
  #calculate total of unique peptide abundance by protein
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>%
  rename_at(all_of(samples_df$File.ID), ~ samples_df$Sample_ID)%>%
  filter(Number.of.PSMs>1, gene_callers_id!="", PSM.Ambiguity!="Rejected") %>% #include only peptides with at least two PSMs
  mutate(gene_callers_id=as.character(gene_callers_id),
         Unique.Peptides = case_when(grepl(";",Positions.in.Master.Proteins)~ 0, TRUE ~ 1),
         Razor.Peptides = case_when(grepl(";",Positions.in.Master.Proteins)~ 1, TRUE ~ 0)) %>% 
  group_by(gene_callers_id) %>% 
  summarize(Number.of.Peptides=n(), Number.of.PSMs=sum(Number.of.PSMs), Number.of.Unique.Peptide=sum(Unique.Peptides), Number.of.Razor.Peptide=sum(Razor.Peptides)) %>% 
  filter(#Number.of.Peptides>1, 
         Number.of.Unique.Peptide> 0) %>% # remove proteins with less than 2 peptides and without unique peptides
  as.data.frame()

row.names(protein_metadata)<- protein_metadata$gene_callers_id

#summarize abundance of peptides by protein
protein_abund<- read.table("data/SO289_clust99_RC_PeptideGroups.txt", sep ="\t", header = TRUE) %>% 
  dplyr::rename(gene_callers_id=Master.Protein.Accessions) %>%
  tidyr::separate_rows(gene_callers_id, sep= "; ") %>%
  #calculate total of unique peptide abundance by protein
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>%
  rename_at(all_of(samples_df$File.ID), ~ samples_df$Sample_ID)%>%
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
#Protein annotations
############################
interpro_annotations<- lapply(c("InterPro", "NCBIfam","Pfam","SUPERFAMILY","PANTHER","ProSitePatterns",
                                "TIGRFAM","Phobius","SignalP_GRAM_NEGATIVE","SignalP_GRAM_POSITIVE" ),
                              function(src){
                                annotation<- paste0(src,"_ann")
                                accession<- paste0(src,"_acc")
                                annotation_df <- read.table(paste0("data/SO289-detected_proteins_",src,".tsv"),quote = "",  sep="\t", h=T) %>% 
                                  dplyr::rename(!!sym(annotation):= function.,
                                                !!sym(accession):= accession) %>% 
                                  select(-c("source")) %>% 
                                  group_by(gene_callers_id)%>%
                                  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
                                  summarise_all(funs(paste(unique(.), collapse='|')))
                                
                                return(annotation_df)
                              })  

interpro_annotations_df <- interpro_annotations%>% purrr::reduce(full_join, by = "gene_callers_id")           

blastp_hits_df <- read.table("data/SO289-detected_proteins_blastp_refseq.txt", h=T, sep ="\t") %>% 
  mutate(gene_callers_id =as.character(gene_callers_id))


CAZYme_annotation_df <- read.table("data/SO289-detected_proteins_dbcan.txt", h=T) %>% 
  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
  dplyr::rename(CAZYme=annotation)

protease_annotation_df <- read.table("data/SO289-detected_proteins_merops.txt", h=T, sep ="\t") %>% 
  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
  dplyr::rename(MEROPS=function.)


protein_annotations <- interpro_annotations_df %>% 
  left_join(CAZYme_annotation_df %>%  select(gene_callers_id,CAZYme), by = "gene_callers_id") %>% 
  left_join(protease_annotation_df %>%  select(gene_callers_id,MEROPS), by = "gene_callers_id") %>% 
  left_join(blastp_hits_df %>%  select(gene_callers_id,blastp_acc,blastp_ann), by = "gene_callers_id")


############################
#Protein taxonomy
############################
#import refseq taxonomy
protein_taxonomy_refseq<- read.table("data/SO289-detected_proteins_kaiju_taxa_refseq.tsv", h=T, sep="\t") %>% 
  tidyr::separate(taxa, into = c("Domain","Phylum", "Order","Class", "Family","Genus", "Species"), sep = ";") %>% 
  mutate(gene_callers_id=as.character(gene_callers_id)) %>% mutate_if(is.character, str_trim)

#import nr taxonomy
protein_taxonomy_nr<- read.table("data/SO289-detected_proteins_kaiju_taxa_nr.tsv", h=T, sep="\t") %>% 
  tidyr::separate(taxa, into = c("Domain","Phylum", "Order","Class", "Family","Genus", "Species"), sep = ";") %>% 
  mutate(Order=gsub(" ","",Order),Class=gsub(" ","",Class) ) %>% 
  mutate(Order = case_when(Order=="NA" ~ paste0(Class, "_uncl"),
                          TRUE~Order)) %>% 
  mutate_if(is.character, str_trim) %>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

#complement missing refseq taxonomy with results from nr
protein_taxonomy_kaiju<- protein_taxonomy_refseq %>% filter(NCBI_taxID==0) %>% 
  select(gene_callers_id) %>% 
  left_join(protein_taxonomy_nr, by = c("gene_callers_id")) %>% 
  rbind(protein_taxonomy_refseq %>% filter(NCBI_taxID!=0))

#complement still missing values from blastp hits
protein_taxonomy_blastp<- blastp_hits_df %>% 
tidyr::separate(blastp_taxa, into = c("Domain","Phylum", "Order","Class", "Family","Genus", "Species"), sep = ",") %>% 
  mutate(gene_callers_id=as.character(gene_callers_id)) %>% 
  mutate_if(is.character, str_trim)

protein_taxonomy<- protein_taxonomy_kaiju %>% filter(is.na(NCBI_taxID))%>% select(gene_callers_id) %>% 
  left_join(protein_taxonomy_blastp %>% select(-starts_with("blastp")), by = c("gene_callers_id")) %>% 
  rbind(protein_taxonomy_kaiju %>% filter(!is.na(NCBI_taxID)) %>% select(-c("NCBI_taxID")))

rm(protein_taxonomy_blastp, protein_taxonomy_kaiju, protein_taxonomy_nr, protein_taxonomy_refseq)
rm(protease_annotation_df, blastp_hits_df, CAZYme_annotation_df, interpro_annotations_df)
############################
#protein localization estimates
############################
protein_localization <- read.table("data/SO289-detected_proteins_deeploc.txt", h=T, sep ="\t") %>% 
  mutate(gene_callers_id = as.character(gene_callers_id)) %>% 
  filter(gene_callers_id %in% row.names(prot.dat.log2_norm))

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


############################
#Check differences in protein composition between fraction
############################
prot.dat.MDS<- prot.dat.log2_norm
prot.dat.MDS[is.na(prot.dat.MDS)]<- 0
prot.dat.MDS<- prot.dat.MDS[!rowSums(prot.dat.MDS)==0,]

# Running NMDS in vegan (metaMDS)
prot_NMDS <-  metaMDS(t(prot.dat.MDS),
                      distance = "euclidean",
                      k = 2,
                      maxit = 999, 
                      trymax = 999,
                      wascores = FALSE,
                      autotransform = FALSE,
                      tidy= "sites",
                      na.rm = TRUE)

# Perform K-means clustering with 4 clusters
kmeans_result <- kmeans(t(prot.dat.MDS), centers=2, nstart = 25)

#plot
prot_NMDS.scores <- as.data.frame(scores(prot_NMDS)) %>% 
  tibble::rownames_to_column(var ="Sample_ID") %>% 
  left_join(samples_df, 
            by = "Sample_ID", copy = TRUE) %>% 
  mutate(Cluster=kmeans_result$cluster)

prot_NMDS.scores %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, shape = Fraction, 
             colour = as.factor(Cluster), label = Sample_ID))+
  geom_point(size =5, colour = "black")+
  geom_point(size =4)+
  geom_text(nudge_y = -0.05, size =4)+
  scale_color_manual(values = c("#009E73", "#F0E442", "#0072B2", 
                                "#D55E00"))+
  theme_bw()

#test whether the differences between the runs are significant
prot_distmat <- 
  vegdist(t(prot.dat.MDS), method = "euclidean",na.rm = TRUE)

df <- samples_df %>% 
  select(Region,Fraction, Sample_ID) %>% 
  mutate(Cluster=kmeans_result$cluster)

adonis_all <- adonis2(prot_distmat ~ Fraction+Cluster, df,
                      permutations=999)

adonis_all

