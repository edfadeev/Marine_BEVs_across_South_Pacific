require(dplyr)
require(tibble)
require(stringr)
<<<<<<< HEAD
=======
require(ggplot2)
require(DEqMS)
require(vegan)

#load ggplot theme and colours
source("source/ggplot_parameters.R")
>>>>>>> e48d2ba4b277485019716045c8d9a68ab236f549

############################
#import metadata
############################
samples_df<- read.table("data/samples_meta.txt", header = TRUE) %>% 
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
         Razor.Peptides = case_when(grepl(";",Positions.in.Master.Proteins)~ 1, TRUE ~ 0),
         Modification = case_when(grepl("Oxidation",Modifications)~ 1, TRUE ~ 0)) %>% 
  group_by(gene_callers_id) %>% 
  summarize(Number.of.Peptides=n(), Number.of.PSMs=sum(Number.of.PSMs), Number.of.Unique.Peptide=sum(Unique.Peptides), Number.of.Razor.Peptide=sum(Razor.Peptides), Oxidation=sum(Modification)) %>% 
  filter(Number.of.Peptides>1, 
         Number.of.Unique.Peptide>0) %>% # remove proteins with less than 2 peptides and without unique peptides
  left_join(read.table("data/SO289-detected-proteins-calls.txt") %>%
              mutate(Length=nchar(aa_sequence),gene_callers_id=as.character(gene_callers_id)) %>% 
              select(gene_callers_id, contig, partial, Length, aa_sequence),
            by="gene_callers_id") %>% as.data.frame()

row.names(protein_metadata)<- protein_metadata$gene_callers_id

write.table(protein_metadata, "data/protein_metadata.txt", col.names = TRUE, row.names=TRUE, quote = FALSE)

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

<<<<<<< HEAD

write.table(protein_abund, "data/protein_abund.txt", col.names = TRUE, row.names=TRUE, quote = FALSE)

=======
>>>>>>> e48d2ba4b277485019716045c8d9a68ab236f549
############################
#Protein annotations
############################
#functional annotation using InterProScan
interpro_annotations<- lapply(c("InterPro", "Pfam","NCBIfam"
                                #,"SUPERFAMILY", "PANTHER","ProSitePatterns",
                                #"Phobius","SignalP_GRAM_NEGATIVE","SignalP_GRAM_POSITIVE" 
                                ),
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
#merge Interpro annotations into a table
interpro_annotations_df <- interpro_annotations%>% purrr::reduce(full_join, by = "gene_callers_id")           

#CAZymes characterization using DBcan
CAZYme_annotation_df <- read.table("data/SO289-detected_proteins_dbcan.txt", h=T) %>% 
  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
  dplyr::rename("CAZYme_acc"="annotation", "CAZYme_ann"="EC_number")

#protease characterization using MEROPS db
protease_annotation_df <- read.table("data/SO289-detected_proteins_merops.txt", h=T, sep ="\t") %>% 
  mutate(gene_callers_id =as.character(gene_callers_id)) %>% 
  dplyr::rename("MEROPS_acc"="accession", "MEROPS_ann"="e_value") %>% select(gene_callers_id,MEROPS_acc,MEROPS_ann)

#BLASTp against bacterial refseq
blastp_hits_df <- read.table("data/SO289-detected_proteins_blastp_refseq.txt", h=T, sep ="\t") %>% 
  mutate(gene_callers_id =as.character(gene_callers_id)) %>%  select(gene_callers_id,blastp_acc,blastp_ann, blastp_taxa, blastp_species)

#protein localization estimates
protein_localization <- read.table("data/SO289-detected_proteins_deeploc.txt", h=T, sep ="\t") %>% 
  mutate(gene_callers_id = as.character(gene_callers_id))

#merge all annotations
protein_annotations <- interpro_annotations_df %>% 
  merge(CAZYme_annotation_df, by = "gene_callers_id", all=TRUE) %>% 
  merge(protease_annotation_df, by = "gene_callers_id", all=TRUE) %>% 
  merge(blastp_hits_df %>% select(-c(blastp_taxa, blastp_species)), by = "gene_callers_id", all=TRUE) %>% 
  merge(protein_localization, by = "gene_callers_id", all=TRUE)


write.table(protein_annotations, "data/protein_annotations.txt", col.names = TRUE, row.names=FALSE, quote = FALSE, sep ="\t")

rm(protease_annotation_df, CAZYme_annotation_df, interpro_annotations_df, interpro_annotations, protein_localization)

############################
#Protein taxonomy
############################
#import refseq taxonomy
protein_taxonomy_refseq<- read.table("data/SO289-detected_proteins_kaiju_taxa_refseq.tsv", h=T, sep="\t") %>% 
  tidyr::separate(taxa, into = c("Domain","Phylum", "Order","Class", "Family","Genus", "Species"), sep = ";",convert=TRUE) %>% 
  mutate(gene_callers_id=as.character(gene_callers_id)) %>% mutate_if(is.character, str_trim)

#import nr taxonomy
protein_taxonomy_nr<- read.table("data/SO289-detected_proteins_kaiju_taxa_nr.tsv", h=T, sep="\t") %>% 
  tidyr::separate(taxa, into = c("Domain","Phylum", "Order","Class", "Family","Genus", "Species"), sep = ";", convert=TRUE) %>% 
  mutate_if(is.character, str_trim) %>% 
  mutate(gene_callers_id=as.character(gene_callers_id))

#complement missing Kaiju refseq taxonomy with results from nr
protein_taxonomy_kaiju<- protein_taxonomy_refseq %>% 
  filter(NCBI_taxID %in% c(0,1,2)) %>% select(gene_callers_id) %>% #select only still unclassified proteins
  left_join(protein_taxonomy_nr, by = c("gene_callers_id")) %>%
  rbind(protein_taxonomy_refseq %>% filter(!NCBI_taxID %in% c(0,1,2)))

#complement still missing values from blastp hits
protein_taxonomy_blastp<- blastp_hits_df %>% 
  filter(!grepl("MULTISPECIES:", blastp_ann)) %>% #exclude hits affiliated to multiple species
  tidyr::separate(blastp_taxa, into = c("Domain","Kingdom", "Phylum", "Class", "Order", "Family","Genus"), sep = ",", convert=TRUE) %>% 
  rename(Species=blastp_species) %>% 
  mutate(gene_callers_id=as.character(gene_callers_id)) %>% 
  mutate_if(is.character, str_trim) %>% select(-Kingdom)

protein_taxonomy_combined<- protein_taxonomy_kaiju %>% 
  filter(NCBI_taxID %in% c(0,1,2))%>% select(gene_callers_id) %>% #select only still unclassified proteins
  left_join(protein_taxonomy_blastp %>% select(-starts_with("blastp")), by = c("gene_callers_id")) %>% 
  rbind(protein_taxonomy_kaiju %>% filter(!NCBI_taxID %in% c(0,1,2)) %>% select(-c("NCBI_taxID")))

<<<<<<< HEAD
protein_taxonomy<- protein_taxonomy_combined %>% 
  filter(Phylum=="NA"|is.na(Phylum)) %>% select(gene_callers_id, Domain) %>% 
  left_join(protein_taxonomy_blastp, by =c("gene_callers_id", "Domain")) %>% 
  select(c("gene_callers_id", "Domain","Phylum", "Order","Class", "Family","Genus", "Species")) %>% 
  rbind(protein_taxonomy_combined %>% filter(Phylum!="NA")) %>% 
  mutate(Phylum=case_when(Phylum=="NA"~ paste0(Domain, "_uncl"), TRUE ~ Phylum))

write.table(protein_taxonomy, "data/protein_taxonomy.txt", col.names = TRUE, row.names=FALSE, quote = FALSE, sep ="\t")
=======
rm(protein_taxonomy_blastp, protein_taxonomy_kaiju, protein_taxonomy_nr, protein_taxonomy_refseq)
rm(protease_annotation_df, blastp_hits_df, CAZYme_annotation_df, interpro_annotations_df)
############################
#protein localization estimates
############################
protein_localization <- read.table("data/SO289-detected_proteins_deeploc.txt", h=T, sep ="\t") %>% 
  mutate(gene_callers_id = as.character(gene_callers_id))

############################
#Explore viral proteins in BEVs fraction
############################
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
  tidyr::spread(Fraction, N_prot)

############################
#summarize number of non viral proteins
############################
protein_abund_no_vir<-protein_abund[!row.names(protein_abund)%in%vir_gcids,]

protein_abund_per_sample <- protein_abund_no_vir %>%
  reshape2::melt(variable.name = "Sample_ID", value.name = "Abund") %>% 
  filter(Abund>0) %>% 
  left_join(samples_df, by = "Sample_ID") %>% 
  group_by(Region,Station_ID,Fraction) %>% 
  summarize(N_prot=n()) 
>>>>>>> e48d2ba4b277485019716045c8d9a68ab236f549

rm(blastp_hits_df, protein_taxonomy_blastp, protein_taxonomy_combined, protein_taxonomy_kaiju, protein_taxonomy_nr, protein_taxonomy_refseq)


<<<<<<< HEAD
#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
=======

############################
#Transform abundances
############################
#log2 transformation of protein abundances
protein_abund.log2 <- log2(as.matrix(protein_abund_no_vir))

#median normalization of the data
prot.dat.log2_norm = equalMedianNormalization(protein_abund.log2)
>>>>>>> e48d2ba4b277485019716045c8d9a68ab236f549
