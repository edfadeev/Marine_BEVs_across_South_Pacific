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

DEqMS_results<- read.table("data/DEqMS_results_regions.txt", h=T)


############################
# Explore enriched proteins of Cyanobacteria
############################
DEqMS.results_cyano<- DEqMS_results %>% 
  left_join(protein_annotations %>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_taxonomy %>% unique(), by ="gene_callers_id") %>%
  filter(grepl("Synechococcales", Order)) %>% 
  left_join(protein_metadata, by ="gene_callers_id")


DEqMS.results_cyano %>% 
  mutate(InterPro_ann=case_when(InterPro_ann=="-"~NCBIfam_ann, is.na(InterPro_ann) ~ blastp_ann, TRUE~InterPro_ann)) %>% 
  mutate(InterPro_ann=gsub("\\[.*|MULTISPECIES:","",InterPro_ann)) %>% View()
group_by(Enr.frac, InterPro_ann) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p) %>% 
  mutate(Cells=case_when(is.na(Cells)~0, TRUE~Cells),
         EVs=case_when(is.na(EVs)~0, TRUE~EVs)) %>% 
  #filter(EVs>Cells) %>% 
  View()


#characterize the porins
DEqMS.results_cyano_Porins<- DEqMS.results_cyano %>% 
  mutate(InterPro_ann=case_when(InterPro_ann=="-"~NCBIfam_ann, is.na(InterPro_ann) ~ blastp_ann, TRUE~InterPro_ann)) %>% 
  mutate(InterPro_ann=gsub("\\[.*|MULTISPECIES:","",InterPro_ann)) %>% 
  filter(grepl("Porin",InterPro_ann, ignore.case = TRUE))

DEqMS.results_cyano_Porins_seqs<- DEqMS.results_cyano_Porins%>% 
  select(gene_callers_id, aa_sequence) %>% unique() %>% 
  mutate(gene_callers_id=paste0("gcid_",gene_callers_id))

seqinr::write.fasta(as.list(DEqMS.results_cyano_Porins_seqs$aa_sequence), names=DEqMS.results_cyano_Porins_seqs$gene_callers_id, 
                    as.string=FALSE, file.out="data/Cyano_porins.fa")


#transfer the files to LISC and run psiblast against TCDB
# module load ncbiblast
# makeblastdb -in $WORKDIR/11_PROTEIN/TCDB_class/TCDB_1_B.fasta -out $WORKDIR/11_PROTEIN/TCDB_class/TCDB_1_B -dbtype prot
# grep ">" TCDB_1_B.fasta | sed 's/>//g' - | sed 's/ /'$'\t''/' - | sed 's/ /'$'\t''/' - > TCDB_1_B.ann
# psiblast -db $WORKDIR/11_PROTEIN/TCDB_class/TCDB_1_B -query $WORKDIR/11_PROTEIN/TCDB_class/Cyano_porins.fa -outfmt "6 qseqid sseqid pident length evalue bitscore" -out Cyano_porins_psiblastp.out -num_threads 4


#import psiblast results
blastp.out<- read.table("data/Cyano_porins_psiblastp.out", col.names =c("gene_callers_id","sseqid","pident","length","evalue","bitscore")) %>% 
  filter(pident>30  & evalue<0.001 & bitscore >50 ) %>% 
  mutate(gene_callers_id=gsub("gcid_","", gene_callers_id))
TCDB_annotation<- read.table("data/TCDB_1_B.ann", h=F, sep="\t")
names(TCDB_annotation)<- c("sseqid", "TCID", "Annotation")

#best hit by identitiy
Cyano_porin_ann<- blastp.out %>% 
  left_join(TCDB_annotation, by ="sseqid") %>%
  group_by(gene_callers_id) %>% 
  slice_max(pident, n = 1) 

#merge annotation
DEqMS.results_cyano_Porins<- DEqMS.results_cyano_Porins %>% 
  left_join(Cyano_porin_ann, by = "gene_callers_id") %>% 
  mutate(InterPro_ann=paste(TCID,Annotation, sep="_")) %>%
  mutate(InterPro_ann=case_when(InterPro_ann=="NA_NA"~"Unclassified porin", TRUE~InterPro_ann)) %>% 
  select(-c("sseqid","pident","length","evalue","bitscore", "TCID", "Annotation"))

#merge and summarize
DEqMS.results_cyano_total_by_frac<- DEqMS.results_cyano %>% 
  mutate(InterPro_ann=case_when(InterPro_ann=="-"~NCBIfam_ann, is.na(InterPro_ann) ~ blastp_ann, TRUE~InterPro_ann)) %>% 
  mutate(InterPro_ann=gsub("\\[.*|MULTISPECIES:","",InterPro_ann)) %>% 
  filter(!grepl("Porin",InterPro_ann, ignore.case = TRUE)) %>% 
  rbind(DEqMS.results_cyano_Porins) %>% 
  group_by(Enr.frac, InterPro_ann) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p) %>% 
  mutate(Cells=case_when(is.na(Cells)~0, TRUE~Cells),
         EVs=case_when(is.na(EVs)~0, TRUE~EVs)) 

#plot results
DEqMS.results_cyno.p<- DEqMS.results_cyano %>% 
  mutate(InterPro_ann=case_when(InterPro_ann=="-"~Pfam_ann, is.na(InterPro_ann) ~ blastp_ann, TRUE~InterPro_ann)) %>% 
  mutate(InterPro_ann=gsub("\\[.*|MULTISPECIES:","",InterPro_ann)) %>% 
  filter(!grepl("Porin",InterPro_ann, ignore.case = TRUE)) %>% 
  rbind(DEqMS.results_cyano_Porins) %>% 
  #filter(InterPro_ann %in% c(DEqMS.results_cyano_total_by_frac %>% filter(EVs>Cells & EVs>4) %>% pull(InterPro_ann) )) %>% 
  group_by(Enr.group, Enr.frac, InterPro_ann) %>% 
  summarize(log2fold_mean = mean(logFC), log2fold_median = median(logFC), log2fold_min = min(logFC), log2fold_max = max(logFC), log2fold_se = se(logFC), count=n())

total_cyano_prot<- DEqMS.results_cyno.p %>% group_by(Enr.frac,Enr.group) %>% summarize(Total_p=sum(count))

DEqMS.results_cyno.p %>%
  filter(InterPro_ann %in% c(DEqMS.results_cyno.p %>% 
                               filter(count>2 & Enr.frac=="EVs") %>% 
                               pull(InterPro_ann))) %>% 
  mutate(Function=case_when(InterPro_ann=="1.B.23.1.2_Hypothetical protein slr0042 - Synechocystis sp. (strain PCC 6803)." ~"Iron uptake porin",
                            InterPro_ann=="1.B.23.1.9_Carbohydrate-selective porin OprB OS=Fischerella sp. JSC-11 GN=FJSC11DRAFT_0273 PE=4 SV=1" ~"Carbohydrate-selective porin",
                            InterPro_ann=="1.B.23.1.1_SOMA - Synechococcus sp. (strain PCC 6301) (Anacystis nidulans)."~ "SomA porin",
                            InterPro_ann=="Flavodoxin/nitric oxide synthase" ~ "Flavoprotein",
                            InterPro_ann=="Phycobilisome, alpha/beta subunit"~"Phycobilisome",
                            InterPro_ann=="Ferritin/DPS protein domain"~"Ferredoxin", TRUE~ InterPro_ann)) %>% 
  filter(!grepl("hypothetical", Function)) %>% 
  left_join(total_cyano_prot, by =c("Enr.frac", "Enr.group")) %>% 
  mutate(Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN")),
         Prop=count/Total_p) %>% 
  ggplot(aes(y=Function , x=log2fold_mean, fill = Enr.group, label = count))+ 
  geom_point(aes(size = Prop), shape =21)+
  geom_text(size = 5, nudge_y = -0.2)+
  geom_errorbar(aes(xmin = log2fold_mean-log2fold_se, xmax = log2fold_mean +log2fold_se), 
                width = 0.2) + 
  xlim(-10,10)+
  facet_grid(~Enr.group)+  
  geom_vline(aes(xintercept=0), linetype="dashed")+
  scale_size_continuous(range = c(1, 20))+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  theme_EF+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=90))


#save the plot
ggsave("./Figures/cyano_proteins_regions.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)
