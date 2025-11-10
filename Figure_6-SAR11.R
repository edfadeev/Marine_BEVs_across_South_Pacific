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

DEqMS_results<- read.table("data/DEqMS_results_regions.txt", h=T)%>% 
  mutate(gene_callers_id=as.character(gene_callers_id))



############################
# Explore enriched proteins of SAR11
############################
DEqMS.results_SAR11<- DEqMS_results %>% 
  left_join(protein_annotations %>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_taxonomy %>% unique(), by ="gene_callers_id") %>%
  filter(grepl("Pelagibacterales", Order)) %>% 
  left_join(protein_metadata, by ="gene_callers_id") %>% 
  mutate(Function=case_when(InterPro_ann=="-"~gsub("\\[.*|MULTISPECIES:","",blastp_ann), 
                            is.na(InterPro_ann) ~ blastp_ann, 
                            TRUE~InterPro_ann)) %>% 
  mutate(Function=case_when(grepl("Porin",Function, ignore.case = TRUE)~"Porin", TRUE~Function))
                                
                                
DEqMS.results_SAR11 %>% 
  group_by(Enr.frac, Function) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p) %>% 
  mutate(Cells=case_when(is.na(Cells)~0, TRUE~Cells),
         EVs=case_when(is.na(EVs)~0, TRUE~EVs)) %>% 
  #filter(EVs>Cells) %>% 
  View()


#characterize the porins
DEqMS.results_SAR11_ABC<- DEqMS.results_SAR11 %>% 
  filter(grepl("ABC|substrate-binding|solute-binding|Periplasmic binding protein",Function, ignore.case = TRUE))

DEqMS.results_SAR11_ABC_seqs<- DEqMS.results_SAR11_ABC%>% 
  select(gene_callers_id, aa_sequence) %>% unique() %>% 
  mutate(gene_callers_id=paste0("gcid_",gene_callers_id))

seqinr::write.fasta(as.list(DEqMS.results_SAR11_ABC_seqs$aa_sequence), names=DEqMS.results_SAR11_ABC_seqs$gene_callers_id, 
                    as.string=FALSE, file.out="data/SAR11_ABC.fa")



#import psiblast results
TCDB_annotation<- read.csv("data/TCDB_2_3_A.ann", h=F, sep="\t", fill=TRUE)
names(TCDB_annotation)<- c("sseqid", "TCID", "Annotation")

SAR11_ABC_annotations<- read.table("data/SAR11_ABC_psiblastp.out", col.names =c("gene_callers_id","sseqid","pident","length","evalue","bitscore")) %>% 
  left_join(TCDB_annotation, by ="sseqid", relationship = "many-to-many") %>% 
  filter(pident>30  & evalue<0.001 & bitscore >50) %>% 
  mutate(gene_callers_id=gsub("gcid_","", gene_callers_id)) %>% 
  group_by(gene_callers_id) %>% 
  slice_max(pident, n = 1)
  

#merge annotation
DEqMS.results_SAR11_ABC_ann<- DEqMS.results_SAR11_ABC %>% 
  left_join(SAR11_ABC_annotations, by = "gene_callers_id") %>% 
  mutate(Function=Annotation) %>%
  select(-c("sseqid","pident","length","evalue","bitscore", "TCID", "Annotation")) %>% 
  mutate(Function=case_when(is.na(Function)~gsub("\\[.*|MULTISPECIES:","",blastp_ann), 
                            TRUE~Function)) %>% 
  mutate(Function=case_when(grepl("Sugar",Function,ignore.case = TRUE)~ "Sugar binding protein",
                            grepl("General L-amino acid",Function,ignore.case = TRUE)~ "General L-amino acids binding protein",
                            grepl("Glutamate",Function,ignore.case = TRUE)~ "Glutamate/glutamine/aspartate/asparagine binding protein",
                            grepl("spermidine",Function,ignore.case = TRUE)~ "Spermidine/putrescine binding protein",
                            grepl("glycine",Function,ignore.case = TRUE)~ "Glycine betaine binding protein",
                            grepl("Iron",Function,ignore.case = TRUE)~ "Iron binding protein",
                            grepl("binding",Function,ignore.case = TRUE)~ "Putative periplasmic binding proteins",
                            TRUE ~Function))


DEqMS.results_SAR11_total_by_frac<- DEqMS.results_SAR11 %>% 
  filter(!grepl("ABC|substrate-binding|solute-binding|Periplasmic binding protein",Function, ignore.case = TRUE)) %>% 
  rbind(DEqMS.results_SAR11_ABC_ann) %>% 
  group_by(Enr.frac,  Function) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p) %>% 
  mutate(Cells=case_when(is.na(Cells)~0, TRUE~Cells),
         EVs=case_when(is.na(EVs)~0, TRUE~EVs)) %>% View()

#plot results
DEqMS.results_SAR11.p<- DEqMS.results_SAR11 %>% 
  filter(!grepl("ABC|substrate-binding|solute-binding|Periplasmic binding protein",Function, ignore.case = TRUE)) %>% 
  rbind(DEqMS.results_SAR11_ABC_ann) %>% 
  mutate(Function=case_when(grepl("Porin",Function, ignore.case = TRUE)~"Porin", 
                                grepl("NusA",Function, ignore.case = TRUE) ~"Transcription factor NusA",
                                grepl("SUF system",Function, ignore.case = TRUE) ~"FeS cluster assembly",
                                grepl("PBP",Function, ignore.case = TRUE) ~"Phosphate binding protein PstS",
                                TRUE~Function)) %>% 
  group_by(Enr.group, Enr.frac, Function) %>% 
  summarize(log2fold_mean = mean(logFC), log2fold_median = median(logFC), log2fold_min = min(logFC), log2fold_max = max(logFC), log2fold_se = se(logFC), count=n())

total_SAR11_prot<- DEqMS.results_SAR11.p %>% group_by(Enr.frac,Enr.group) %>% summarize(Total_p=sum(count))

DEqMS.results_SAR11.p %>%
  filter(Function %in% c(DEqMS.results_SAR11.p %>% 
                               filter(count>1 & Enr.frac=="EVs") %>% 
                               pull(Function))) %>% 
  left_join(total_SAR11_prot, by =c("Enr.frac", "Enr.group")) %>% 
  mutate(Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN")),
         Prop=count/Total_p) %>% 
  ggplot(aes(y=Function , x=log2fold_mean, fill = Enr.group, label = count))+ 
  geom_point(aes(size = Prop), shape =21)+
  geom_text(size = 5, nudge_y = -0.2)+
  geom_errorbar(aes(xmin = log2fold_mean-log2fold_se, xmax = log2fold_mean +log2fold_se), 
                width = 0.2) + 
  xlim(-8,8)+
  facet_grid(~Enr.group)+  
  geom_vline(aes(xintercept=0), linetype="dashed")+
  scale_size_continuous(range = c(1, 20))+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  theme_EF+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=90))


#save the plot
ggsave("./Figures/SAR11_proteins_regions.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)
