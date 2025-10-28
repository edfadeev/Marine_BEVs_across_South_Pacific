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

############################
#Protein enrichment analysis between fractions in each region
############################
#carry out enrichment tests based on regions
enrichment_tests_list <- lapply(c("WEST","GYRE", "TRAN"), function(x) {
  
  samples_meta_sub<- samples_df %>% 
    mutate(Enr.group = factor(Region, levels =c("WEST","GYRE", "TRAN"))) %>% 
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
  
  prot.dat.log2_norm.filter <- prot.dat.log2_norm[,c(EV_sample_IDs,Cell_sample_IDs)]
  
  # Filter proteins that were observed in at least two samples in each fraction
  prot.dat.log2_norm.filter <- prot.dat.log2_norm.filter[rowSums(!is.na(prot.dat.log2_norm.filter[, c(EV_sample_IDs)]))>1 &
                                                           rowSums(!is.na(prot.dat.log2_norm.filter[, c(Cell_sample_IDs)]))>1,]
  
  #run linear model
  fit1<- lmFit(prot.dat.log2_norm.filter, design)
  fit2 <- contrasts.fit(fit1,contrasts = contrast)
  fit3<- eBayes(fit2)
  
  #correct bias of variance estimate based on number of PSMs per protein
  fit3$count <- protein_metadata[rownames(fit3$coefficients), "Number.of.PSMs"]

  fit4 = spectraCounteBayes(fit3)
  
  #results
  DEqMS.results <- outputResult(fit4,coef_col = 1) %>% 
                 dplyr::rename("gene_callers_id"="gene") %>% 
                 mutate(Enr.frac = case_when(logFC>1 & sca.adj.pval<0.1  
                                             ~ "EVs",
                                             -1> logFC & sca.adj.pval<0.1  
                                             ~ "Cells", TRUE ~ "Not.enr"),
                        log.sca.pval = -log10(sca.P.Value),
                        Enr.group =x) 
  
  return(DEqMS.results)
})

############################
#enrichment results 
############################
DEqMS_results<- bind_rows(enrichment_tests_list) %>% 
  filter(Enr.frac!="Not.enr")

#total of different proteins
DEqMS_results %>%  
  select(Enr.frac, gene_callers_id) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())

#total of different orders
DEqMS_results %>% select(Enr.frac, gene_callers_id) %>% unique() %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  select(Enr.frac, Class, Order) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_orders = n())

#proteins per order
DEqMS_results %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  select(Enr.frac, Phylum, Class, Order, gene_callers_id) %>% unique() %>% 
  group_by(Enr.frac, Phylum, Class, Order) %>% 
  summarize(N_p = n()) %>% 
  tidyr::spread(Enr.frac, N_p) %>% View()

#proteins per order per region
DEqMS_results %>% 
  left_join(protein_taxonomy %>% unique(), by = "gene_callers_id") %>% 
  select(Enr.group, Enr.frac, Phylum, Class, Order, gene_callers_id) %>% unique() %>% 
  group_by(Enr.group, Enr.frac, Phylum, Class, Order) %>% 
  summarize(N_p = n()) %>% 
  tidyr::spread(Enr.frac, N_p) %>% View()

#total of proteins epr region
DEqMS_results_totals<- DEqMS_results %>% 
  select(Enr.group, Enr.frac, gene_callers_id) %>% 
  group_by(Enr.group, Enr.frac) %>% 
  summarize(Total_p = n())

############################
#taxonomy overview of enriched proteins
############################
#plot taxa of the enriched proteins
DEqMS_results_Class_tax<- DEqMS_results %>% 
  left_join(protein_taxonomy , by = "gene_callers_id") %>% 
  group_by(Enr.group, Enr.frac, Class, Order) %>% 
  summarize(N_p = n()) 

DEqMS_results_Class_totals<- DEqMS_results_Class_tax %>% 
  group_by(Enr.group, Enr.frac) %>% summarize(Total_p = sum(N_p)) 

DEqMS_results_Class_tax %>% 
  left_join(DEqMS_results_Class_totals, by = c("Enr.group", "Enr.frac")) %>% 
  mutate(Prop=100*N_p/Total_p, Order=gsub("^ ","",Order),Class=gsub("^ ","",Class) ) %>% 
  mutate(Order=gsub("Candidatus |Candidatus","",Order),
          Taxa=case_when(Prop<1 ~"Other taxa < 1%", 
                        is.na(Order) | Order=="NA" |Order=="NA_uncl" | Order=="Unknown" ~"unclassified",
                        TRUE~Order)) %>%  
  mutate(Taxa = factor(Taxa, levels=c( "Caulobacterales", 
                                       "Hyphomicrobiales" , 
                                       "Hyphomonadales","Kordiimonadales", "Minwuiales","Parvularculales","Pelagibacterales" ,"Puniceispirillales",
                                       "Rhodobacterales" , "Rhodospirillales" , "Sphingomonadales" , "Alphaproteobacteria_uncl", "Burkholderiales", 
                                       "Chitinophagales", "Cytophagales" , "Synechococcales", "Flavobacteriales" ,
                                       "Alteromonadales" , "Cellvibrionales" , "Chromatiales", "Kangiellales", "Lysobacterales" ,"Moraxellales",
                                       "Oceanospirillales" , "Pseudomonadales" , "Sphingobacteriales", "Gammaproteobacteria_uncl", 
                                       "Other taxa < 1%", "Viral", "unclassified")),
         Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN","UP")),
         Enr.frac=ifelse(Enr.frac=="EVs","BEVs",Enr.frac)) %>%  
  ggplot(aes(x=Enr.frac, y= Prop, fill = Taxa))+
  geom_col()+
  scale_fill_manual(values = tol21rainbow)+
  facet_grid(cols=vars(Enr.group),scales="free_x",space="free_x",switch="x") +
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/prot_enr_Tax_regions.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)


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
# psiblast -db $WORKDIR/11_PROTEIN/TCDB_class/TCDB_1_B -query $WORKDIR/11_PROTEIN/TCDB_class/Cyano_porins.fa -outfmt "6 qseqid sseqid pident length evalue bitscore" -out Cyano_porins_psiblastp.out -num_threads 4
# grep ">" TCDB_1_B.fasta | sed 's/>//g' - | sed 's/ /'$'\t''/' - | sed 's/ /'$'\t''/' - > TCDB_1_B.ann



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
  mutate(InterPro_ann=case_when(InterPro_ann=="NA_NA"~"Unknown porin", TRUE~InterPro_ann)) %>% 
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
  mutate(InterPro_ann=case_when(InterPro_ann=="-"~NCBIfam_ann, is.na(InterPro_ann) ~ blastp_ann, TRUE~InterPro_ann)) %>% 
  mutate(InterPro_ann=gsub("\\[.*|MULTISPECIES:","",InterPro_ann)) %>% 
  filter(!grepl("Porin",InterPro_ann, ignore.case = TRUE)) %>% 
  rbind(DEqMS.results_cyano_Porins) %>% 
  filter(InterPro_ann %in% c(DEqMS.results_cyano_total_by_frac %>% filter(EVs>Cells & EVs>4) %>% pull(InterPro_ann) )) %>% group_by(Enr.group, Enr.frac, InterPro_ann) %>% 
  summarize(log2fold_mean = mean(logFC), log2fold_median = median(logFC), log2fold_min = min(logFC), log2fold_max = max(logFC), log2fold_se = se(logFC), count=n())


DEqMS.results_cyno.p %>% 
  ggplot(aes(y=log2fold_mean , x=InterPro_ann, colour = Enr.group, label = count))+ 
  geom_text(size = 5, position = position_dodge(width = 0.5))+
  geom_point(size = 5, position = position_dodge(width = 1))+
  scale_size_continuous(range = c(1, 20))+
  geom_errorbar(aes(ymin = log2fold_mean-log2fold_se, ymax = log2fold_mean +log2fold_se), 
                width = 0.2, position = position_dodge(width = 1)) + 
  ylab("log2 foldchange")+
  scale_color_manual(values = c("#009E73", "#F0E442", "#0072B2", 
                                "#D55E00","#E32356"))+
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 30)) +
  facet_grid(.~Enr.group)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme_EF+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=90))




DEqMS.results_cyno %>% 
  group_by(Enr.frac, Enr.group) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p)

DEqMS.results_cyno.p<- DEqMS.results_cyno %>% 
  filter(grepl("Ferritin|Flavodoxin|Peroxiredoxin|Superoxide|oxidoreductase", InterPro_ann, ignore.case = TRUE)) %>% 
  mutate(InterPro_ann=case_when(grepl("Ferritin", InterPro_ann)~"Ferritins",
                                grepl("Flavodoxin", InterPro_ann)~"Flavodoxins",
                                grepl("Peroxiredoxin", InterPro_ann)~"Peroxiredoxins",
                                grepl("Superoxide", InterPro_ann)~"Superoxide dismutase",
                                grepl("oxidoreductase", InterPro_ann, ignore.case = TRUE)~"Oxidoreductase"),
         Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN"))) %>% 
  group_by(Enr.group, Enr.frac, InterPro_ann) %>% 
  summarize(log2fold_mean = mean(logFC), log2fold_median = median(logFC), log2fold_min = min(logFC), log2fold_max = max(logFC), log2fold_se = se(logFC), count=n())

DEqMS.results_cyno.totals<- DEqMS.results_cyno %>% 
  group_by(Enr.group, Enr.frac) %>% 
  summarize(Total=n())
  
  



#save the plot
ggsave("./Figures/cyano_oxids_regions.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)




############################
# Explore enriched proteins of SAR11
############################
DEqMS.results_SAR11<- DEqMS_results %>% 
  left_join(protein_annotations %>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_taxonomy %>% unique(), by ="gene_callers_id") %>%
  filter(grepl("Pelagibacterales", Order)) %>% 
  left_join(protein_metadata, by ="gene_callers_id")



DEqMS.results_SAR11 %>% filter(grepl("Periplasmic binding", SUPERFAMILY_ann)) %>% 
  select(Enr.frac, gene_callers_id, PANTHER_ann) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(p=n())

DEqMS.results_SAR11 %>% filter(grepl("Periplasmic binding", SUPERFAMILY_ann)) %>% 
  select(Enr.frac, gene_callers_id, PANTHER_ann) %>% unique() %>% 
  group_by(Enr.frac, PANTHER_ann) %>% 
  summarize(p=n()) %>% View()

DEqMS.results_SAR11 %>% filter(grepl("Periplasmic binding", SUPERFAMILY_ann)) %>% 
  select(Enr.frac, Enr.group,gene_callers_id, PANTHER_ann) %>% unique() %>% 
  group_by(Enr.frac, Enr.group, PANTHER_ann) %>% 
  summarize(p=n()) %>% View()


DEqMS.results_SAR11 %>% 
  group_by(Enr.frac, Enr.group) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p)



#SusCD proteins
SusCD_acc <- c("PF14322", "PF12771", "PF07980", "PF12741",
               "TIGR04056", "TIGR04057", "TIGR04056|TIGR04057", "TIGR04057|TIGR04056")

DEqMS.results_Flavobacteriales<- DEqMS_results %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(protein_taxonomy%>% unique(), by = "gene_callers_id") %>% 
  filter(Order =="Flavobacteriales")  %>% 
  left_join(protein_metadata%>% unique(), by ="gene_callers_id") %>%
  left_join(protein_swissprot%>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_uniref%>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_localization%>% unique(),  by="gene_callers_id") %>% 
  left_join(protein_annotations %>% unique(), by = "gene_callers_id") %>% 
  mutate(Sus=case_when(grepl("Susd", Pfam_ann, ignore.case =TRUE)|
                         grepl("SusD", InterPro_ann, ignore.case =TRUE)|
                         grepl("SusD", blastp_ann, ignore.case =TRUE)|
                         grepl("Q8A1G2", blastp_acc.x) ~"SusD",
                       NCBIfam_acc %in% SusCD_acc |grepl("SusC", InterPro_ann, ignore.case =TRUE)|
                         grepl("SUSC", PANTHER_ann, ignore.case =TRUE)|
                         grepl("SusC", blastp_ann, ignore.case =TRUE)|
                         grepl("T2KPJ3|T2KMI3|T2KM18|Q8A1G1", blastp_acc.x)~"SusC"))


DEqMS.results_Flavobacteriales %>% filter(!is.na(Sus)) %>% 
  select(Enr.frac, Sus, gene_callers_id) %>% unique() %>%  
  group_by(Enr.frac, Sus) %>% 
  summarize(p=n()) %>% 
spread(Enr.frac, p)

DEqMS.results_Flavobacteriales %>% filter(!is.na(Sus)) %>% 
  select(Enr.frac, Sus, gene_callers_id) %>% unique() %>%  
  group_by(Enr.frac, Sus) %>% 
  summarize(p=n()) %>% View()


DEqMS.results_Flavobacteriales %>% filter(!is.na(Sus)) %>% 
  group_by(Enr.frac, Enr.group, Sus) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p)

DEqMS.results_Flavobacteriales %>% 
  group_by(Enr.frac, Enr.group, Sus) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p)


DEqMS.results_TonB<- DEqMS_results %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(protein_taxonomy%>% unique(), by = "gene_callers_id") %>% 
  #filter(Class== "Gammaproteobacteria")  %>% 
  left_join(protein_metadata%>% unique(), by ="gene_callers_id") %>%
  left_join(protein_swissprot%>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_uniref%>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_localization%>% unique(),  by="gene_callers_id") %>% 
  left_join(protein_annotations %>% unique(), by = "gene_callers_id") %>% 
  mutate(TonB=case_when(grepl("Susd", Pfam_ann, ignore.case =TRUE)|
                        grepl("SusD", InterPro_ann, ignore.case =TRUE)|
                        grepl("SusD", blastp_ann, ignore.case =TRUE)|
                        grepl("Q8A1G2", blastp_acc.x)|
                        NCBIfam_acc %in% SusCD_acc |grepl("SusC", InterPro_ann, ignore.case =TRUE)|
                        grepl("SUSC", PANTHER_ann, ignore.case =TRUE)|
                        grepl("SusC", blastp_ann, ignore.case =TRUE)|
                        grepl("T2KPJ3|T2KMI3|T2KM18|Q8A1G1", blastp_acc.x)|
                        NCBIfam_acc %in% TonBs_acc |
                        TIGRFAM_acc %in% TonBs_acc |
                        grepl("TonB", InterPro_ann)| 
                        grepl("TonB", ProSitePatterns_ann)|
                        grepl("TONB", PANTHER_ann)|
                        grepl("TONB", blastp_ann)~ "TonB-dependent transport systems")) 




DEqMS.results_TonB %>% 
  select(Enr.frac, Enr.group, gene_callers_id, TonB) %>% unique() %>% 
  group_by(Enr.frac, Enr.group, TonB) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p) %>% 
  mutate(Prop_Cells=Cells/sum(Cells, na.rm = TRUE),
         Prop_EVs=EVs/sum(EVs, na.rm = TRUE)) 






DEqMS.results_Alteromonas %>% 
  group_by(Enr.frac, deeploc) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p) %>% 
  mutate(Prop_Cells=Cells/sum(Cells, na.rm = TRUE),
         Prop_EVs=EVs/sum(EVs, na.rm = TRUE)) 

DEqMS.results_Rhodobacterales<- DEqMS_results %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(protein_taxonomy%>% unique(), by = "gene_callers_id") %>% 
  filter(Order =="Rhodobacterales")  %>% 
  left_join(protein_metadata%>% unique(), by ="gene_callers_id") %>%
  left_join(protein_swissprot%>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_uniref%>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_localization%>% unique(),  by="gene_callers_id") %>% 
  left_join(protein_annotations %>% unique(), by = "gene_callers_id") 


DEqMS.results_Rhodobacterales %>% 
  filter(is.na(InterPro_ann)) %>% #View()
  select(Enr.frac, Enr.group, gene_callers_id) %>% unique() %>% 
  group_by(Enr.group, Enr.frac) %>% 
  summarize(p=n())

DEqMS.results_Rhodobacterales %>% 
  group_by(Enr.frac, deeploc) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p) %>% 
  mutate(Prop_Cells=Cells/sum(Cells, na.rm = TRUE),
         Prop_EVs=EVs/sum(EVs, na.rm = TRUE)) 

############################
#iron related porins in cyanobacteria
############################
#filter cyanobacterial proteins
DEqMS_results_cyano<- DEqMS_results %>%
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  filter(grepl("Cyanophyceae", Class)) 

#summarize number of protein
DEqMS_results_cyano%>% 
  group_by(Enr.group, Enr.frac, Class, Order) %>% 
  summarize(N_p = n()) 


DEqMS_results_cyano%>%
  left_join(protein_annotations, by = "gene_callers_id") %>% 
  mutate(Fe_porin=case_when(grepl("porin", InterPro_ann, ignore.case = TRUE)|
                        grepl("iron uptake porin", NCBIfam_ann, ignore.case = TRUE)|
                        grepl("porin", Pfam_ann, ignore.case = TRUE)|
                        grepl("iron uptake porin", blastp_ann, ignore.case = TRUE) ~ "Yes",
                        TRUE ~ "No")) %>% 
  group_by(Enr.group, Enr.frac, Fe_porin) %>% 
  summarize(Total_p = n()) %>% 
  tidyr::spread(Enr.frac, Total_p)


############################
#TonB-related proteins
############################
#SusCD proteins
SusCD_acc <- c("PF14322", "PF12771", "PF07980", "PF12741",
               "TIGR04056", "TIGR04057", "TIGR04056|TIGR04057", "TIGR04057|TIGR04056")
#TonB proteins
TonBs_acc<- c("TIGR01352", #TonB family C-terminal domain
              "TIGR01776",
              "TIGR01778", #copper
              "TIGR01779",#vitamin B12
              "TIGR01782", #polysaccharide
              "TIGR01783",#siderophore
              "TIGR01785",# haem and haemoglobin
              "TIGR01786", #haem and haemoglobin
              "PTHR32552", "PTHR30442", "PTHR30069"
)

#summarize number of TonB transporters
DEqMS_results_TonB<- DEqMS_results%>% 
  left_join(protein_annotations, by = "gene_callers_id", ) %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% unique() %>% 
  mutate(TonB=case_when(NCBIfam_acc %in% SusCD_acc | Pfam_acc%in% SusCD_acc |PANTHER_acc %in% TonBs_acc|
                        grepl("Sus", InterPro_ann, ignore.case =TRUE)|
                          grepl("Sus", blastp_ann, ignore.case =TRUE) ~"SusCD transport system",
                          NCBIfam_acc %in% TonBs_acc |
                          TIGRFAM_acc %in% TonBs_acc |
                          grepl("TonB", InterPro_ann)| 
                          grepl("TonB", ProSitePatterns_ann)|
                          grepl("TONB", PANTHER_ann)|
                          grepl("TONB", blastp_ann)~ "TonB-dependent transport systems")) 


DEqMS_results_TonB %>% filter(NCBIfam_acc=="TIGR01783" | TIGRFAM_acc=="TIGR01783" | grepl("siderophore", blastp_ann, ignore.case=TRUE)) %>% 
  group_by(Enr.group, Enr.frac) %>% 
  summarize(N_p = n())


DEqMS_results_TonB%>% 
  group_by(Enr.group, Enr.frac, TonB) %>% 
  summarize(N_p = n()) %>% 
  left_join(DEqMS_results_totals, by = c("Enr.group", "Enr.frac")) %>% 
  mutate(Prop=N_p/Total_p) %>% 
  #filter(grepl("Gammaproteobacteria", Class)) %>% 
  #filter(grepl("Alteromonadales|Cellvibrionales", Order)) %>% 
  select(Enr.group, Enr.frac, Prop, TonB) %>% 
  tidyr::spread(Enr.group, Prop)%>%  
  View()

DEqMS_results_TonB%>% filter(Enr.frac %in% c("Cells","EVs")) %>% 
  group_by(Enr.group, Enr.frac, Class, Order, TonB) %>% 
  summarize(count = n()) %>% View()

DEqMS_results_TonB%>% 
  group_by(Enr.group, Enr.frac, Class, TonB) %>% 
  summarize(N_p = n()) %>% 
  left_join(DEqMS_results_totals, by = c("Enr.group", "Enr.frac")) %>% 
  mutate(Prop=N_p/Total_p) %>% 
  filter(grepl("Gammaproteobacteria", Class)) %>% 
  #filter(grepl("Alteromonadales|Cellvibrionales", Order)) %>% 
  select(Enr.group, Enr.frac, Class, Prop, TonB) %>% 
  tidyr::spread(Enr.group, Prop)%>%  
  View()

############################
#SusCD proteins
############################
DEqMS_results_Sus<- DEqMS_results%>% 
  left_join(protein_annotations, by = "gene_callers_id", relationship = "many-to-many") %>%  unique() %>% 
  mutate(Sus=case_when(grepl("Susd", Pfam_ann, ignore.case =TRUE)|
                          grepl("SusD", InterPro_ann, ignore.case =TRUE)|
                          grepl("SusD", blastp_ann, ignore.case =TRUE) ~"SusD",
                        NCBIfam_acc %in% SusCD_acc |grepl("SusC", InterPro_ann, ignore.case =TRUE)|
                          grepl("SUSC", PANTHER_ann, ignore.case =TRUE)|
                          grepl("SusC", blastp_ann, ignore.case =TRUE) ~"SusC"))


#how many unique proteins
DEqMS_results_Sus %>% filter(!is.na(Sus)) %>% 
  distinct(gene_callers_id, Enr.frac, Sus) %>% group_by(Enr.frac, Sus) %>% 
  summarise(length(gene_callers_id))


DEqMS_results_Sus %>% filter(!is.na(Sus)) %>% 
  group_by(Enr.group, Enr.frac, Sus) %>%
  summarise(length(gene_callers_id))


DEqMS_results_Sus %>% filter(!is.na(Sus)) %>% 
  group_by(Enr.group, Enr.frac, Sus) %>% 
  summarise(n(), mean(logFC), se(logFC)) %>% View()



#total Flavobacterial enriched proteins
DEqMS_results_Flavo<- DEqMS_results %>% left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  filter(grepl("Flavobacteriales", Order)) %>% 
  group_by(Enr.group, Enr.frac) %>% summarize(Total_p = n()) 


DEqMS_results_Sus%>% filter(!is.na(Sus)) %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>%
  group_by(Enr.group, Enr.frac, Order, Sus) %>% 
  filter(grepl("Flavobacteriales", Order)) %>% 
  summarize(count = n()) %>% 
  left_join(DEqMS_results_Flavo, by = c("Enr.group", "Enr.frac")) %>% 
  mutate(Prop=count/Total_p) %>% 
  group_by(Enr.group, Enr.frac) %>% 
  summarize(Total_prop = sum(Prop)) %>% 
  tidyr::spread(Enr.group, Total_prop)%>%  
  View()

DEqMS_results_Flavo_all<- DEqMS_results %>% left_join(protein_taxonomy, by = "gene_callers_id") %>% unique() %>% 
  filter(grepl("Flavobacteriales", Order)) %>% left_join(protein_annotations, by = "gene_callers_id") 


############################
#Plot enrichment of relevant proteins 
############################
DEqMS_results_plot<- DEqMS_results%>% filter(Enr.frac %in% c("Cells","EVs")) %>% left_join(protein_annotations, by = "gene_callers_id") %>% 
  mutate(Prot_group=case_when(NCBIfam_acc %in% SusCD_acc | Pfam_acc%in% SusCD_acc |
                                grepl("Sus", InterPro_ann, ignore.case =TRUE)|
                                grepl("Sus", blastp_ann, ignore.case =TRUE) ~"SusCD transport system",
                                NCBIfam_acc %in% TonBs_acc |
                                TIGRFAM_acc %in% TonBs_acc |
                                grepl("TonB", InterPro_ann)| 
                                grepl("TonB", ProSitePatterns_ann)|
                                grepl("TONB", PANTHER_ann)|
                                grepl("TONB", blastp_ann)~ "TonB-dependent transport systems",
                              grepl("porin", InterPro_ann, ignore.case = TRUE)|
                                grepl("iron uptake porin", NCBIfam_ann, ignore.case = TRUE)|
                                grepl("porin", Pfam_ann, ignore.case = TRUE)|
                                grepl("iron uptake porin", blastp_ann, ignore.case = TRUE) ~ "Porins")) %>%  
  filter(!is.na(Prot_group)) %>% 
  group_by(Enr.group, Enr.frac, Prot_group) %>% 
  summarize(log2fold_mean = mean(logFC), log2fold_median = median(logFC), log2fold_min = min(logFC), log2fold_max = max(logFC), log2fold_se = se(logFC), count=n()) 


DEqMS_results_plot %>% left_join(DEqMS_results_totals, by =c("Enr.group", "Enr.frac")) %>% 
  mutate(prot_prop= count/Total_p,
         Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN","UP"))) %>% 
  ggplot(aes(y=log2fold_mean , x=Prot_group, colour = Enr.group, label = count))+ 
  #geom_point(size = 2, position = position_dodge(width = 1))+
  geom_text(position = position_dodge(width = 1))+
  geom_point(shape = 21, position = position_dodge(width = 1), aes(size = prot_prop,  colour = Enr.group))+
  scale_size_continuous(range = c(1, 20))+
  geom_errorbar(aes(ymin = log2fold_mean-log2fold_se, ymax = log2fold_mean +log2fold_se), 
                width = 0.2, position = position_dodge(width = 1)) + 
  ylab("log2 foldchange")+
  scale_color_manual(values = c("#009E73", "#F0E442", "#0072B2", 
                                "#D55E00"))+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme_EF+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=90))


#save the plot
ggsave("./Figures/Selected_prot_group_enr.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       #height = 90, 
       scale = 2,
       dpi = 300)

############################
#export results for excel 
############################
DEqMS_results %>% filter(Enr.frac!="Not.enr") %>% 
  left_join(protein_annotations %>% select(gene_callers_id,InterPro_acc,InterPro_ann,NCBIfam_acc,NCBIfam_ann,Pfam_acc,Pfam_ann), by = "gene_callers_id", relationship = "many-to-many") %>% unique() %>% 
  left_join(protein_taxonomy, by = "gene_callers_id") %>% 
  select(Enr.group, Enr.frac, logFC, sca.adj.pval, starts_with("InterPro_"), starts_with("NCBIfam_"),
         starts_with("Pfam_"), Domain, Phylum, Class,Order,Family,Genus) %>% 
  openxlsx::write.xlsx(., colNames = TRUE, 
                       file= 'Tables/Table_S3-DEqMS_results_Fractions.xlsx')

#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
