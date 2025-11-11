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
# Explore enriched proteins of TonB receptors
############################
#TonB proteins
TonBs_acc<- c("TIGR01352", #TonB family C-terminal domain
              "TIGR01776",
              "TIGR01778", #copper
              "TIGR01779",#vitamin B12
              "TIGR01782", #polysaccharide
              "TIGR01783",#siderophore
              "TIGR01785",# haem and haemoglobin
              "TIGR01786", #haem and haemoglobin
              "TIGR04057", "TIGR04056") #SusC proteins

#extract only enriched TonB-related proteins
DEqMS.results_TonB<- DEqMS_results %>% 
  left_join(protein_annotations %>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_taxonomy %>% unique(), by ="gene_callers_id") %>%
  left_join(protein_metadata, by ="gene_callers_id") %>% 
  mutate(Function=case_when(NCBIfam_acc %in% TonBs_acc ~ NCBIfam_ann, 
                            grepl("TonB", InterPro_ann, ignore.case = TRUE)~InterPro_ann,
                            grepl("TONB", blastp_ann, ignore.case = TRUE)~ gsub("\\[.*|MULTISPECIES:","",blastp_ann),
                            TRUE~"Other")) %>% 
  filter(Function!="Other")

#summarize how many are enriched in each fraction/region
DEqMS.results_TonB %>%
  group_by(Enr.group, Enr.frac) %>% 
  summarize(TonB_p = n()) %>% 
  left_join(DEqMS_results %>% group_by(Enr.group, Enr.frac) %>% 
              summarize(Total_p = n()), by =c("Enr.group", "Enr.frac")) %>% 
  mutate(Prop_TonB=TonB_p/Total_p)
  

############################
#Outer-membrane protein enrichment analysis between fractions in each region
############################
OM_prot_gcids<- protein_annotations %>% 
                              filter(deeploc%in% c("Outer Membrane")) %>% pull(gene_callers_id)

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
  
  prot.dat.log2_norm.filter <- prot.dat.log2_norm[OM_prot_gcids,]
  prot.dat.log2_norm.filter <- prot.dat.log2_norm.filter[,c(EV_sample_IDs,Cell_sample_IDs)]
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
                                ~ "BEVs",
                                -1> logFC & sca.adj.pval<0.1  
                                ~ "Cells", TRUE ~ "Not.enr"),
           log.sca.pval = -log10(sca.P.Value),
           Enr.group =x) 
  
  return(DEqMS.results)
})


############################
#Explore results
############################
DEqMS_OM_results<- bind_rows(enrichment_tests_list) %>% 
  filter(Enr.frac!="Not.enr")

write.table(DEqMS_OM_results, "data/DEqMS_OM_results_regions.txt", col.names =T, row.names = F, quote = F)

#total of different proteins
DEqMS_OM_results %>% select(Enr.frac, gene_callers_id) %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())


#extract only enriched TonB-related proteins
DEqMS_OM_TonB<- DEqMS_OM_results %>% 
  left_join(protein_annotations %>% unique(), by ="gene_callers_id") %>% 
  left_join(protein_taxonomy %>% unique(), by ="gene_callers_id") %>%
  left_join(protein_metadata, by ="gene_callers_id") %>% 
  mutate(Function=case_when(NCBIfam_acc %in% TonBs_acc ~ NCBIfam_ann, 
                            grepl("TonB", InterPro_ann, ignore.case = TRUE)~InterPro_ann,
                            grepl("TONB", blastp_ann, ignore.case = TRUE)~ gsub("\\[.*|MULTISPECIES:","",blastp_ann),
                            TRUE~"Other")) %>% 
  filter(Function!="Other")

#total number of TonB in each fraction
DEqMS_OM_TonB %>% select(Enr.frac, gene_callers_id) %>% unique() %>% 
  group_by(Enr.frac) %>% 
  summarize(Total_p = n())

#overlap between fractions
enr_tonB_overlaps<-VennDiagram::calculate.overlap(x=list("BEVs"=DEqMS_OM_TonB %>% filter(Enr.frac=="BEVs") %>% pull(gene_callers_id) %>% unique(),
                                      "Cells"=DEqMS_OM_TonB %>% filter(Enr.frac=="Cells") %>% pull(gene_callers_id)%>% unique()))
summary(enr_tonB_overlaps)

#Taxonomy of enriched TonB
DEqMS_OM_TonB_tax<- DEqMS_OM_TonB %>% 
  group_by(Enr.frac, Phylum, Class, Order) %>% unique() %>% 
  summarize(N_p = n()) 

############################
#characterize the ligands of the TonB receptors
############################
DEqMS_OM_TonB_seqs<- DEqMS_OM_TonB%>% 
  select(gene_callers_id, aa_sequence) %>% unique() %>% 
  mutate(gene_callers_id=paste0("gcid_",gene_callers_id))

seqinr::write.fasta(as.list(DEqMS_OM_TonB_seqs$aa_sequence), names=DEqMS_OM_TonB_seqs$gene_callers_id, 
                    as.string=FALSE, file.out="data/TonB_prot.fa")

#import psiblast results
blastp.out<- read.table("data/TonB_psiblastp.out", col.names =c("gene_callers_id","sseqid","pident","length","evalue","bitscore")) %>% 
  filter(pident>30  & length>100 & evalue<0.001 & bitscore >50 ) %>% 
  mutate(gene_callers_id=gsub("gcid_","", gene_callers_id))
TCDB_annotation<- read.table("data/TCDB_1_B.ann", h=F, sep="\t")
names(TCDB_annotation)<- c("sseqid", "TCID", "Annotation")


TCDB_TonB_desc<- read.csv("data/TCDB_1B_description.txt", sep ="\t") %>% 
  filter(TCID!="") %>% 
  mutate(Name=iconv(Name, from = "latin1", to = "UTF-8", sub = ""))


#best hit by identitiy
TonB_ann<- blastp.out %>% 
  left_join(TCDB_annotation, by ="sseqid") %>%
  group_by(gene_callers_id) %>% 
  slice_max(pident, n = 1) %>% 
  left_join(TCDB_TonB_desc, by="TCID") %>% 
  mutate(TonB_category=case_when(grepl("SusC",Name, ignore.case=TRUE)~"SusC-like",
                                 grepl("siderophore",Name, ignore.case=TRUE)~"Siderophores",
                                 grepl("catechol|chelin",Name, ignore.case=TRUE)~"Siderophores",
                                 grepl("bactin",Name, ignore.case=TRUE)~"Siderophores",
                                 grepl("cobalamin",Name, ignore.case=TRUE)~"Cobalamine",
                                 grepl("FecA",Name, ignore.case=TRUE)~"Ferric citrate",
                                 grepl("CfrA|RagA",Name, ignore.case=TRUE)~"Iron",
                                 grepl("Ferrioxamine",Name, ignore.case=TRUE)~"Ferrioxamine",
                                 grepl("heme",Name, ignore.case=TRUE)~"Heme",
                                 grepl("Transferrin",Name, ignore.case=TRUE)~"Transferrin",
                                 grepl("oligosaccharide",Name, ignore.case=TRUE)~"Oligosaccharides",
                                 grepl("salicin",Name, ignore.case=TRUE)~"Glucosides",
                                 grepl("collagenase",Name, ignore.case=TRUE)~"Collagenases",
                                TRUE~"Unknown ligands")) %>% 
  mutate(TonB_category=case_when(grepl("FecA3 ",Proteins, ignore.case=TRUE)~"Nickel",TRUE~TonB_category)) 


#add SusD enriched proteins
DEqMS_OM_results_SusD<- DEqMS_OM_results %>% 
  left_join(protein_annotations %>% unique(), by ="gene_callers_id") %>%
  left_join(protein_taxonomy %>% unique(), by ="gene_callers_id") %>%
  left_join(protein_metadata, by ="gene_callers_id") %>% 
  mutate(TonB_category=case_when(grepl("SusD", InterPro_ann, ignore.case = TRUE)~"SusD-like",
                            TRUE~"Other")) %>% 
  filter(TonB_category!="Other") %>% 
  select(Enr.group, Enr.frac, Class, TonB_category, logFC)

#merge and summarize
DEqMS_OM_TonB_total_by_frac<- DEqMS_OM_TonB %>% 
  left_join(TonB_ann,  by ="gene_callers_id") %>% 
  mutate(TonB_category=case_when(is.na(TonB_category)~"Unknown ligands",TRUE~TonB_category)) %>%
  select(Enr.group, Enr.frac, TonB_category, logFC) %>% 
  rbind(DEqMS_OM_results_SusD %>% select(Enr.group, Enr.frac, TonB_category, logFC)) %>% 
  group_by(Enr.frac, TonB_category) %>% 
  summarize(p=n()) %>% 
  spread(Enr.frac, p) %>% 
  mutate(Cells=case_when(is.na(Cells)~0, TRUE~Cells),
         BEVs=case_when(is.na(BEVs)~0, TRUE~BEVs)) 

#plot results
DEqMS_OM_TonB.p<- DEqMS_OM_TonB %>% 
  left_join(TonB_ann,  by ="gene_callers_id") %>% 
  mutate(TonB_category=case_when(is.na(TonB_category)~"Unknown ligands",TRUE~TonB_category)) %>% 
  select(Enr.group, Enr.frac, TonB_category, logFC) %>% 
  rbind(DEqMS_OM_results_SusD %>% select(Enr.group, Enr.frac, TonB_category, logFC)) %>% 
  group_by(Enr.group, Enr.frac,  TonB_category) %>% 
  summarize(log2fold_mean = mean(logFC), log2fold_median = median(logFC), log2fold_min = min(logFC), log2fold_max = max(logFC), log2fold_se = se(logFC), count=n())

total_TonB_prot<- DEqMS_OM_TonB.p %>% group_by(Enr.frac,Enr.group) %>% summarize(Total_p=sum(count))

DEqMS_OM_TonB.p %>%
  left_join(total_TonB_prot, by =c("Enr.frac", "Enr.group")) %>% 
  mutate(Enr.group = factor(Enr.group, levels =c("WEST","GYRE", "TRAN")),
         TonB_category= factor(TonB_category, levels =c("Unknown ligands", "Cobalamine","Collagenases", "Nickel", 
                                                        "Oligosaccharides", "Glucosides","SusC-like", "SusD-like",
                                                         "Ferric citrate", "Ferrioxamine", "Heme", "Iron", "Siderophores", "Transferrin")),
         Prop=count/Total_p) %>% 
  ggplot(aes(y=TonB_category , x=log2fold_mean, fill = Enr.group, label = count))+ 
  geom_point(aes(size = Prop), shape =21)+
  geom_text(size = 5, nudge_y = -0.2)+
  geom_errorbar(aes(xmin = log2fold_mean-log2fold_se, xmax = log2fold_mean +log2fold_se), 
                width = 0.2) + 
  xlim(-7,7)+
  facet_grid(~Enr.group)+  
  geom_vline(aes(xintercept=0), linetype="dashed")+
  scale_size_continuous(range = c(10, 30))+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  theme_EF+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=90))


#save the plot
ggsave("./Figures/TonB_proteins_regions.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)


#add taxonomy pie charts 

#plot results
DEqMS_OM_TonB %>% 
  left_join(TonB_ann,  by ="gene_callers_id") %>% 
  mutate(TonB_category=case_when(is.na(TonB_category)~"Unknown ligands",TRUE~TonB_category)) %>% 
  select(Enr.group, Enr.frac, Class, TonB_category, logFC) %>% 
  rbind(DEqMS_OM_results_SusD) %>% 
  mutate(TonB_category=case_when(is.na(TonB_category)~"Unknown ligands",TRUE~TonB_category),
         Class=case_when(is.na(Class)~"Unclassified",TRUE~Class)) %>% 
  dplyr::count(Enr.group, Enr.frac, Class, TonB_category) %>%
  group_by(Enr.group, Enr.frac, TonB_category) %>% # Re-group by the desired summation level
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x="", y=freq, fill=Class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_manual(values=c(cbbPalette, tol21rainbow))+
  facet_wrap(Enr.group~Enr.frac+TonB_category)
  

#save the plot
ggsave("./Figures/TonB_proteins_regions_pies.pdf",
       plot = last_plot(),
       units = "mm",
       #width = 180,
       #height = 90, 
       scale = 2,
       dpi = 300)

  

  