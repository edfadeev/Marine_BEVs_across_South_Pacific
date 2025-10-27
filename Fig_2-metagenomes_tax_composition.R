require(pavian)
require(dplyr)

#load ggplot theme and colours
source("source/ggplot_parameters.R")

#import metadata
samples_meta_df<- read.table("data/samples_meta.txt", sep ="\t", header = TRUE, row.names = 1) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) 
k_reports<- readRDS("./data/kraken_reports.rds")

##################################################
# Generate summary and proportions from kreports #
##################################################
#reads overview
reports_overview<- summarize_reports(k_reports) %>% 
  tibble::rownames_to_column("Station_ID") %>% 
  mutate(Station_ID =gsub("\\.kreport*","",Station_ID))

#merge taxonomy table
merged_report<- merge_reports(k_reports)
class(merged_report)<- "data.frame"

merged_read_counts<- merged_report%>% 
  select(Name,TaxRank, TaxLineage, ends_with("cladeReads")) %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  mutate(Station_ID=gsub("\\.kreport.*","", variable)) %>% 
  dplyr::rename(Reads=value) %>% 
  select(-variable) %>% 
  left_join(reports_overview, by="Station_ID") %>% 
  left_join(samples_meta_df %>% filter(Fraction  =="Cells"), by="Station_ID") 

##########################################
# Plot Bacterial composition             #
##########################################
BAC_read_proportion<- merged_read_counts %>% 
  filter(grepl("d_Bacteria",TaxLineage)) %>% 
  mutate(Proportion= Reads/bacterial_reads) 


#remove below 0.01% reads
top_class <- BAC_read_proportion %>% 
  mutate(TaxRank= case_when(Name=="Cyanobacteria/Melainabacteria group" ~ "C", TRUE ~ TaxRank)) %>% 
  mutate(Name =case_when(Name=="Cyanobacteria/Melainabacteria group" ~ "Cyanophyceae", TRUE ~ Name)) %>% 
  filter(TaxRank=="C" & Proportion>0.005) %>% pull(Name) %>% unique() %>% sort()

BAC_class_proportion <- BAC_read_proportion %>%  
  mutate(TaxRank= case_when(Name=="Cyanobacteria/Melainabacteria group" ~ "C", TRUE ~ TaxRank)) %>% 
  mutate(Name =case_when(Name=="Cyanobacteria/Melainabacteria group" ~ "Cyanophyceae", TRUE ~ Name)) %>% 
  filter(TaxRank=="C") %>% 
  mutate(Name = case_when(Proportion< 0.005 ~ "Other classes < 0.5%", TRUE ~ Name))

BAC_class_proportion <- BAC_class_proportion %>% 
  mutate(Name = factor(Name,levels=c(top_class,"Other classes < 0.5%", "Unclassified")))

#complement all the unclassifed reads 
BAC_class_uncl <- BAC_class_proportion %>%  
  select(Name, Region,  Station_ID, Proportion) %>% 
  group_by(Region,Station_ID) %>% 
  summarize(Total = sum(Proportion)) %>% 
  mutate(Proportion = 1- Total, Name = "Unclassified") %>% 
  select(Name, Region, Station_ID, Proportion)

BAC_read_total<- merged_read_counts %>% select(Station_ID, bacterial_reads) %>% unique() %>% 
  rstatix::get_summary_stats(bacterial_reads, show=c("min","max","sd","mean"))

#plot
BAC_class_proportion %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) %>% 
  select(Name, Region, Station_ID, Proportion) %>% 
  rbind(BAC_class_uncl) %>% 
  ggplot(aes(x= Station_ID, y= Proportion*100, fill = Name)) +
  geom_col()+
  labs(y ="Proportion of bacterial reads (%)")+
  guides(fill = guide_legend(title="Class", ncol = 3))+
  scale_fill_manual(values = tol21rainbow)+
  #geom_text(aes(x=Station_ID, label=Total_reads), y=100)+
  theme_EF+
  theme(axis.line = element_blank(),
        #axis.text.x= element_blank(),
        axis.text.x= element_text(angle=90),
        axis.title.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5))+
  facet_grid(cols=vars(Region),scales="free_x",space="free_x",switch="x")


#save the plot
ggsave("./Figures/metaG_class_comp.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)