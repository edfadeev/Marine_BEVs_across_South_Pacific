library(flowWorkspace)
library(ncdfFlow)
library(flowAI)
library(ggcyto)
library(ggpmisc)
library(tidyr)

################################################################
###Import FCS files                                          ###
################################################################
#define directory with the FCS files
path.CC <- "data/FCM/"

#create file list and metadata table
fcsFiles <- list.files(path = paste0(path.CC,"FSC_files"), 
                       full = TRUE, recursive = TRUE, pattern = "*.fcs")

#create a flowSet
fs <- read.ncdfFlowSet(fcsFiles)

#quality control using flowAI
fs_comp_clean <- flow_auto_qc(fs)

################################################################
###Define gating configurations                              ###
################################################################
# Initialize a single frame
data.1frame <- fs_comp_clean[[1]]
# fill the single frame with the exprs data from each frame in the flow set
exprs(data.1frame) <- fsApply(fs_comp_clean, function(x) {
  x <- exprs(x)
  return(x)
})

autoplot(data.1frame, x = "FSC-A", "SSC-A")+
  scale_x_logicle() + scale_y_logicle()

autoplot(data.1frame, x = "SSC-A", "Sybr Green 1-A")+
  scale_x_logicle() + scale_y_logicle()

#transform the combined data
chnls <- c("FSC-A", "SSC-A", "FSC-H", "SSC-H","Sybr Green 1-A","Sybr Green 1-H", "PI-A", "PI-H")
trans <- estimateLogicle(data.1frame, channels = chnls)
inv.trans <- inverseLogicleTransform(trans)
data.1frame <- transform(data.1frame, trans)
autoplot(data.1frame, x = "FSC-A", "SSC-A")

#select population
population <- openCyto::gate_flowclust_2d(data.1frame, xChannel = "Sybr Green 1-A", 
                                          yChannel = "SSC-A", K = 3)

autoplot(data.1frame, "SSC-A", "Sybr Green 1-A")+ geom_gate(population)+theme_bw()+theme(legend.position = "none")

#inverse transform the gate to fit all samples
population_gate <- transform(population, inv.trans)

#save the plot
ggsave("./Figures/FCM_gating_plot.pdf",
       plot = last_plot(),
       units = "mm",
       width = 90,
       height = 90, 
       scale = 2,
       dpi = 300)

################################################################
###Process all samples separately                            ###
################################################################
#generate gating set
gs_clean <- GatingSet(fs_comp_clean)

#select stained cells
gs_pop_add(gs_clean, population_gate, parent = "root", name="Population")
recompute(gs_clean)

autoplot(gs_clean, y="Sybr Green 1-A", x="SSC-A", "Population")+
  scale_x_logicle() + scale_y_logicle()

#get the cell counts
nodes <- c("root","Population")
counts_raw <- as(gs_pop_get_stats(gs_clean, nodes, "count"),"data.frame") %>% 
  spread(pop, count) %>% 
  mutate(Prop= Population/root)


################################################################
###Calculate cell abundance                                  ###
################################################################
#import measured flow rates
flow_rate <- read.csv(paste0(path.CC,"flow_rates.csv"), h=T)              

#plot regressions
flow_rate %>% 
  ggplot(aes(Flowrate, Volume)) + 
  geom_point(size =3)+
  geom_smooth(method = "lm")+
  scale_y_continuous(labels=scales::scientific_format())+
  stat_poly_eq(formula = y~x,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,
               colour = "black")+
  theme_bw()


#produce linear model for the flow rate
flow_rate.lm <- lm(Volume ~ Flowrate, data = flow_rate)


#import metadata
meta_data <- read.table(paste0(path.CC,"sample_data.txt")) %>% 
  mutate(Volume = predict(flow_rate.lm, .), #add volumes according to the flow rate
         Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1")))
counts_all <- counts_raw %>% 
  separate(sample, into = c("Cruise","x","Sample_ID"), sep ="_") %>% 
  mutate(Sample_ID = gsub("\\.fcs","", Sample_ID)) %>% 
  select(-c(x)) %>% 
  left_join(meta_data, by = "Sample_ID") %>% 
  mutate(Cell_conc = 1000*Dilution_factor*Population/Volume,
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) %>% 
  filter(!is.na(Station_ID)) 

#save counts
write.table(counts_all,"data/cell_counts.txt")

#plot
counts_all %>% mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
                      Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                                 "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                                 "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                                 "SO289_3", "SO289_1"))) %>% 
  ggplot(mapping=aes(x=Station_ID,y=Cell_conc,fill=Region)) +
  geom_col(color="black") +
  scale_y_continuous(labels=scales::scientific_format())+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", 
                               "#D55E00"))+
  facet_grid(cols=vars(Region),scales="free_x",space="free_x",switch="x") +
  theme_EF +
  theme(axis.title=element_text(size=12,face="bold"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=90),
        axis.text=element_text(size=10,color="black"),
        axis.line = element_line(colour = 'gray50', size = 1),
        strip.placement='outside',
        strip.background.x=element_blank(),
        strip.text=element_text(size=12,color="black",face="bold"),
        panel.spacing.x=unit(10,"pt"),
        panel.border = element_blank())

counts_all %>% 
  #mutate(Region = case_when(Region %in% c("UP","TRAN")~"EAST", TRUE~ Region)) %>%
  group_by(Region) %>% 
  as.data.frame() %>% 
  rstatix::get_summary_stats(Cell_conc, show = c("min","max","mean","se"))

counts_all %>% 
 # mutate(Region = case_when(Region %in% c("UP","TRAN")~"EAST", TRUE~ Region)) %>%
  group_by(Region) %>% 
  as.data.frame() %>% 
  rstatix::t_test(data =., Cell_conc~Region)


###################
#print session info and clean the workspace
###################
sessionInfo()
rm(list = ls())
gc()
dev.off(dev.list()["RStudioGD"])
