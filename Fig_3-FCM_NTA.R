###########################################
#Load libraries
###########################################
library(tidyverse)
library(ggpmisc)
library(ggpattern)

source("source/ggplot_parameters.R")

#processggpattern#process NTA output
source("source/process_NTA_results.R")

###########################################
#Plot size distribution
###########################################
tidy_data_filt_total<- tidy_data_filt %>% 
  group_by(Station_ID) %>% 
  summarize(Part.n=n())


EVs_size.p<- tidy_data_filt %>% left_join(tidy_data_filt_total, by="Station_ID") %>% 
  ggplot(aes(x=Station_ID, y=Size.nm, group=Station_ID, fill=Region)) +
  #geom_violin(outliers = TRUE)+ 
  geom_boxplot(outliers = FALSE)+
  #geom_text(aes(label=Part.n), y= 20)+
  #ylim(0,250)+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2",  "#D55E00"))+
  theme_EF+
  theme(axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5))+
  facet_grid(cols=vars(Region),scales="free_x",space="free_x",switch="x")

#save the plot
ggsave("./Figures/NTA_size_dist.pdf",
       plot = EVs_size.p,
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)


#calculate mean of the triplicates selected based on ParticleData
EVs_mean_conc<- SummaryData_df%>% 
  filter(Type=="EVs") %>% 
  left_join(bind_rows(ks_tests), by = "Station_ID") %>%
  filter(Replicate ==replicate) %>%  #keep only the consistent replicates 
  #filter(Type=="EVs", Replicate %in% c("C","D","E")) %>% 
  group_by(Region, Station_ID, Type, Bin.centre..nm.) %>%
  summarise(Mean.conc= mean(Part.conc),
            SD.conc = sd(Part.conc))

###########################################
#Plot total concentration
###########################################
#calculate totals in each replicate and calculate mean
EVs_total_conc<- SummaryData_df%>% 
  filter(Type=="EVs") %>% 
  left_join(bind_rows(ks_tests), by = "Station_ID") %>%
  filter(Replicate ==replicate) %>%  #keep only the consistent replicates 
  group_by(Region, Station_ID, Type, Replicate) %>%
  summarise(Total.conc= sum(Part.conc)) %>% 
  group_by(Region, Station_ID, Type) %>%
  summarise(Mean.conc= 1000*mean(Total.conc),
            SD.conc = 1000*sd(Total.conc))

#plot
EVs_total_conc.p<- EVs_total_conc %>% ggplot(aes(x= Station_ID, y = Mean.conc, group = Type, fill = Region))+
  geom_col(position = "dodge")+ 
  geom_errorbar(aes(ymin = Mean.conc-SD.conc, ymax = Mean.conc+SD.conc), colour = "gray50", width=0.5)+
  labs(y="Concentration (particle L-1)", x = "Station")+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  theme_EF+
  theme(axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5))+
  facet_grid(cols=vars(Region),scales="free_x",space="free_x",switch="x")


###########################################
#Cell counts
###########################################
counts_all<- read.table("data/FCM_cell_counts.txt") %>% 
                mutate(Cell_conc=Cell_conc*1000,
                       Region = factor(Region, levels =c("WEST","GYRE", "TRAN","UP")),
                       Type="Cells")

Cells_conc.p<- counts_all %>% ggplot(aes(x= Station_ID, y = Cell_conc, fill = Region))+
  geom_col(position = "dodge")+ 
  labs(y="Concentration (cells L-1)", x = "Station")+
  #ggbreak::scale_y_break(c(4e8, 6e8), scales= c(1,3))+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  theme_EF+
  theme(axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5))+
  facet_grid(cols=vars(Region),scales="free",space="free_x",switch="x")



FCM_NTA_plot<- counts_all %>% select("Region","Station_ID", "Type", "Cell_conc") %>% 
  dplyr::rename(Concentration=Cell_conc) %>% 
  rbind(EVs_total_conc %>% select("Region","Station_ID", "Type", "Mean.conc") %>% 
          dplyr::rename(Concentration=Mean.conc)) %>% 
  mutate(Region = factor(Region, levels =c("WEST","GYRE", "TRAN","UP"))) %>% 
  ggplot(aes(x= Station_ID, y = Concentration, fill = Region, group =Type))+
  #geom_col(position = "dodge")+ 
  geom_col_pattern(aes(pattern=Type),
                   #pattern = 'stripe',
                   position = position_dodge(width = .8), width=.7, #pattern_density = 0.5,
                   colour ="black")+
  ggbreak::scale_y_break(c(6e8, 1e9), scales= c(1,3))+
  labs(y="Concentration (# L-1)", x = "Station")+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  scale_pattern_manual(values=c('stripe', 'wave'))+
  facet_grid(cols=vars(Region),scales="free",space="free_x",switch="x")+
  theme_EF+
  theme(axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5))

#save the plot
ggsave("./Figures/FCM_NTA_plot.pdf",
       plot = FCM_NTA_plot,
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)


###########################################
#Estimate effect of density gradient separation on concentrations
###########################################
#before density gradient
SummaryData_df_no_dens_total<- SummaryData_df_no_dens%>% 
  group_by(Region,Station_ID,Replicate) %>% 
  summarize(Part.conc=sum(Part.conc)) %>% 
  group_by(Region,Station_ID) %>% 
  summarize(Total.conc= mean(Part.conc), Total.SD= sd(Part.conc))

#after density gradient
SummaryData_df_total<- SummaryData_df %>% 
  group_by(Region,Station_ID,Replicate) %>% 
  summarize(Part.conc=sum(Part.conc)) %>% #sum all fractions to compare the totals
  group_by(Region,Station_ID) %>% 
  summarize(Mean.conc= mean(Part.conc), SD.conc= sd(Part.conc))


#plot abundances
merge(SummaryData_df_total,SummaryData_df_no_dens_total,
      by=c("Region","Station_ID")) %>% 
  select(Region, Station_ID, Mean.conc, Total.conc) %>% 
  reshape2::melt() %>% 
  mutate(Region = factor(Region, levels =c("WEST","GYRE", "TRAN","UP"))) %>% 
  ggplot(aes(x= Station_ID, y = 1000*value, fill = Region, group =variable))+
  #geom_col(position = "dodge")+ 
  geom_col_pattern(aes(pattern=variable),
                   #pattern = 'stripe',
                   position = position_dodge(width = .8), width=.7, #pattern_density = 0.5,
                   colour ="black")+
  #ggbreak::scale_y_break(c(8e8, 1e9), scales= c(1,3))+
  labs(y="Concentration (# L-1)", x = "Station")+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  scale_pattern_manual(values=c('stripe', 'wave'))+
  facet_grid(cols=vars(Region),scales="free",space="free_x",switch="x")+
  theme_EF+
  theme(axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5))

#save the plot
ggsave("./Figures/density_gradient.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

#calculate ratios
merge(SummaryData_df_total,
      SummaryData_df_no_dens_total,
      by=c("Region","Station_ID")) %>% 
  mutate(Ratio=Mean.conc/Total.conc) %>% 
  View()

#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()