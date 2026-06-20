###########################################
#Load libraries
###########################################
library(tidyverse)
library(ggpmisc)
library(ggpattern)

# load custom ggplot theme and colour palettes
source("R/utils.R")

###########################################
#process NTA and FCM data 
###########################################
run_source_if_missing("data/derived/NTA_summary_data.rds", "R/process_NTA.R")
run_source_if_missing("data/derived/FCM_cell_counts.txt", "R/process_FCM.R")

###########################################
#Plot size distribution
###########################################
NTA_dfs<- readRDS("data/derived/NTA_summary_data.rds") 

tidy_data_filt_total<- NTA_dfs$tidy_data_filt |> 
  filter(Type=="BEVs") |>
  group_by(Region) |> 
  summarize(Part.n=n())

#plot size distribution
EVs_size.p<- NTA_dfs$tidy_data_filt |> 
  filter(Type=="BEVs") |>
  left_join(tidy_data_filt_total, by = "Region") |> 
  ggplot(aes(x=Region, y=Size.nm, group=Region, fill=Region)) +
  #geom_violin(draw_quantiles = c(0.25, 0.75))+
  geom_boxplot(outlier.shape = NA, width = 0.5)+
  geom_text(data = tidy_data_filt_total, aes(x = Region, y = Inf, label = paste0("n=", Part.n)), vjust = 1.5, size = 3, inherit.aes = FALSE) +
  #ylim(0,250)+
  #geom_hline(yintercept =70)+
  #geom_hline(yintercept =130)+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2",  "#D55E00"))+
  facet_grid(cols=vars(Region),scales="free_x",space="free_x",switch="y")+
    theme_EF+
  theme(axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",  
        plot.title = element_text(hjust = 0.5))

#save the plot
ggsave("./Figures/NTA_size_dist.pdf",
       plot = EVs_size.p,
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

#calculate mean of the triplicates selected based on ParticleData
EVs_mean_conc<- NTA_dfs$SummaryData_df|> 
  filter(Type=="BEVs") |> 
  group_by(Region, Station_ID, Type, Bin.centre..nm.) |>
  summarise(Mean.conc= mean(Part.conc),
            SD.conc = sd(Part.conc))

###########################################
#Plot total concentration
###########################################
#calculate totals in each replicate and calculate mean
EVs_total_conc<- NTA_dfs$SummaryData_df|> 
  filter(Type=="BEVs") |> 
  group_by(Region, Station_ID, Type, Replicate) |>
  summarise(Total.conc= sum(Part.conc)) |> 
  group_by(Region, Station_ID, Type) |>
  summarise(Mean.conc= 1000*mean(Total.conc), #convert to particles per L
            SD.conc = 1000*sd(Total.conc))

#plot
EVs_total_conc.p<- EVs_total_conc |> ggplot(aes(x= Station_ID, y = Mean.conc, group = Type, fill = Region))+
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
#import cell counts
counts_all<- read.table("data/derived/FCM_cell_counts.txt") |> 
                mutate(Region = factor(Region, levels =c("WEST","GYRE", "TRAN","UP")),
                       Type="Cells")

#plot cell counts
Cells_conc.p<- counts_all |> ggplot(aes(x= Station_ID, y = Cell_conc, fill = Region))+
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


#combined figure cells and BEVs abundance
FCM_NTA_combined <- counts_all |> select("Region","Station_ID", "Type", "Cell_conc") |> 
  dplyr::rename(Concentration=Cell_conc) |> 
  rbind(EVs_total_conc |> select("Region","Station_ID", "Type", "Mean.conc") |> 
          dplyr::rename(Concentration=Mean.conc)) |> 
  mutate(Region = factor(Region, levels =c("WEST","GYRE", "TRAN","UP"))) 


FCM_NTA_plot<- FCM_NTA_combined |>
  ggplot(aes(x= Type, y = Concentration, fill = Region, group =interaction(Region,Type)))+
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  labs(y="Concentration (# L-1)", x = "Fraction")+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  scale_colour_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  facet_grid(~Region)+
  scale_y_log10()+
  theme_EF+
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), strip.text.x = element_blank())


cowplot::plot_grid(EVs_size.p, FCM_NTA_plot, ncol = 1, align = "hv")

#save the plot
ggsave("./Figures/FCM_NTA_plot.pdf",
       plot = last_plot(),
       units = "mm",
       width = 90,
       height = 90, 
       scale = 2,
       dpi = 300)

#cells BEVs proportions
cell_BEVs_prop<- counts_all |> select("Region","Station_ID", "Type", "Cell_conc") |> 
  dplyr::rename(Concentration=Cell_conc) |> 
  rbind(EVs_total_conc |> select("Region","Station_ID", "Type", "Mean.conc") |> 
          dplyr::rename(Concentration=Mean.conc)) |> 
  mutate(Region = factor(Region, levels =c("WEST","GYRE", "TRAN","UP"))) |> 
  spread(Type, Concentration) |> 
  mutate(Proportion=BEVs/Cells)


###########################################
#Estimate effect of density gradient separation on concentrations
###########################################
SummaryData_df_total<- NTA_dfs$SummaryData_df |> 
  mutate(den_sep=case_when(Type %in% c("BEVs", "Soluble", "High-density") ~ "DG",
                           TRUE ~ "Total")) |>
  group_by(Region, Station_ID, Type, den_sep, Replicate) |> 
  summarize(Part.conc=sum(Part.conc)) |>
  group_by(Region, Station_ID, Type, den_sep) |> 
  summarize(Mean.conc= mean(Part.conc)) |> 
  group_by(Region,Station_ID, den_sep) |> 
  summarize(Total.conc=sum(Mean.conc)) 

#plot abundances
SummaryData_df_total |> 
  mutate(Region = factor(Region, levels =c("WEST","GYRE", "TRAN","UP"))) |> 
  ggplot(aes(x= Station_ID, y = 1000*Total.conc, fill = Region, group =den_sep))+
  #geom_col(position = "dodge")+ 
  geom_col_pattern(aes(pattern=den_sep),
                   #pattern = 'stripe',
                   position = position_dodge(width = .8), width=.7, #pattern_density = 0.5,
                   colour ="black")+
  scale_y_log10()+
  labs(y="Concentration (# L-1)", x = "Station")+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  scale_pattern_manual(values=c('stripe', 'wave'))+
  facet_grid(cols=vars(Region),scales="free",space="free_x",switch="x")+
  labs(pattern="Type", fill="Region")+
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
       height = 135, 
       scale = 1,
       dpi = 300)

#calculate ratios
SummaryData_df_total|> 
  spread(den_sep, Total.conc) |>
  mutate(Ratio=DG/Total) |> 
  ungroup() |> 
  summarise(Mean_ratio= mean(Ratio), SE_ratio=se(Ratio))


#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()

