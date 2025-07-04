require(dplyr)
require(ggplot2)

#load ggplot theme and colours
source("source/ggplot_parameters.R")

#import metadata 
sample_meta<- read.table("data/samples_meta.txt", h=T)

############################
#Figure 1 - Phytoplankton composition
############################
#import HPLC chemtax data
SO289_chemtax<- read.table("data/SO289_HPLC_chemtax.txt", h = T, sep = "\t")

#summarize taxonomy
SO289_chemtax_sum<- SO289_chemtax %>% 
  mutate(Station_ID = gsub("Stn", "SO289_", Station)) %>% 
  left_join(sample_meta, by ="Station_ID", multiple="first") %>% 
  filter(Depth == 5, Station_ID !="SO289_43wdh") %>% 
  select(Region, Station_ID, Diatoms,Dinoflagellates,Cryptophytes,Chrysophytes,
         Pelagophytes,Haptophytes_3,Haptophytes_4,Prasinophytes,Chlorophytes,
         Synechococcus,Prochloroccocus) %>%
  reshape2::melt(id.vars =c("Region", "Station_ID")) %>% 
  mutate(Station_ID=as.numeric(gsub("SO289_", "", Station_ID))) %>% 
  dplyr::rename(Proportion= value, Taxa = variable) %>% 
  mutate(Taxa=case_when(Taxa %in% c("Haptophytes_3","Haptophytes_4")~"Haptophytes", TRUE~Taxa),
         Region = case_when(Station_ID< 6 ~ "UP",
                            Station_ID<17 & Station_ID>5~ "TRAN",
                            Station_ID<33 & Station_ID>16~ "GYRE", 
                            Station_ID>32~ "WEST")) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels =c(rev(1:45))))
  
#plot
SO289_chemtax_sum %>% 
  ggplot(aes(x=Station_ID, y = Proportion, fill = Taxa))+
  geom_col()+
  ylab("Chemtax_composition")+
  scale_fill_manual(values = tol21rainbow)+
  theme_EF+
  theme(axis.text.x= element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom", 
    plot.title = element_text(hjust = 0.5))+
  facet_grid(cols=vars(Region),scales="free_x",space="free_x",switch="x")

#save the plot
ggsave("./Figures/Fig_1-phytoplankton.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

############################
#Figure S1 - Chlorophyll a and nutrients
############################
#import HPLC data
SO289_HPLC<- read.table("data/SO289_HPLC_pigments.txt", h = T, sep = "\t")

#import nutrients data
SO289_NPFe<- read.table("data/SO289NPFe_HL.csv", h = T, sep = ",", dec =".", row.names = 1)


#calculate total Chl. a
Chla_conc <- SO289_HPLC %>% 
  mutate(Station_ID = gsub("Stn", "SO289_", Station)) %>% 
  left_join(sample_meta, by ="Station_ID", multiple="first") %>% 
  filter(Depth == 5, Station_ID !="SO289_43wdh") %>% 
  mutate(Chla_corr=(Chl_a+Div_a)*0.001) %>% #transform to mg/m-3
  dplyr::rename(Longitude=Lon) %>% 
  select(Longitude, Chla_corr) %>% 
  mutate(Longitude=case_when(Longitude>0~ -180-(180-Longitude), TRUE~Longitude),
       Region = case_when(-72 < Longitude ~ "UP",
                          -72 > Longitude & -100 < Longitude~ "TRAN",
                          -100 > Longitude & -150 < Longitude~ "GYRE", 
                          -150 > Longitude~ "WEST")) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")))

#plot
Chla_conc %>% 
  ggplot(aes(x=Longitude, y=Chla_corr))+
  geom_point(aes(fill = Region), shape=21, size = 4)+
  geom_line(alpha=0.3)+
  labs(y="Chlorophyll a [mg/m-3]")+
  ggbreak::scale_y_break(c(0.5, 2.5), scales= c(1,2))+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  #ylim(0,0.5)+
  xlim(-188, -71)+
  theme_EF+
  theme(legend.position = "none")

ggsave("./Figures/Fig_S1-chla_conc.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

#plot phosphate
plot_P<- SO289_NPFe %>% select(Longitude, Phosphate) %>% 
  reshape2::melt(id.vars=c("Longitude")) %>% 
  mutate(Longitude=case_when(Longitude>0~ -180-(180-Longitude), TRUE~Longitude),
         Region = case_when(-72 < Longitude ~ "UP",
                            -72 > Longitude & -100 < Longitude~ "TRAN",
                            -100 > Longitude & -150 < Longitude~ "GYRE", 
                            -150 > Longitude~ "WEST")) %>%
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP"))) %>% 
  filter(!is.na(value))%>%
  ggplot(aes(x=Longitude, y=value))+
  geom_point(aes(fill = Region), shape=21, size = 2)+
  geom_line(alpha=0.3)+
  labs(y="Phosphate [uM]")+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  ylim(0,0.6)+
  xlim(-188, -71)+
  theme_EF+
  theme(legend.position = "none")

ggsave("./Figures/Fig_S1-P.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

#plot iron
plot_Fe<- SO289_NPFe %>% select(Longitude, dFe) %>% 
  reshape2::melt(id.vars=c("Longitude")) %>% 
  mutate(Longitude=case_when(Longitude>0~ -180-(180-Longitude), TRUE~Longitude),
         Region = case_when(-72 < Longitude ~ "UP",
                            -72 > Longitude & -100 < Longitude~ "TRAN",
                            -100 > Longitude & -150 < Longitude~ "GYRE", 
                            -150 > Longitude~ "WEST")) %>%
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP"))) %>% 
  filter(!is.na(value))%>%
  ggplot(aes(x=Longitude, y=value))+
  geom_point(aes(fill = Region), shape=21, size = 2)+
  geom_line(alpha=0.3)+
  labs(y="Dissolved Fe [nM]")+
  scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2", "#D55E00"))+
  ylim(0,1.5)+
  xlim(-188, -71)+
  theme_EF+
  theme(legend.position = "none")

ggsave("./Figures/Fig_S1-Fe.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)

#plot nitrogen
plot_N<- SO289_NPFe %>% select(Longitude, TON) %>% 
  reshape2::melt(id.vars=c("Longitude")) %>% 
  mutate(Longitude=case_when(Longitude>0~ -180-(180-Longitude), TRUE~Longitude),
         Region = case_when(-72 < Longitude ~ "UP",
                            -72 > Longitude & -100 < Longitude~ "TRAN",
                            -100 > Longitude & -150 < Longitude~ "GYRE", 
                            -150 > Longitude~ "WEST")) %>%
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP"))) %>% 
  filter(!is.na(value))%>%
  ggplot(aes(x=Longitude, y=value))+
  geom_point(aes(fill = Region), shape=21, size = 2)+
  geom_line(alpha=0.3)+
  labs(y = "Nitrate [uM]")+
  #facet_grid(variable~., scales = "free_y")+
  ylim(0,0.5)+
  xlim(-188, -71)+
  theme_EF+
  theme(legend.position = "none")

ggsave("./Figures/Fig_S1-N.pdf",
       plot = last_plot(),
       units = "mm",
       width = 180,
       height = 90, 
       scale = 2,
       dpi = 300)