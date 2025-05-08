###########################################
#Load libraries
###########################################
library(tidyverse)
library(ggpmisc)
#source("R_scripts/extra_R_functions.R")

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}


station_IDs <- read.table("data/samples_meta.txt", h=T) %>% 
  mutate(Region = factor(Region, levels = c("WEST","GYRE","TRAN","UP")),
         Station_ID = factor(Station_ID, levels = c("SO289_44", "SO289_43", "SO289_41",  "SO289_39", "SO289_37", "SO289_34",
                                                    "SO289_33", "SO289_32", "SO289_30", "SO289_27", "SO289_23", "SO289_20", 
                                                    "SO289_17", "SO289_16", "SO289_13", "SO289_12", "SO289_9", "SO289_6", 
                                                    "SO289_3", "SO289_1"))) %>% 
  select(Station_name, Station_ID, Region) %>% unique()

#define working directory with the files
path.NTA <- "C:/Users/fadee/ucloud/Projects/Marine_EVs/SO289/01_Results/NTA/DensityGradient/"


###########################################
#Particle size distributions
###########################################
#list all ParticleData files, import them and aggregate to a data frame
ParticleData_files <- list.files(path=path.NTA, pattern = "*_ParticleData.csv",
                                 full.names = TRUE, recursive = T)

particleData_list <- lapply(ParticleData_files, function(x) {
  
  dat <- read.csv(x,header=TRUE, sep=",", dec = ".", fileEncoding = "utf8") %>% 
    filter(Included.in.distribution.=="True") %>% 
    select(Particle.ID,Size.nm) %>% 
    mutate_all(as.numeric)
  
  # Add a column with the sampling time
  label <- unlist(strsplit(basename(x), split ="_"))
  
  dat<- dat %>% 
    mutate(StationID = label[1],
           Media = label[2],
           Dye = label[3],
           Dilution = label[4],
           Method = label[5],
           stamp = label[6]) %>% 
    separate(StationID, into=c("Station_name","Fraction"), sep ="-") %>% 
    separate(stamp, into =c("Flow_rate","Date","Stamp"), sep =" ") 
  return(dat)
})


#aggregate all the samples into a single dataframe
tidy_data <- bind_rows(particleData_list) %>%
  mutate(Tech_run = paste0(Date,Stamp),
         Type = factor(case_when(Fraction %in% c("DEFG","HI")~ "EVs",
                                 Fraction =="JP"~ "High-density",
                                 Fraction =="ABC"~ "Soluble"),levels = c("Soluble", "EVs", "High-density"))) %>% 
  left_join(station_IDs, by = "Station_name", keep = FALSE) %>% 
  group_by(Region, Station_ID,Type, Fraction) %>%
  mutate(Replicate = factor(Tech_run, 
                            levels = unique(Tech_run),
                            labels = c("A","B","C","D","E")[1:length(unique(Tech_run))])) %>% 
  filter(Size.nm < 250, Size.nm > 20)


#Kolmogorov-Smirnov Tests to identify the most similar technical replicates in each samples
#maximum of three replicates per sample were used 
ks_tests<- lapply(station_IDs$Station_ID, function(x) {
  #empty dataframe
  dk <- data.frame(Station_ID=character(10),pair=NA, D=numeric(10), p=numeric(10) )
  
  for(j in 1:10){  
    pair<- c("A_B","A_C","A_D","A_E","B_C","B_D","B_E","C_D","C_E","D_E")[j]
    k <- ks.test(x= tidy_data %>% filter(Type=="EVs" & Station_ID== x & Replicate == str_split_1(pair,"_")[1]) %>% pull("Size.nm"), 
                 y = tidy_data %>% filter(Type=="EVs" & Station_ID== x & Replicate == str_split_1(pair,"_")[2])%>% pull("Size.nm"))
    dk$Station_ID<- x
    dk$pair[j] <- c("A_B","A_C","A_D","A_E","B_C","B_D","B_E","C_D","C_E","D_E")[j]
    dk$D[j]       <- k$statistic
    dk$p[j]       <- k$p.value
    
  }
  
  replicates<- dk %>% filter(p>0.05) %>% tidyr::separate(pair, into = c("Rep_1","Rep_2"), "_") %>% 
    select(Station_ID, Rep_1,Rep_2) %>% 
    reshape2::melt("Station_ID") %>% 
    mutate(replicate=factor(value, c("A","B","C","D","E"))) %>% 
    group_by(Station_ID, replicate) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    head(n = 3)
  
  return(replicates)
})


tidy_data_filt<- tidy_data%>% 
  filter(Type=="EVs") %>% 
  left_join(bind_rows(ks_tests), by = "Station_ID") %>%
  filter(Replicate ==replicate) 




###########################################
#Import summary concentration
###########################################
#list all Experiment files
Summary_files <- list.files(path=path.NTA,
                            pattern = "_Summary.csv",
                            full.names = TRUE, recursive = F)


#generate a list of all experiments 
SummaryData_list <- lapply(Summary_files, function(x) {
  
  dat <- read.csv(x,header=TRUE, sep=",", dec = ".",
                  fileEncoding = "utf8", skip = 93, nrows = 1000) %>% 
    select(Bin.centre..nm.,Concentration..particles...ml.) %>% 
    mutate_all(as.numeric)
  
  # Add a column with the sampling time
  label <- unlist(strsplit(basename(x), split ="_"))
  
  dat<- dat %>% 
    mutate(StationID = label[1],
           Media = label[2],
           Dye = label[3],
           Dilution = label[4],
           Method = label[5],
           stamp = label[6]) %>% 
    separate(StationID, into=c("Station_name","Fraction"), sep ="-") %>% 
    separate(stamp, into =c("Flow_rate","Date","Stamp"), sep =" ") 
  return(dat)
})


#aggregate all the samples into a single dataframe
SummaryData_df <- bind_rows(SummaryData_list) %>%
  mutate(Tech_run = paste0(Date,Stamp),
         Type = factor(case_when(Fraction %in% c("DEFG","HI")~ "EVs",
                                 Fraction =="JP"~ "High-density",
                                 Fraction =="ABC"~ "Soluble"),levels = c("Soluble", "EVs", "High-density"))) %>% 
  left_join(station_IDs, by = "Station_name", keep = FALSE) %>% 
  group_by(Region, Station_ID,Type, Fraction) %>%
  mutate(Replicate = factor(Tech_run, 
                            levels = unique(Tech_run),
                            labels = c("A","B","C","D","E")[1:length(unique(Tech_run))])) %>% 
  filter(Bin.centre..nm. < 250, Bin.centre..nm. > 20) %>% 
  group_by(Region, Station_ID, Type, Fraction, Replicate, Bin.centre..nm.) %>%
  summarize(Conc=Concentration..particles...ml./720) %>% #correct for x60 (vivaflow)  and x30 (vivaspin) dilutions and x2.5 concentration for NTA
  group_by(Region, Station_ID, Type, Replicate, Bin.centre..nm.) %>%
  summarize(Part.conc=sum(Conc))#


###########################################
#Process NTA results of samples before density gradient
###########################################
#define working directory with the files
path.NTA <- "C:/Users/fadee/ucloud/Projects/Marine_EVs/SO289/01_Results/NTA/No_densityGradient/"

#list all ParticleData files, import them and aggregate to a data frame
ParticleData_files <- list.files(path=path.NTA, pattern = "*_ParticleData.csv",
                                 full.names = TRUE, recursive = T)

particleData_list <- lapply(ParticleData_files, function(x) {
  
  dat <- read.csv(x,header=TRUE, sep=",", dec = ".", fileEncoding = "utf8") %>% 
    filter(Included.in.distribution.=="True") %>% 
    select(Particle.ID,Size.nm) %>% 
    mutate_all(as.numeric)
  
  # Add a column with the sampling time
  label <- unlist(strsplit(basename(x), split ="_"))
  
  dat<- dat %>% 
    mutate(Station_name = label[1],
           Media = label[2],
           Dye = label[3],
           Dilution = label[4],
           Method = label[5],
           stamp = label[6]) %>% 
    #separate(StationID, into=c("Station_name","Fraction"), sep ="-") %>% 
    separate(stamp, into =c("Flow_rate","Date","Stamp"), sep =" ") 
  return(dat)
})


#aggregate all the samples into a single dataframe
tidy_data <- bind_rows(particleData_list) %>%
  mutate(Tech_run = paste0(Date,Stamp)) %>% 
  group_by(Station_name) %>%
  mutate(Replicate = factor(Tech_run, 
                            levels = unique(Tech_run),
                            labels = c("A","B","C","D","E")[1:length(unique(Tech_run))])) %>% 
  left_join(station_IDs, by = "Station_name") %>% 
  filter(Size.nm < 250, Size.nm > 20)


#Kolmogorov-Smirnov Tests to identify the most similar technical replicates in each samples
#maximum of three replicates per sample were used 
ks_tests<- lapply(station_IDs$Station_ID, function(x) {
  #empty dataframe
  dk <- data.frame(Station_ID=character(10),pair=NA, D=numeric(10), p=numeric(10) )
  
  for(j in 1:10){  
    pair<- c("A_B","A_C","A_D","A_E","B_C","B_D","B_E","C_D","C_E","D_E")[j]
    k <- ks.test(x= tidy_data %>% filter(Station_ID== x & Replicate == str_split_1(pair,"_")[1]) %>% pull("Size.nm"), 
                 y = tidy_data %>% filter(Station_ID== x & Replicate == str_split_1(pair,"_")[2])%>% pull("Size.nm"))
    dk$Station_ID<- x
    dk$pair[j] <- c("A_B","A_C","A_D","A_E","B_C","B_D","B_E","C_D","C_E","D_E")[j]
    dk$D[j]       <- k$statistic
    dk$p[j]       <- k$p.value
    
  }
  
  replicates<- dk %>% filter(p>0.05) %>% tidyr::separate(pair, into = c("Rep_1","Rep_2"), "_") %>% 
    select(Station_ID, Rep_1,Rep_2) %>% 
    reshape2::melt("Station_ID") %>% 
    mutate(replicate=factor(value, c("A","B","C","D","E"))) %>% 
    group_by(Station_ID, replicate) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    head(n = 3)
  
  return(replicates)
})


tidy_data_filt_no_dens<- tidy_data%>% 
  mutate(Type="Total") %>% 
  left_join(bind_rows(ks_tests), by = "Station_ID") %>%
  filter(Replicate ==replicate) 


###########################################
#Import summary concentration
###########################################
#list all Experiment files
Summary_files <- list.files(path=path.NTA,
                            pattern = "_Summary.csv",
                            full.names = TRUE, recursive = F)


#generate a list of all experiments 
SummaryData_list <- lapply(Summary_files, function(x) {
  
  dat <- read.csv(x,header=TRUE, sep=",", dec = ".",
                  fileEncoding = "utf8", skip = 93, nrows = 1000) %>% 
    select(Bin.centre..nm.,Concentration..particles...ml.) %>% 
    mutate_all(as.numeric)
  
  # Add a column with the sampling time
  label <- unlist(strsplit(basename(x), split ="_"))
  
  dat<- dat %>% 
    mutate(StationID = label[1],
           Media = label[2],
           Dye = label[3],
           Dilution = label[4],
           Method = label[5],
           stamp = label[6]) %>% 
    separate(StationID, into=c("Station_name","Fraction"), sep ="-") %>% 
    separate(stamp, into =c("Flow_rate","Date","Stamp"), sep =" ") 
  return(dat)
})


#aggregate all the samples into a single dataframe
SummaryData_df_no_dens <- bind_rows(SummaryData_list) %>%
  mutate(Type="Total") %>% 
  mutate(Tech_run = paste0(Date,Stamp)) %>% 
  left_join(station_IDs, by = "Station_name", keep = FALSE) %>% 
  group_by(Region, Station_ID,Type, Fraction) %>%
  mutate(Replicate = factor(Tech_run, 
                            levels = unique(Tech_run),
                            labels = c("A","B","C","D","E")[1:length(unique(Tech_run))])) %>% 
  filter(Bin.centre..nm. < 250, Bin.centre..nm. > 20) %>% 
  group_by(Region, Station_ID, Type, Fraction, Replicate, Bin.centre..nm.) %>%
  summarize(Conc=Concentration..particles...ml./60) %>% #correct for x60 (vivaflow)
  group_by(Region, Station_ID, Type, Replicate, Bin.centre..nm.) %>%
  summarize(Part.conc=sum(Conc))#
