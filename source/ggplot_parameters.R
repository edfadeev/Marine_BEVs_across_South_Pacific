require(ggplot2)

# 8 colours
cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7")

#21 colours
tol21rainbow<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", "#DD7788")

#theme
theme_EF <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        title = element_text(size=25, face ="bold"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        strip.text = element_text(size = 22, face ="bold"),
        #axis
        axis.title=element_text(size=20,face="bold"),
        #axis.title.x=element_blank(),
        #axis.text.x = element_text(angle=90),
        axis.text=element_text(size=18,color="black"),
        axis.line = element_line(colour = 'black', linewidth = 0.8),
        strip.placement='outside',
        strip.background.x=element_blank(),
        panel.spacing.x=unit(5,"pt"))
