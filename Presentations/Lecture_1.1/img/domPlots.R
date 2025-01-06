


codom_phen <- c(10, 5, 0)
dom_phen <- c(10, 10, 0)
gen <- c("AA","Aa","aa")

df<-data.frame( Codominant = codom_phen,
                Dominant = dom_phen,
                genotype = gen)
library(reshape2)

df_m <- melt(df, id = "genotype")
df_m$genotype <- factor(df_m$genotype,
                        levels = c("AA","Aa","aa")) 

library(ggplot2)

ggplot(data= df_m[df_m$variable=="Codominant",], aes(x = genotype, 
                       y = value, 
                       group = variable,
                       col = variable,
                       shape = variable))+
  geom_line(lwd=1)+
  geom_point(size =8)+
  ylab("Phenotype")+
  xlab("Genotype")+
  scale_colour_manual(values = c("#E69F00", "#56B4E9"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 10),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold"))



ggplot(data= df_m, aes(x = genotype, 
                       y = value, 
                       group = variable,
                       col = variable,
                       shape = variable))+
  geom_line(lwd=1)+
  geom_point(size =8)+
  ylab("Phenotype")+
  xlab("Genotype")+
  scale_colour_manual(values = c("#E69F00", "#56B4E9"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 10),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold"))


df$arbitrary_1 <- c(0, 7.5, 10)
df$arbitrary_2 <- c(0, 1.5, 10)

df_m <- melt(df, id = "genotype")
df_m$genotype <- factor(df_m$genotype,
                        levels = c("AA","Aa","aa")) 
df_m$variable <- factor(df_m$variable,
                        levels = c("Codominant","Dominant","arbitrary_1", "arbitrary_2"),
                        labels = c("Codominant","Dominant","Partial Dominance", "Partial Recessivity")) 

colorBlindPalette   <- c( "#E69F00", "#56B4E9", "#009E73", 
                          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(ggplot2)


ggplot(data= df_m, aes(x = genotype, 
                       y = value, 
                       group = variable,
                       col = variable,
                       shape = variable))+
  geom_line(lwd=1)+
  geom_point(size =8)+
  ylab("Phenotype")+
  xlab("Genotype")+
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#D55E00","#009E73"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 10),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold"))

