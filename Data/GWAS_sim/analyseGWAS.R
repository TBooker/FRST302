rm(list= ls())

phen = read.csv("~/UBC/Teaching/FRST302/Data/GWAS_sim/GWAS/gwasDemo.fam", sep = " ", header = F)
names(phen)<-c("id","id2","dummy","dummy2","dummy3","phen")
phen$pop <- c(rep("p1", nrow(phen)/2),rep("p2", nrow(phen)/2))


library(ggplot2)
phen$Island <- factor(phen$pop, levels = c("p1", "p2"), labels = c("1","2"))
min(phen$phen)
phenPlot <- ggplot(data = phen, aes(x = phen+10/3, fill = Island))+ ## Adjust the phenotypes so there are no negative heights on the plot
  geom_density(alpha = 0.7)+
  xlab("Height (m)")+
  ylab("Density")+
  scale_fill_manual("Island", values = c("#E69F00", "#56B4E9"))+
  xlim(c(0,25))+
  theme_bw()

ggsave("~/UBC/Teaching/FRST302/Presentations/Lecture_4.1/img/treePhens.png",
       plot = phenPlot,
       width = 6,
       height = 5)

t.test(phen$phen+10/3~phen$pop)


snp <- read.csv("~/UBC/Teaching/FRST302/Data/GWAS_sim/GWAS_temp//singleSNP.csv")
sum(snp$alleleCount)/2000


snp_mod <- lm(snp$phen~snp$alleleCount) 

ggplot(data=snp, aes(x=alleleCount, y = phen))+
  geom_jitter(height = FALSE, width = 0.1)+
  scale_y_continuous("Height (m)")+
  scale_x_continuous("Number of Copies of the Allele",
                     breaks = c(0,1,2))+
  theme_bw()


ggplot(data=snp, aes(x=alleleCount, y = phen))+
  geom_jitter(height = FALSE, width = 0.1)+
  scale_y_continuous("Height (m)")+
  geom_abline(slope = coef(snp_mod)[["snp$alleleCount"]], 
              intercept = coef(snp_mod)[["(Intercept)"]],
              col = "red", lwd= 1, lty = 2)+
  scale_x_continuous("Number of Copies of the Allele",
                     breaks = c(0,1,2))+
  theme_bw()
summary(snp_mod)

#Exported jpegs with width=800, height=350



## Now do that simple model on all SNPs

library(vcfR)
vcf<- read.vcfR("~/UBC/Teaching/FRST302/Data/GWAS_sim/GWAS/gwasDemo.vcf.gz")
gtMat <- extract.gt(vcf, element = 'GT', as.numeric = TRUE)

manual_gwas <- list()
for (snp in 1:nrow(gtMat)){
  snp_mod <- lm(phen$phen ~ gtMat[snp,])
  # Give the slope estimate and the p-value as output of the model
  snp_mod_out <-data.frame(snp=rownames(gtMat)[snp], 
                           simpleEst = snp_mod$coefficients[2],
                           uncor_pVal=summary(snp_mod)$coefficients[,4][2])
  manual_gwas[[snp]] = snp_mod_out
}

simpleGWAS <- do.call(rbind, manual_gwas)



## Now read in the GEMMA GWAS
gwasRes <- read.csv("~/UBC/Teaching/FRST302/Data/GWAS_sim/GWAS/output/gwasDemo.assoc.txt", sep = "\t")
names(gwasRes)[3] <- "position"

gwasRes[gwasRes$position==10000000,]
snpDat <- read.csv("~/UBC/Teaching/FRST302/Data/GWAS_sim/GWAS/gwasDemo.muts.csv")
snpPosData <- merge(gwasRes, snpDat, on ="position", all.x=T)


snpPosDataGWAS <- merge(gwasRes, snpPosData, on = "id", all.x=T)
snpPosDataGWAS$effect[is.na(snpPosDataGWAS$effect)] <- 0
library(dplyr)
snpPosDataGWAS[snpPosDataGWAS$effect^2>0.01,]
snpPosDataGWAS <- snpPosDataGWAS %>% mutate(effect = ifelse(is.na(effect), 0, effect))

names(snpPosDataGWAS)[2] <- "snp"

bothGWAS <- merge(snpPosDataGWAS, simpleGWAS, on = "snp", all.x=T)
justTruePos <- bothGWAS[abs(bothGWAS$effect)>0.0,]


hist(bothGWAS$uncor_pVal)
hist(bothGWAS$p_lrt)

unCorPlot<- ggplot(data = bothGWAS, aes(x= position/1e6, y= -log10(uncor_pVal)))+
 # geom_point()+
 # geom_point(data = justTruePos, aes(x =position/1e6, y = -log10(uncor_pVal)), col = "red",size = 3)+
  xlab("Position in Chromosome (Mbp)")+ 
 #geom_hline(yintercept = -log10(0.05/nrow(gtMat)))+
  ylab(expression(-log[10]*"(p-value)"))+
  theme_bw()

ggsave("~/UBC/Teaching/FRST302/Presentations/Lecture_4.1/img/uncorPlot_noData.png",
       plot = unCorPlot,
       width = 6,
       height = 3)

unCorPlot<- unCorPlot + geom_point()
  
ggsave("~/UBC/Teaching/FRST302/Presentations/Lecture_4.1/img//uncorPlot_Data.png",
       plot = unCorPlot,
       width = 6,
       height = 3)

unCorPlot<- unCorPlot + geom_hline(yintercept = -log10(0.05),
                                   lty = 2, col = "#D55E00")


ggsave("~/UBC/Teaching/FRST302/Presentations/Lecture_4.1/img//uncorPlot_DataHline.png",
       plot = unCorPlot,
       width = 6,
       height = 3)


unCorPlot<- unCorPlot +  geom_point(data = justTruePos, aes(x =position/1e6, y = -log10(uncor_pVal)), fill = "#009E73",shape=21,size = 3)



ggsave("~/UBC/Teaching/FRST302/Presentations/Lecture_4.1/img//uncorPlot_DataHlineTrue.png",
       plot = unCorPlot,
       width = 6,
       height = 3)




corPlot <- ggplot(data = bothGWAS, aes(x= position/1e6, y= -log10(p_lrt)))+
  geom_point()+
  geom_point(data = justTruePos, aes(x =position/1e6, y = -log10(p_lrt)), fill = "#009E73",shape=21,size = 3)+
  xlab("Position in Chromosome (Mbp)")+
  ylab(expression(-log[10]*"(p-value)"))+
  theme_bw()

ggsave("~/UBC/Teaching/FRST302/Presentations/Lecture_4.1/img/corPlot_unCorrectedLine.png",
       plot = corPlot,
       width = 6,
       height = 3)

corPlot_adj <-  corPlot +  geom_hline(yintercept = -log10(0.05/nrow(bothGWAS)), col = "orange", lty=2, lwd=2)

ggsave("~/UBC/Teaching/FRST302/Presentations/Lecture_4.1/img/corPlot_correctedLine.png",
       plot = unCorPlot,
       width = 6,
       height = 3)

justfuncMuts <- snpPosDataGWAS[snpPosDataGWAS$effect!=0,]
justSigMuts <- justfuncMuts[justfuncMuts$p_wald<0.05,]

ggplot(data= justSigMuts, aes(x=beta^2/effect^2))+
  geom_density(bins =30)+
  geom_vline(xintercept=1)




mat <- scan('~/UBC/Teaching/FRST302/Data/GWAS_sim/GWAS/output/gwasDemo.cXX.txt')
mat <- matrix(mat, ncol = 1000, byrow = TRUE)

relMat.df <- reshape2::melt(mat, c("x", "y"), value.name = "z")

ggplot(data=relMat.df,aes(x=x,y=y,fill=z))+
  geom_tile()


mat2 <- matrix(mat[900:1000,1:100], ncol = 100, byrow = TRUE)

relMat.df <- reshape2::melt(mat2, c("x", "y"), value.name = "z")



ggplot(data=relMat.df,aes(x=x,y=y,fill=z))+
  geom_tile()
