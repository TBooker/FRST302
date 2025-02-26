## Make pedigree diagrams for FRST302

library(kinship2)

data(sample.ped)
sample.ped[1:10,]

datped2 <- read.csv("~/UBC/Teaching/FRST302/Presentations/Lecture_4.2/demoPed.csv")

tryout <- try({
  ped2 <- with(pedCoded, pedigree(id, father, mother, sex))
})
fixped2 <- with(datped2, fixParents(id, father, mother, sex))
ped2 <- with(fixped2, pedigree(id, dadid, momid, sex))
plot(ped2)

 kin2 <- kinship(ped2)
heatmap(kin2)

df <- as.data.frame(as.table(kin2))
 
library(ggplot2)

ggplot(data = df, aes(x = Var2, y = Var1, fill = Freq*2))+
  geom_tile()
library(ASRgenomics)
kinship.heatmap(kin2)
 

library(pedSimulate)
nSNP = 10

hist(rbeta(1000, 0.1,0.1))
AF = rbeta(1000, 0.1,0.1)

mut.rate = runif(nSNP, 0, 10^-5)
ped = data.frame(ID=1:5, SIRE=c(0,0,1,1,3), DAM=c(0,0,2,2,4))

pedCoded = data.frame(ped=1, id=1:5, father=c(0,0,1,1,3), mother=c(0,0,2,2,4), sex = c(1,2,1,2,1), affected=0,avail=0)

gen = simulateGen(ped, AF, mut.rate)
n=nrow(gen)
x= as.matrix(
  as.data.frame(gen)[(colSums(as.data.frame(gen)) != 0)&(colSums(as.data.frame(gen)) != 2*n)])




p_hat = apply(x, 2, sum)/(2*n)
w = apply(rbind(x,p_hat), 2, function(x) (x-2*x[length(x)])/sqrt(2*x[length(x)]*(1-x[length(x)])))[1:n,]

A = w %*% t(w) /m
heatmap(A)

n=5; m=3
set.seed(10)
p = runif(m, min=0.2, max=0.5) ### allele frequency were draw from a uniform distribution
x_A1 = t(replicate(n, rbinom(m, 1, p)))   ### for the plink ped file
x_A2 = t(replicate(n, rbinom(m, 1, p)))
x = x_A1 + x_A2
colnames(x) = paste("SNP",1:m,sep="")
rownames(x) = paste("indi", 1:n, sep="")
x
