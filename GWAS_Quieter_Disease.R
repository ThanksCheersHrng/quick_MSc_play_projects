'GWAS Quieter Disease'

#Save Screen 
save.image("quiet_disease.R")

setwd("C:/Users/user/Desktop/PROJECT/PCA_FA_ICA_code")

load("quiet_disease.R")

#Packages to Load 

library(LEA)
library(caret)
library(psych)
library(stats)
library(fastICA)

##Data Setup 
#Load and Clean 

#test<-ped2geno("hapmap3_r2_b36_fwd.consensus.qc.poly.ped")
read1<-read.geno("hapmap3_r2_b36_fwd.consensus.qc.poly.geno")
#for some unknown reason, my LEA:read.geno() function is producing read1 as a list instead of a matrix! 
#very frustration
#the .geno file in wd() has also mysteriously shrunk from 1.6 million KB (last I checked) to 600 thousand KB today! 
#it does seem to have lost information, as well, reading in 1440616 items in the past, and recently, just 530572. 


#instead, I'm starting from .ped scratch: 
test<-ped2geno("hapmap3_r2_b36_fwd.consensus.qc.poly.ped")
'
	- number of detected individuals:	1184
	- number of detected loci:		1440616
'
#this refreshsed the .geno to its original 1.6 million KB. 
read1<-read.geno("hapmap3_r2_b36_fwd.consensus.qc.poly.geno")
#and read1 came in as a large matrix this time! 


max.col <-NA
for(i in 1:ncol(read1)) max.col[i]<-max(read1[,i])
clean1<-read1[,max.col==2]
#range(clean1) #confirms the range is now 0 to 2

#to save space, periodically remove now-cleaned data frames
rm(test, read1, max.col) #even though it took ages to get read1 back, 
#it also takes ages to save.image whenever read1 is still present. it's enormous!! 

var.col <- NA 
for (i in 1:ncol(clean1))
  var.col[i] = var(clean1[,i])
clean2 <- clean1[,var.col!=0]

rm(clean1, var.col)

#this time, let's select a subset of SNPs from the middle of the genome; 
#if these are ordered by location on chromosomes, we may have been working 
#on telomeres all this time! (this is why I had to re-load read1)
#due to their un-expressed-ness, I don't know if telomeres would 
#vary in a comparable manner to verbose genes. 

dim(clean2) # 1184 464817 # Ostensibly could easily save 45 different cleaned SNP datasets. 

corr.col <- findCorrelation(cor(clean2[,c(15000:25000)]), cutoff=0.9)
 
clean3 <- clean2[,c(15000:25000)]
clean3 <- clean3[,-corr.col]

save.image("quiet_disease.R")

corr.col <- findCorrelation(cor(clean2[,c(115000:125000)]), cutoff=0.9)

clean1 <- clean2[,c(115000:125000)]
clean1 <- clean1[,-corr.col] 

corr.col <- findCorrelation(cor(clean2[,c(65000:75000)]), cutoff=0.9)

clean4 <- clean2[,c(65000:75000)]
clean4 <- clean4[,-corr.col]

corr.col <- findCorrelation(cor(clean2[,c(285000:295000)]), cutoff=0.9)

clean2 <- clean2[,c(285000:295000)]
clean2 <- clean2[,-corr.col]
#save clean1 through clean2 and you have 4 cleaned data sets of <10,000 SNPs to play with. 

rm(corr.col)

#Randomly Assigning Disease Status 
#I've used clean 1 but I could equivalently replace with clean2,3, or 4 and get different results. 

#each time, these will be programmed with a 20% disease prevalence 
disease.index <- sample(1:nrow(clean1), nrow(clean1)/5, replace=FALSE)

disease.status <- rep(0,nrow(clean1))

for (i in disease.index)
  disease.status[i] = 1

mock_disease_data <- cbind(disease.status, clean1)

save.image("quiet_disease.R")


##Population Structure (Dimensionality Reduction)
#PCA        ### (on a further reduced SNP dataset)
pca7<-principal(mock_disease_data[,1:2000],nfactors=7,scores=TRUE)

pc1 <- pca7$scores[,1]
pc2 <- pca7$scores[,2]
pc3 <- pca7$scores[,3]
pc4 <- pca7$scores[,4]
pc5 <- pca7$scores[,5]
pc6 <- pca7$scores[,6]
pc7 <- pca7$scores[,7]


pca_sources <- data.frame(pc1, pc2, pc3, pc4, pc5,
                          pc6, pc7)

#FA 7 times since each q-factor model is separately generated...
'
fac1<-factanal(mock_disease_data[,1:2000],1,scores="regression")
fac2<-factanal(mock_disease_data[,1:2000],2,scores="regression")
fac3<-factanal(mock_disease_data[,1:2000],3,scores="regression")
fac4<-factanal(mock_disease_data[,1:2000],4,scores="regression")
fac5<-factanal(mock_disease_data[,1:2000],5,scores="regression")
fac6<-factanal(mock_disease_data[,1:2000],6,scores="regression")
fac7<-factanal(mock_disease_data[,1:2000],7,scores="regression")

fa_sources <- data.frame(fac1$scores,fac2$scores,
                         fac3$scores,fac4$scores,fac5$scores,
                         fac6$scores,fac7$scores)

names(fa_sources) <- c("fa1.1",
                       "fa2.1", "fa2.2",
                       "fa3.1","fa3.2", "fa3.3",
                       "fa4.1","fa4.2","fa4.3","fa4.4",
                       "fa5.1","fa5.2","fa5.3","fa5.4","fa5.5",
                       "fa6.1","fa6.2","fa6.3","fa6.4","fa6.5",
                       "fa6.6",
                       "fa7.1","fa7.2","fa7.3","fa7.4","fa7.5",
                       "fa7.6","fa7.7") 
'
#ICA 
ica<-fastICA(mock_disease_data[,1:2000],7)
#plot(test.ica$S, col = 6, pch = 17)

ica_sources <- data.frame(ica$S)
#names(ica_scores) #X5 instead of ic5 

save.image("quiet_disease.R")

##Create data frame with pop structure sources 

dat <- as.data.frame(mock_disease_data[,1:2000]) 
snp_names <- as.character(1:2000)
colnames(dat) <- snp_names

dis <- mock_disease_data[,1]


data_full <- cbind(dis,dat,pca_sources,#fa_sources,
                   ica_sources)
data_full <- as.data.frame(data_full)

save.image("quiet_disease.R")

'Single SNP GLMs with varying model structures'

#this function returns SNP coefficients for the model structure with no population structure
coef_glm0 <- function(x) coef(glm(dis ~ x, data = data_full, family = binomial))[2]

#apply for every column (MARGIN = 2 for columns) 
coef0_list <- apply(dat, MARGIN = 2, coef_glm0)

sum0 <- sum(abs(coef0_list))
sum0 

#PCA: 
coef_glm1.1 <- function(x) coef(glm(dis ~ x + pc1, 
                                    data = data_full, family = binomial))[2]
coef_glm1.2 <- function(x) coef(glm(dis ~ x + pc1+pc2, 
                                    data = data_full, family = binomial))[2]
coef_glm1.3 <- function(x) coef(glm(dis ~ x + pc1+pc2+pc3, 
                                    data = data_full, family = binomial))[2]
coef_glm1.4 <- function(x) coef(glm(dis ~ x + pc1+pc2+pc3+pc4, 
                                    data = data_full, family = binomial))[2]
coef_glm1.5 <- function(x) coef(glm(dis ~ x + pc1+pc2+pc3+pc4+pc5, 
                                    data = data_full, family = binomial))[2]
coef_glm1.6 <- function(x) coef(glm(dis ~ x + pc1+pc2+pc3+pc4+pc5+pc6, 
                                    data = data_full, family = binomial))[2]
coef_glm1.7 <- function(x) coef(glm(dis ~ x + pc1+pc2+pc3+pc4+pc5+pc6+pc7, 
                                    data = data_full, family = binomial))[2]
'
#FA: 
coef_glm2.1 <- function(x) coef(glm(dis ~ x + fa1.1, 
                                    data = data_full, family = binomial))[2]
coef_glm2.2 <- function(x) coef(glm(dis ~ x + fa2.1+fa2.2, 
                                    data = data_full, family = binomial))[2]
coef_glm2.3 <- function(x) coef(glm(dis ~ x + fa3.1+fa3.2+fa3.3, 
                                    data = data_full, family = binomial))[2]
coef_glm2.4 <- function(x) coef(glm(dis ~ x + fa4.1+fa4.2+fa4.3+fa4.4, 
                                    data = data_full, family = binomial))[2]
coef_glm2.5 <- function(x) coef(glm(dis ~ x + fa5.1+fa5.2+fa5.3+fa5.4+fa5.5, 
                                    data = data_full, family = binomial))[2]
coef_glm2.6 <- function(x) coef(glm(dis ~ x + fa6.1+fa6.2+fa6.3+fa6.4+fa6.5+fa6.6, 
                                    data = data_full, family = binomial))[2]
coef_glm2.7 <- function(x) coef(glm(dis ~ x + fa7.1+fa7.2+fa7.3+fa7.4+fa7.5+fa7.6+fa7.7, 
                                    data = data_full, family = binomial))[2]
'
#ICA: 
coef_glm3.1 <- function(x) coef(glm(dis ~ x + X1, 
                                    data = data_full, family = binomial))[2]
coef_glm3.2 <- function(x) coef(glm(dis ~ x + X1+X2, 
                                    data = data_full, family = binomial))[2]
coef_glm3.3 <- function(x) coef(glm(dis ~ x + X1+X2+X3, 
                                    data = data_full, family = binomial))[2]
coef_glm3.4 <- function(x) coef(glm(dis ~ x + X1+X2+X3+X4, 
                                    data = data_full, family = binomial))[2]
coef_glm3.5 <- function(x) coef(glm(dis ~ x + X1+X2+X3+X4+X5, 
                                    data = data_full, family = binomial))[2]
coef_glm3.6 <- function(x) coef(glm(dis ~ x + X1+X2+X3+X4+X5+X6, 
                                    data = data_full, family = binomial))[2]
coef_glm3.7 <- function(x) coef(glm(dis ~ x + X1+X2+X3+X4+X5+X6+X7, 
                                    data = data_full, family = binomial))[2]


#Apply and get sums for each model structure! 
#PCA
coef1.1_list <- apply(dat, MARGIN = 2, coef_glm1.1)
sum1.1 <- sum(abs(coef1.1_list))
coef1.2_list <- apply(dat, MARGIN = 2, coef_glm1.2)
sum1.2 <- sum(abs(coef1.2_list))
coef1.3_list <- apply(dat, MARGIN = 2, coef_glm1.3)
sum1.3 <- sum(abs(coef1.3_list))
coef1.4_list <- apply(dat, MARGIN = 2, coef_glm1.4)
sum1.4 <- sum(abs(coef1.4_list))
coef1.5_list <- apply(dat, MARGIN = 2, coef_glm1.5)
sum1.5 <- sum(abs(coef1.5_list))
coef1.6_list <- apply(dat, MARGIN = 2, coef_glm1.6)
sum1.6 <- sum(abs(coef1.6_list))
coef1.7_list <- apply(dat, MARGIN = 2, coef_glm1.7)
sum1.7 <- sum(abs(coef1.7_list))

save.image("quiet_disease.R")

'
#FA
coef2.1_list <- apply(dat, MARGIN = 2, coef_glm2.1)
sum2.1 <- sum(abs(coef2.1_list))
coef2.2_list <- apply(dat, MARGIN = 2, coef_glm2.2)
sum2.2 <- sum(abs(coef2.2_list))
coef2.3_list <- apply(dat, MARGIN = 2, coef_glm2.3)
sum2.3 <- sum(abs(coef2.3_list))
coef2.4_list <- apply(dat, MARGIN = 2, coef_glm2.4)
sum2.4 <- sum(abs(coef2.4_list))
coef2.5_list <- apply(dat, MARGIN = 2, coef_glm2.5)
sum2.5 <- sum(abs(coef2.5_list))
coef2.6_list <- apply(dat, MARGIN = 2, coef_glm2.6)
sum2.6 <- sum(abs(coef2.6_list))
coef2.7_list <- apply(dat, MARGIN = 2, coef_glm2.7)
sum2.7 <- sum(abs(coef2.7_list))

save.image("quiet_disease.R")
'

#ICA
coef3.1_list <- apply(dat, MARGIN = 2, coef_glm3.1)
sum3.1 <- sum(abs(coef3.1_list))
coef3.2_list <- apply(dat, MARGIN = 2, coef_glm3.2)
sum3.2 <- sum(abs(coef3.2_list))
coef3.3_list <- apply(dat, MARGIN = 2, coef_glm3.3)
sum3.3 <- sum(abs(coef3.3_list))
coef3.4_list <- apply(dat, MARGIN = 2, coef_glm3.4)
sum3.4 <- sum(abs(coef3.4_list))
coef3.5_list <- apply(dat, MARGIN = 2, coef_glm3.5)
sum3.5 <- sum(abs(coef3.5_list))
coef3.6_list <- apply(dat, MARGIN = 2, coef_glm3.6)
sum3.6 <- sum(abs(coef3.6_list))
coef3.7_list <- apply(dat, MARGIN = 2, coef_glm3.7)
sum3.7 <- sum(abs(coef3.7_list))

save.image("quiet_disease.R")

#mini table

n_sources <- c(1, 2, 3, 4,5,6,7)
PC_sums <- c(sum1.1,sum1.2,sum1.3,sum1.4,sum1.5,sum1.6,sum1.7)
#FA_sums <- c(sum2.1, sum2.2, sum2.3, sum2.4,sum2.5, sum2.6,sum2.7)
IC_sums <- c(sum3.1, sum3.2, sum3.3, sum3.4,sum3.5,sum3.6,sum3.7)

display_results <- data.frame(n_population_structure_terms = n_sources, 
                              PC_SNP_coef_total = PC_sums,
                              #FA_SNP_coef_total = FA_sums,
                              IC_SNP_coef_total = IC_sums)
sum0
display_results 


###A quick visual comparison of coefficient sums

bar_stac <- data.frame(sources=rep(c("PCA",#"FA",
                                     "ICA"), each = 7),
                       n_source=rep(n_sources, times = 2),
                       sum = c(PC_sums, #FA_sums, 
                               IC_sums))

ggplot(bar_stac, aes(col=sources, y=sum, x=n_source)) + 
  geom_point(shape = 18, size = 5) +
  geom_line(aes(group = sources))+
  #geom_smooth(aes(group = sources), method = "loess", se = FALSE)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))

