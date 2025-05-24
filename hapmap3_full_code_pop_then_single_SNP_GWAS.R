'GWAS Main1'

#Save Screen 
save.image("hapmap3_full_code.R")

setwd("C:/Users/user/Desktop/PROJECT/PCA_FA_ICA_code")

load("hapmap3_full_code.R")

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

max.col <-NA
for(i in 1:ncol(read1)) max.col[i]<-max(read1[,i])
clean1<-read1[,max.col==2]
range(clean1) #confirms the range is now 0 to 2

#to save space, periodically remove now-cleaned data frames
rm(test, read1, max.col)

var.col <- NA 
for (i in 1:ncol(clean1))
  var.col[i] = var(clean1[,i])
clean2 <- clean1[,var.col!=0]

rm(clean1, var.col)

corr.col <- findCorrelation(cor(clean2[,c(1:12000)]), cutoff=0.9)

clean3 <- clean2[,c(1:12000)]
clean3 <- clean3[,-corr.col]

rm(clean2, corr.col,i)

#Randomly Assigning Disease Status 

disease.index <- sample(1:nrow(clean3), nrow(clean3)/2, replace=FALSE)

disease.status <- rep(0,nrow(clean3))

for (i in disease.index)
  disease.status[i] = 1

mock_disease_data <- cbind(disease.status, clean3)

mock_disease_data[1:10,1:5] #just to check it looks right. 
dim(mock_disease_data) #1184 samples, one affection status variable and 9170 SNPs remaining after cleaning
sum(mock_disease_data[,1]==1) #592 out of 1184 is roughly half the population


rm(clean3, disease.index, disease.status, i)


##### Global Environment now contains
#####  disease.status (length = 89) and mock_disease_data (dim = 89 by 10127)

save.image("hapmap3_full_code.R")


#Dim redux on whole pop 

corr_matrix <- cor(mock_disease_data[,-1])

##Population Structure (Dimensionality Reduction)
#PCA 
pca25<-principal(mock_disease_data[,2:1150],nfactors=25,scores=TRUE)
#plot(test.pca$scores[,1],test.pca$scores[,2], col = 4, pch = 15)

pc1 <- pca25$scores[,1]
pc2 <- pca25$scores[,2]
pc3 <- pca25$scores[,3]
pc4 <- pca25$scores[,4]
pc5 <- pca25$scores[,5]
pc6 <- pca25$scores[,6]
pc7 <- pca25$scores[,7]
pc8 <- pca25$scores[,8]
pc9 <- pca25$scores[,9]
pc10 <- pca25$scores[,10]
pc11 <- pca25$scores[,11]
pc12 <- pca25$scores[,12]
pc13 <- pca25$scores[,13]
pc14 <- pca25$scores[,14]
pc15 <- pca25$scores[,15]
pc16 <- pca25$scores[,16]
pc17 <- pca25$scores[,17]
pc18 <- pca25$scores[,18]
pc19 <- pca25$scores[,19]
pc20 <- pca25$scores[,20]
pc21 <- pca25$scores[,21]
pc22 <- pca25$scores[,22]
pc23 <- pca25$scores[,23]
pc24 <- pca25$scores[,24]
pc25 <- pca25$scores[,25]

pca_sources <- data.frame(pc1, pc2, pc3, pc4, pc5,
                          pc6, pc7, pc8, pc9, pc10,
                          pc11, pc12, pc13, pc14, pc15,
                          pc16, pc17, pc18, pc19, pc20,
                          pc21, pc22, pc23, pc24, pc25)

#FA #since each q-factor model is separately generated...
fac1<-factanal(mock_disease_data[,2:1150],1,scores="regression")
fac2<-factanal(mock_disease_data[,2:1150],2,scores="regression")
fac3<-factanal(mock_disease_data[,2:1150],3,scores="regression")
fac4<-factanal(mock_disease_data[,2:1150],4,scores="regression")
fac5<-factanal(mock_disease_data[,2:1150],5,scores="regression")
fac6<-factanal(mock_disease_data[,2:1150],6,scores="regression")
  # honestly it became painful watching the computer try to run EVERY factor analysis, 
  # so I skipped a few at this stage....
#fac7<-factanal(mock_disease_data[,2:1150],7,scores="regression")
fac8<-factanal(mock_disease_data[,2:1150],8,scores="regression")
#fac9<-factanal(mock_disease_data[,2:1150],9,scores="regression")
fac10<-factanal(mock_disease_data[,2:1150],10,scores="regression")
#fac11<-factanal(mock_disease_data[,2:1150],11,scores="regression")
#fac12<-factanal(mock_disease_data[,2:1150],12,scores="regression")
#fac13<-factanal(mock_disease_data[,2:1150],13,scores="regression")
#fac14<-factanal(mock_disease_data[,2:1150],14,scores="regression")
#fac15<-factanal(mock_disease_data[,2:1150],15,scores="regression")
fac16<-factanal(mock_disease_data[,2:1150],16,scores="regression")
#fac17<-factanal(mock_disease_data[,2:1150],17,scores="regression")
#fac18<-factanal(mock_disease_data[,2:1150],18,scores="regression")
#fac19<-factanal(mock_disease_data[,2:1150],19,scores="regression")
#fac20<-factanal(mock_disease_data[,2:1150],20,scores="regression")
##fac21<-factanal(mock_disease_data[,2:1150],21,scores="regression")
#fac22<-factanal(mock_disease_data[,2:1150],22,scores="regression")
#fac23<-factanal(mock_disease_data[,2:1150],23,scores="regression")
#fac24<-factanal(mock_disease_data[,2:1150],24,scores="regression")
##fac25<-factanal(mock_disease_data[,2:1150],25,scores="regression")


fa_sources <- data.frame(fac1$scores,fac2$scores,fac3$scores,fac4$scores,fac5$scores,
                         fac6$scores,#fac7$scores,
                         fac8$scores,#fac9$scores,
                         fac10$scores,
                         #fac11$scores,fac12$scores,fac13$scores,fac14$scores,fac15$scores,
                         fac16$scores#,#fac17$scores,fac18$scores,fac19$scores,fac20$scores,
                        # fac21$scores,#fac22$scores,fac23$scores,fac24$scores,
                        # fac25$scores
                        )
names(fa_sources) <- c("fa1.1",
                       "fa2.1", "fa2.2",
                       "fa3.1","fa3.2", "fa3.3",
                       "fa4.1","fa4.2","fa4.3","fa4.4",
                       "fa5.1","fa5.2","fa5.3","fa5.4","fa5.5",
                       "fa6.1","fa6.2","fa6.3","fa6.4","fa6.5",
                       "fa6.6",
                       #"fa7.1","fa7.2","fa7.3","fa7.4","fa7.5",
                       #"fa7.6","fa7.7",
                       "fa8.1","fa8.2","fa8.3","fa8.4","fa8.5",
                       "fa8.6","fa8.7","fa8.8",
                       #"fa9.1","fa9.2","fa9.3","fa9.4","fa9.5",
                       #"fa9.6","fa9.7","fa9.8","fa9.9",
                       "fa10.1","fa10.2","fa10.3","fa10.4","fa10.5",
                       "fa10.6","fa10.7","fa10.8","fa10.9","fa10.10",
                       #
                       "fa16.1","fa16.2","fa16.3","fa16.4","fa16.5",
                       "fa16.6","fa16.7","fa16.8","fa16.9","fa16.10",
                       "fa16.11","fa16.12","fa16.13","fa16.14","fa16.15",
                       "fa16.16"#,
                       #
                      # "fa21.1","fa21.2","fa21.3","fa21.4","fa21.5",
                      # "fa21.6","fa21.7","fa21.8","fa21.9","fa21.10",
                      # "fa21.11","fa21.12","fa21.13","fa21.14","fa21.15",
                      # "fa21.16","fa21.17","fa21.18","fa21.19","fa21.20",
                      # "fa21.21",
                       #
                      # "fa25.1","fa25.2","fa25.3","fa25.4","fa25.5",
                      # "fa25.6","fa25.7","fa25.8","fa25.9","fa25.10",
                      # "fa25.11","fa25.12","fa25.13","fa25.14","fa25.15",
                      # "fa25.16","fa25.17","fa25.18","fa25.19","fa25.20",
                      # "fa25.21","fa25.22","fa25.23","fa25.24","fa25.25"
                      ) 

#ICA 
ica<-fastICA(mock_disease_data[,2:1150],25)
#plot(test.ica$S, col = 6, pch = 17)

ica_sources <- data.frame(ica$S)
names(ica_scores) #so now I know to refer to ic's as X, i.e. X9 is the 9th independent component. 


save.image("hapmap3_full_code.R")

##Create data frame with pop structure sources 
#dim redux has only occurred for the first 1149 SNPs. 
#this experiment could be conveniently repeated for the next 1100 or so.

dat <- as.data.frame(mock_disease_data[,2:1150]) #this is the first 1149 SNPs, no disease (i.e. inputs)
snp_names <- as.character(1:1149)
colnames(dat) <- snp_names
fix(dat) #can use this to quickly check if column names are in place. 

dis <- mock_disease_data[,1]


###In order to run glm on a single data frame, 
'I need to add the disease status and scores (above) 
to the data frame called dat, with appropriate column names.'

data_full <- cbind(dis,dat,pca_sources,fa_sources, ica_sources)
data_full <- as.data.frame(data_full)


save.image("hapmap3_full_code.R")

'Single SNP GLMs with varying model structures'

#Model Structure 0 : no corrective components 
glm0.1<-glm(dis ~ dat$"3", data = data_full, family = binomial)

glm(dis ~ dat$"3", data = data_full, family = binomial)$coefficients[2] 
glm(dis ~ dat$"4", data = data_full, family = binomial)$coefficients[2]

scoef1 <- coef(glm(dis ~ dat$"1", data = data_full, family = binomial))[2] 
scoef2 <- coef(glm(dis ~ dat$"2", data = data_full, family = binomial))[2]

coef_SNPs.0 <- c(scoef1,scoef2)

sum(coef_SNPs.0)


glm0_func <- function(x) glm(dis ~ x, data = data_full, family = binomial)

glm0_func(dat$"3") #works 

####start here 

#this function returns SNP coefficients for the model structure with no population structure
coef_glm0 <- function(x) coef(glm(dis ~ x, data = data_full, family = binomial))[2]

coef_glm0(dat$"3") 

#apply for every column (MARGIN = 2 for columns) 
coef0_list <- apply(dat, MARGIN = 2, coef_glm0)

sum0 <- sum(abs(coef0_list))
#sum0 will give sum of coefficients of SNPs when not controlled with any pop struc (i.e. glm0)

#model structure naming system: 
#coef_glmA.B, coefA.B_list, and sumA.B 
#A is 1 for PCA, 2 for FA, and 3 for ICA 
#B is the number of added components/factors

#just recall what those term names are:
names(data_full[,1185:1255])
#faSIZE.FACTOR
#X12 = ic12 
names(data_full[,1150:1185])
#pcCOMPONENT

##First I'll set up each function, for model structures where B=1, 2, 3, 4, 5, 6, 8, 10, or 16

#A=1 for PCA: 
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
coef_glm1.8 <- function(x) coef(glm(dis ~ x + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8, 
                                    data = data_full, family = binomial))[2]
coef_glm1.10 <- function(x) coef(glm(dis ~ x + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, 
                                    data = data_full, family = binomial))[2]
coef_glm1.16 <- function(x) coef(glm(dis ~ x + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16, 
                                    data = data_full, family = binomial))[2]

#A=2 for FA: 
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
coef_glm2.8 <- function(x) coef(glm(dis ~ x + fa8.1+fa8.2+fa8.3+fa8.4+fa8.5+fa8.6+fa8.7+fa8.8, 
                                    data = data_full, family = binomial))[2]
coef_glm2.10 <- function(x) coef(glm(dis ~ x + fa10.1+fa10.2+fa10.3+fa10.4+fa10.5+fa10.6+fa10.7+fa10.8+fa10.9+fa10.10, 
                                     data = data_full, family = binomial))[2]
coef_glm2.16 <- function(x) coef(glm(dis ~ x + fa16.1+fa16.2+fa16.3+fa16.4+fa16.5+fa16.6+fa16.7+fa16.8+fa16.9+fa16.10+fa16.11+fa16.12+fa16.13+fa16.14+fa16.15+fa16.16, 
                                     data = data_full, family = binomial))[2]

#A=3 for ICA: 
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
coef_glm3.8 <- function(x) coef(glm(dis ~ x + X1+X2+X3+X4+X5+X6+X7+X8, 
                                    data = data_full, family = binomial))[2]
coef_glm3.10 <- function(x) coef(glm(dis ~ x + X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, 
                                     data = data_full, family = binomial))[2]
coef_glm3.16 <- function(x) coef(glm(dis ~ x + X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16, 
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
coef1.8_list <- apply(dat, MARGIN = 2, coef_glm1.8)
sum1.8 <- sum(abs(coef1.8_list))
coef1.10_list <- apply(dat, MARGIN = 2, coef_glm1.10)
sum1.10 <- sum(abs(coef1.10_list))
coef1.16_list <- apply(dat, MARGIN = 2, coef_glm1.16)
sum1.16 <- sum(abs(coef1.16_list))

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
coef2.8_list <- apply(dat, MARGIN = 2, coef_glm2.8)
sum2.8 <- sum(abs(coef2.8_list))
coef2.10_list <- apply(dat, MARGIN = 2, coef_glm2.10)
sum2.10 <- sum(abs(coef2.10_list))
coef2.16_list <- apply(dat, MARGIN = 2, coef_glm2.16)
sum2.16 <- sum(abs(coef2.16_list))

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
coef3.8_list <- apply(dat, MARGIN = 2, coef_glm3.8)
sum3.8 <- sum(abs(coef3.8_list))
coef3.10_list <- apply(dat, MARGIN = 2, coef_glm3.10)
sum3.10 <- sum(abs(coef3.10_list))
coef3.16_list <- apply(dat, MARGIN = 2, coef_glm3.16)
sum3.16 <- sum(abs(coef3.16_list))


#mini table- we should see the SNP coef 
#sums descend with the addition of each pc 

n_sources <- c(1, 2, 3, 4,5,6,8,10,16)
PC_sums <- c(sum1.1,sum1.2,sum1.3,sum1.4,sum1.5,sum1.6,sum1.8,sum1.10,sum1.16)
FA_sums <- c(sum2.1, sum2.2, sum2.3, sum2.4,sum2.5, sum2.6,sum2.8,sum2.10,sum2.16)
IC_sums <- c(sum3.1, sum3.2, sum3.3, sum3.4,sum3.5,sum3.6,sum3.8,sum3.10,sum3.16)

display_results <- data.frame(n_population_structure_terms = n_sources, 
                              PC_SNP_coef_total = PC_sums,
                              FA_SNP_coef_total = FA_sums,
                              IC_SNP_coef_total = IC_sums)
display_results 


###A quick visual comparison of coefficient sums

bar_stac <- data.frame(sources=rep(c("PCA","FA","ICA"), each = 9),
                       n_source=rep(n_sources, times = 3),
                       sum = c(PC_sums, FA_sums, IC_sums))

ggplot(bar_stac, aes(col=sources, y=sum, x=n_source)) + 
  geom_point(shape = 18, size = 5) +
  geom_line(aes(group = sources))+
  #geom_smooth(aes(group = sources), method = "loess", se = FALSE)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))






#GLM with PC's
#1pc
coef_glm1.1 <- function(x) coef(glm(dis ~ x + pc1, 
                                  data = data_full, family = binomial))[2]

#What if the above function doesn't 'line up' correctly? 
coef_glm1.1(dat$"3")
coef(glm(dis ~ dat$"3" +pc1, data_full, family = binomial))
#it does line up with coef[2] as long as the SNP is the first predictor in the model 

coef1.1_list <- apply(dat, MARGIN = 2, coef_glm1.1)

sum1.1 <- sum(abs(coef1.1_list))

#2pc's
coef_glm1.2 <- function(x) coef(glm(dis ~ x + pc1 + pc2, 
                                    data = data_full, family = binomial))[2]

coef1.2_list <- apply(dat, MARGIN = 2, coef_glm1.2)

sum1.2 <- sum(abs(coef1.2_list))

#3pc's 
coef_glm1.3 <- function(x) coef(glm(dis ~ x + pc1 + pc2 + pc3,
                                   data = data_full, family = binomial))[2]
coef1.3_list <- apply(dat, MARGIN =2, coef_glm1.3)

sum1.3 <- sum(abs(coef1.3_list))

#4pc's 
coef_glm1.4 <- function(x) coef(glm(dis ~ x + pc1 + pc2 + pc3 + pc4,
                                   data = data_full, family = binomial))[2]
coef1.4_list <- apply(dat, MARGIN =2, coef_glm1.4)

sum1.4 <- sum(abs(coef1.4_list))

#mini table- we should see the SNP coef 
#sums descend with the addition of each pc 

PC_names <- c("no pcs", "1 pc", "2 pc", "3 pc", "4 pc")
PC_sums <- c(sum0, sum1.1,sum1.2,sum1.3,sum1.4)

quick_check_pc <- data.frame(model_structre = PC_names, 
                             sum_of_SNP_coefficients = PC_sums)
quick_check_pc 
#concerningly(?) it's not showing any marked improvement, or even pattern.



#GLM with F's 
#1factor 
coef_glm2.1 <- function(x) coef(glm(dis ~ x + fa1.1,
                                    data = data_full, family = binomial))[2]
coef2.1_list <- apply(dat, MARGIN =2, coef_glm2.1)

sum2.1 <- sum(abs(coef2.1_list))

#2factor 
coef_glm2.2 <- function(x) coef(glm(dis ~ x + fa2.1 + fa2.2,
                                    data = data_full, family = binomial))[2]
coef2.2_list <- apply(dat, MARGIN =2, coef_glm2.2)

sum2.2 <- sum(abs(coef2.2_list))

#3factor 
coef_glm2.3 <- function(x) coef(glm(dis ~ x + fa3.1 + fa3.2 + fa3.3,
                                    data = data_full, family = binomial))[2]
coef2.3_list <- apply(dat, MARGIN =2, coef_glm2.3)

sum2.3 <- sum(abs(coef2.3_list))

#4factor 
coef_glm2.4 <- function(x) coef(glm(dis ~ x + fa4.1 + fa4.2 + fa4.3 + fa4.4,
                                    data = data_full, family = binomial))[2]
coef2.4_list <- apply(dat, MARGIN =2, coef_glm2.4)

sum2.4 <- sum(abs(coef2.4_list))



#GLM with IC's
#1ic
coef_glm3.1 <- function(x) coef(glm(dis ~ x + ic1,
                                    data = data_full, family = binomial))[2]
coef3.1_list <- apply(dat, MARGIN =2, coef_glm3.1)

sum3.1 <- sum(abs(coef3.1_list))

#### As the results were unexpectedly pattern-less, 
#### I decided to double check that each function was 
### actually producing the desired coefficient: 
'
> coef_glm3.1(dat$"7")
x 
0.1511531 
> coef(glm(dis ~dat$"7"+ic1, data=data_full, family = binomial))
(Intercept)     dat$"7"         ic1 
-0.2976677   0.1511531  -0.1377120 
'
#Indeed they were. Therefore my best conclusion was that 
#I would need model structures featuring 
#more population structure sources 
#before I would see a downward trend in the sum of SNP coefficients. 

#2ic
coef_glm3.2 <- function(x) coef(glm(dis ~ x + ic1 + ic2,
                                    data = data_full, family = binomial))[2]
coef3.2_list <- apply(dat, MARGIN =2, coef_glm3.2)

sum3.2 <- sum(abs(coef3.2_list))

#3ic
coef_glm3.3 <- function(x) coef(glm(dis ~ x + ic1 + ic2 + ic3,
                                    data = data_full, family = binomial))[2]
coef3.3_list <- apply(dat, MARGIN =2, coef_glm3.3)

sum3.3 <- sum(abs(coef3.3_list))


#3ic
coef_glm3.4 <- function(x) coef(glm(dis ~ x + ic1 + ic2 + ic3 + ic4,
                                    data = data_full, family = binomial))[2]
coef3.4_list <- apply(dat, MARGIN =2, coef_glm3.4)

sum3.4 <- sum(abs(coef3.4_list))


#mini table- we should see the SNP coef 
#sums descend with the addition of each pc 

n_sources <- c("1 source", "2 sources", "3 sources", "4 sources")
PC_sums1 <- c(sum1.1,sum1.2,sum1.3,sum1.4)
FA_sums1 <- c(sum2.1, sum2.2, sum2.3, sum2.4)
IC_sums1 <- c(sum3.1, sum3.2, sum3.3, sum3.4)

display_results <- data.frame(n_population_structure_terms = n_sources, 
                             PC_SNP_coef_total = PC_sums1,
                             FA_SNP_coef_total = FA_sums1,
                             IC_SNP_coef_total = IC_sums1)
display_results 
#ICA shows better sums on average, 
#but none of the sums decreases monotonically
#I expect we have not yet found the best number of sources
#for population structure. 


###A quick visual comparison of coefficient sums

bar_stac <- data.frame(sources=rep(c("PCA","FA","ICA"), each = 4),
                       n_source=rep(n_sources, times = 3),
                       sum = c(PC_sums1 , FA_sums1, IC_sums1))

ggplot(bar_stac, aes(col=sources, y=sum, x=n_source)) + 
  geom_point(shape = 18, size = 5) +
  geom_line(aes(group = sources))+
  #geom_smooth(aes(group = sources), method = "loess", se = FALSE)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))


save.image("hapmap3_full_code.R")

#From the graphs above, it looks as though ICA is consistently better
#at reducing spurious SNP correlations. 
#It is also clear that there is a point of diminishing returns, where
#adding more terms to the glm (in the form of more sources) 
#initially reduces noisy SNPs, but later, amplifies them. 
#the turning point appears to be around roughly the addition of the 7th SNP. 
#much as the tuning parameter in lasso or ridge regression needs to be 
#thoughtuflly selected, so to does the number of sources 

##Lambda1000 


