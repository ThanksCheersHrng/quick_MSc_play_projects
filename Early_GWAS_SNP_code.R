#this function returns SNP coefficients for the model structure with no population structure
coef_glm0 <- function(x) coef(glm(dis ~ x, data = data_full, family = binomial))[2]

coef_glm0(dat$"3") 

#apply for every column (MARGIN = 2 for columns) 
coef0_list <- apply(dat, MARGIN = 2, coef_glm0)

sum0 <- sum(abs(coef0_list))
#sum0 will give sum of coefficients of SNPs when not controlled with any pop struc (i.e. glm0)


#Write the functions... 
#GLMs with PCs as pop struc:
coef_glm1.1 <- function(x) coef(glm(dis ~ x + pc1, 
                                    data = data_full, family = binomial))[2]
coef_glm1.2 <- function(x) coef(glm(dis ~ x + pc1 + pc2, 
                                    data = data_full, family = binomial))[2]
coef_glm1.3 <- function(x) coef(glm(dis ~ x + pc1 + pc2 + pc3,
                                    data = data_full, family = binomial))[2]
coef_glm1.4 <- function(x) coef(glm(dis ~ x + pc1 + pc2 + pc3 + pc4,
                                    data = data_full, family = binomial))[2]
#GLMs with Factors as pop struc: 
coef_glm2.1 <- function(x) coef(glm(dis ~ x + fa1.1,
                                    data = data_full, family = binomial))[2]
coef_glm2.2 <- function(x) coef(glm(dis ~ x + fa2.1 + fa2.2,
                                    data = data_full, family = binomial))[2]
coef_glm2.3 <- function(x) coef(glm(dis ~ x + fa3.1 + fa3.2 + fa3.3,
                                    data = data_full, family = binomial))[2]
coef_glm2.4 <- function(x) coef(glm(dis ~ x + fa4.1 + fa4.2 + fa4.3 + fa4.4,
                                    data = data_full, family = binomial))[2]

#GLMs with ICs as pop struc: 
coef_glm3.1 <- function(x) coef(glm(dis ~ x + ic1,
                                    data = data_full, family = binomial))[2]
coef_glm3.2 <- function(x) coef(glm(dis ~ x + ic1 + ic2,
                                    data = data_full, family = binomial))[2]
coef_glm3.3 <- function(x) coef(glm(dis ~ x + ic1 + ic2 + ic3,
                                    data = data_full, family = binomial))[2]
coef_glm3.4 <- function(x) coef(glm(dis ~ x + ic1 + ic2 + ic3 + ic4,
                                    data = data_full, family = binomial))[2]


#Use apply() on each function and calculate sum of SNP coefficients 

#PCA: 
coef1.1_list <- apply(dat, MARGIN = 2, coef_glm1.1)
sum1.1 <- sum(abs(coef1.1_list))

coef1.2_list <- apply(dat, MARGIN = 2, coef_glm1.2)
sum1.2 <- sum(abs(coef1.2_list))

coef1.3_list <- apply(dat, MARGIN =2, coef_glm1.3)
sum1.3 <- sum(abs(coef1.3_list))

coef1.4_list <- apply(dat, MARGIN =2, coef_glm1.4)
sum1.4 <- sum(abs(coef1.4_list))

#FA: 
coef2.1_list <- apply(dat, MARGIN =2, coef_glm2.1)
sum2.1 <- sum(abs(coef2.1_list))

coef2.2_list <- apply(dat, MARGIN =2, coef_glm2.2)
sum2.2 <- sum(abs(coef2.2_list))

coef2.3_list <- apply(dat, MARGIN =2, coef_glm2.3)
sum2.3 <- sum(abs(coef2.3_list))

coef2.4_list <- apply(dat, MARGIN =2, coef_glm2.4)
sum2.4 <- sum(abs(coef2.4_list))

#ICA: 
coef3.1_list <- apply(dat, MARGIN =2, coef_glm3.1)
sum3.1 <- sum(abs(coef3.1_list))

coef3.2_list <- apply(dat, MARGIN =2, coef_glm3.2)
sum3.2 <- sum(abs(coef3.2_list))

coef3.3_list <- apply(dat, MARGIN =2, coef_glm3.3)
sum3.3 <- sum(abs(coef3.3_list))

coef3.4_list <- apply(dat, MARGIN =2, coef_glm3.4)
sum3.4 <- sum(abs(coef3.4_list))
