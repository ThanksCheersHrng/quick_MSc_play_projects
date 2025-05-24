#exploring relationship between crime and 
#lower socioeconomic status percentage 
#in the BostonHousing dataset 


library(mlbench)
data(BostonHousing, package = "mlbench")
attach(BostonHousing)

par(mar = c(2,2,2,2)) 
#pre-fix margins so plot will work.  

crim_v_lstat <- plot(BostonHousing$lstat, BostonHousing$crim, col = 4)
#such a heavy cluster around the low end for crime 
#that it's hard to see much of a relationship
log_crim_v_lstat <- plot(BostonHousing$lstat, log(BostonHousing$crim), col = 4)
#this clarifies the relationship somewhat
#but seems to lack homoscedasticity 
log_crim_v_log_lstat <- plot(log(BostonHousing$lstat), log(BostonHousing$crim), col = 4)
#this model looks most suited to linearity. 


#checking homoscedasticity on log_crim only:
log_crim_model <- lm(log(crim)~lstat, data = BostonHousing)
plot(predict(log_crim_model), rstudent(log_crim_model))
abline(0,0, col = 5)
#actually that's pretty 
#decent for heterocedasticity 

#check normality of variances: 
qqnorm(rstudent(log_crim_model))
qqline(rstudent(log_crim_model), col = 6)
#here we see some 
#unexpected concavities 
#but overall decent quantiles to match 

#let's compare those with the 
#model where both variables are log'd

log_both_model <- lm(log(crim)~log(lstat), data = BostonHousing)
plot(predict(log_crim_model), rstudent(log_crim_model))
abline(0,0, col = 7)
#there's a slight 'downhilling' 
#of the residuals
#i.e. the residuals all 
#suddenly seem negative 
#for larger predicted values of crime 

qqnorm(rstudent(log_both_model))
qqline(rstudent(log_both_model), col = 8)
#still some concavity. 
#the quantiles plot actually 
#looks even less linear here.

#from these analyses of our assumptions,
#it would appear the variables 
#most suited to a linear model
#are the crim_log_vs_lstat. 

##what could account for that concavity shape? 
#it looks like a cube root graph. 

lstat_cr <- (BostonHousing$lstat)^(1/3)

log_crim_cr_lstat_model <- lm((log(crim))~lstat_cr, data = BostonHousing)
plot(predict(log_crim_cr_lstat_model), rstudent(log_crim_cr_lstat_model))
abline(0,0, col = 9)
#generally seems to reduce the variance
#but there's still a 'downhilling' effect. 

qqnorm(rstudent(log_crim_cr_lstat_model))
qqline(rstudent(log_crim_cr_lstat_model), col = 8)
#this nearly perfects the middle 
#but wrecks the extremes. 

#let's compare the log_crim options 
#as lines- one with cr_lstat
#the other just as _lstat. 

crim_v_lstat <- plot(BostonHousing$lstat, log(BostonHousing$crim), col = 4)
line_lc_l <- abline(log_crim_model$coefficients, col = 10)
line_lc_crl <- abline(log_crim_cr_lstat_model$coefficients, col = 11)

#this clearly illustrates that the residuals 
#are significantly better for the log_crim model
#than for log_crim_cr_lstat 

#perhaps the curvature is due to some sinusoidal function? 

#but here we may lose 
#interpretability in the name 
#of predictability 
#and our log_crim model is 
#already pretty good as a predictor. 
#to confirm this, let's observe the summary 

summary(log_crim_model)
#the high degree of variability around the linear model
#is accounted for by a multiple r-squared of just 0.39

#the gradient coefficient is roughly 0.2, suggesting
#that a 1% increase in population considered
#'lower socioeconomic status' corresponds to an,
#on average, 0.2 % increase in per capita crime. 

#while the effect of lower status population on crime
#is small, it is likely to be a true correlation, 
#with a p-value less than 2e^-16, and a standard error 
#just roughly 5% of the gradient value (0.01/0.18). 

crim_v_lstat <- plot(BostonHousing$lstat, log(BostonHousing$crim), col = 4)
line_lc_l <- abline(log_crim_model$coefficients, col = 10)

#in extremely short, there is a 
#weak positive correlation
#between crime rates and percentage 
#of the population classed as 
#'lower socioeconomic status.'

