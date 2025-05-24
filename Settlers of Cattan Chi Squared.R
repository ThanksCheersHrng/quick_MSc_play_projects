# My husband thinks our Settlers of Cattan dice are unfair. 
# He got the following frequencies, for the numbers, 
#in order: 20, 18, 9, 20, 17, 15 

Obs = c(20, 18, 9, 20, 17, 15)

#H_0 : die equally likely to land on any number 
#H_1 : die NOT equally likely to land on any number 

# Expected frequencies: 16.67 

Exp = rep.int(16.67, 6)

# Test statistic; 

#X^2 = sum[(obs - expect)^2/expect]

SqDiff = (Obs - Exp)^2

Test_stat = sum(SqDiff/Exp)

#Finding Critical Value 

alpha = 0.01 

dfree = 5 #6 sided die

Crit_Val = qchisq(alpha,dfree)

#Final Test 

#There's evidence to reject the null if  
#the following is true: 

as.logical(Test_stat > Crit_Val)

#We have evidence to reject the null hypothesis. 
#We may believe the die is unfair. 
#We'll have to get new dice to play. 

#Upload to git repository
usegit()

