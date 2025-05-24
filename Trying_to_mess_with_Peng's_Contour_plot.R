#Peng said the contour plot of b0 vs b1 is always a series
#of concentric ellipses 
#with the centre being the values with the minimal RSS 
#I wonder if I can break this rule using 
#b0=0 or b1=0 
#I suppose to do this I will just stretch out the axes to include 
#0's and negatives 
#i might also use a heat map 

b0 <- c(-3:3)
b1 <- c(-4:4)
x<- seq(5,6,length=20) + rnorm(20,sd=0.001)
y<- 2.1*x - 0.6 - rnorm(20,sd=0.0006)
yvx <- lm(y~x)
i <- 1:length(b0)
j <- 1:length(b1)
k <- 1:length(x)
error_square <- (y[k] - b0[i] + b1[j]*x[k])^2
error1_square <- (y[k] - b0[1] + b1[1]*x[k])^2
error1_square_sum <- sum(error1_square)

RSS <- c(1,9,2) + c(3,2,1) #to add vectors DO NOT use sum. 
#sum() will add the elements of all vectors. 
z<- outer(x,b0,b1, FUN = lbf(x,b0,b1)) #it's struggling with this 
source("Commands.R") #see this to see the lbf design

error_square(x,y,b0,b1)
#lbf needs debugging. it can't deal with b0, b1, 
#and x being vectors of different lengths!

h<- c(5,6,7,8)
g<- c(0.1,0.01,1)

h*g #returned 0.5, 0.06, 7.00, 0.8 

g*h #returned the same # it's just multiplying elts 1-1

#R doesn't treat vectors and matrices QUITE the same
#matrix multiplication function is %*%, not * 
#also, it turns out I don't need to convert these
#to matrices first. 

g%*%h #non-conformable (both row vectors or both column vectors)
g%*%t(h) #this produced the desired result! 

length(b0)
length(b1)
length(x)
length(y)
y

b1%*%t(x)
dim(b1%*%t(x)) #9 20 

?lm
z<- lm(y~x)
summary(z)

