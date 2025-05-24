#PS4DS- PC Session 4 

#Introduction

sample(1:0, 10, replace = TRUE) #if replace set to false
#which is the default
#get error message 
#because str only returns integers 


#Ex - flip a coin 100 times 

heads_sample <- sample(1:0, 100, replace = TRUE) #h= 1, t=0 
summary(heads_sample) #tries to give quartiles; not helpful 
table(heads_sample) #58 tails, 42 heads
sum(heads_sample) #counts just the ones (heads) so sum is 42 

#Ex - roll a die 10 times 

coin_sample <- sample(1:6, 10, replace = TRUE) 
#I reran this because the first time I got no 6's 
#that either meant I'd have to do 1:7, or my random
#sample just so happened to give no 6's. 
#when I ran it a second time... 
table(coin_sample) #3 sixes, and interestingly no 3s. 
 

#2 Distributions 
#fill in the gaps below 
#they're random, so will be 
#different every run 

qf(0.05, 6,7,lower.tail =F)
#3.865969

pt(1,3)
#0.8044989

dnorm(1)
#0.2419707

dunif(x=1, min = 0, max = 3)
#0.33333333

qunif(.05, min = 0, max = 3)
#0.15

punif(0.05, 0, 3)
#0.01666667

punif(1.5, 0, 3)
#0.5


#3 Simulation 
# I'll try one with a poisson 
# distribution
# note that the mean should 
# approach the reciprocal
# of the rate 

pois_trials <- rpois(500, 0.25)

csum1 <- cumsum(pois_trials) 

rfreq2 <- csum1/(1:500) #not just over 500 
#because we're calcing for each new trial 

par(mar = c(1,1,1,1), mfrow = c(1,1))

plot(rfreq2, main = "Weak Law Large Poisson", 
     ylab = "Relative Frequency", ylim = c(0,0.3))

abline(0.25, 0)

#see the same thing, 
#approaches 0.25 from below non-monotonically. 



#3.2 Central Limit Theorem 
#1 for loop runs the above code 500 times 

#2 xbar <-c() tells the program to 
#concatenate- i.e. create a vector
#based on the for loop 
#the vector should have 500 elements 

#3 500 (trials) x 10 (samples) = 5000 (numbers) 

#4 run code - what do you observe? 

pdf("hist_xbar.pdf") # prefix

plot(for (j in 1: length(n<-c(5,10,20,50)))
{xbar<-c() 
for (i in 1:m) {xbar[i] <- mean(runif(n[j],a,b))}
hist((xbar),main = "Samp_Dist_of_xbar", xlim = c(0, 500), ylim = c(0, 0.5))})

def.off() #postfix

#try the above to get rid of the plot size issue 
# didn't work. needed to fix margins this time. 

par(mar = c(1.2,1,1.2,1))
par(mfrow=c(2,2))
xbar <- c() 
m = 500 ; a = 0 ; b = 1
for (j in 1: length(n<-c(5,10,20,50)))
{xbar<-c() 
 for (i in 1:m) {xbar[i] <- mean(runif(n[j],a,b))}
  hist((xbar))}

#figure margins too large warning came again. 
#fixed the margins with par(mar = c())

#what do you see as n increases? 
# the distribution shifts from sortof uniform to sortof normal.
# strangely, there doesn't seem to be a relation between 
# width of the bars and size of n.  

#note to self the par(mfrow...) makes it go row by row, 
#like English reading 

#####Exercises 
# 1 Generate random sample n=100 from U[-5,5]
# Find min and max from the data. 

Uni_Dist <- runif(100, min=-5, max = 5)
min(Uni_Dist) #-4.944168
max(Uni_Dist) # 4.989253 

#2. For the above sample, find how many numbers
#are negative using a for loop and if statement 

#needed a baby version of the count first: 
cater = c(1, -2, 3, -4, 5)
total <- 0
for (i in 1:5) {if (cater[i] < 0) {total = total + 1}} ; total

#Uni_Dist version: 

total <- 0
for (i in 1:length(Uni_Dist)) 
  {if (Uni_Dist[i] < 0) {total = total + 1}} ; total 

#yay! Returns 54. 

#3. Amend the above to work for any U[a,b], size n. 

# I think I already have. 
Uni_Dist <- runif(190, min=-1, max = 19)
total <- 0
for (i in 1:length(Uni_Dist)) 
{if (Uni_Dist[i] < 0) {total = total + 1}} ; total 
#yup, this time it returned 9, which seems sensible. 

# or do they want us to separately program n, a, and b? 

n11 = 200
a11 = -1
b11 = 20
neg_count_func <- function(n1, a1, b1)
{
Uni_Dist = runif(n1, min=a1, max =b1)
  total <- 0
for (i in 1:length(Uni_Dist)) 
{if (Uni_Dist[i] < 0) {total = total + 1}} ; total } 

neg_count_func(n11, a11, b11)
# returned 12 . also sensible. Hooray! 