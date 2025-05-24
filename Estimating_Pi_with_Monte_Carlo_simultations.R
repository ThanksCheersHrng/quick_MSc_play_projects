# Monte Carlo sampling pi 

set.seed(82)
n = 1000
x<-runif(n,0,1)
y<-runif(n,0,1)
plot(x,y,xlab="x values",ylab="y values",main="random points in unit square",xaxt="n",
     yaxt="n",pch=20,cex=.7) 

inside<-x*x+y*y<1
plot(x,y,xlab="x values",ylab="y values",main="in or out the circle?",xaxt="n",
     yaxt="n",pch=20,cex=.7,col=inside+1)

red <- sum(inside)
piiii <- red*4/n
piiii

## Approximate expectation of standard normal dist

x <- rnorm(1e6,0,1) #1e9 too much for computer 
mean(x)
y <- rnorm(10,0,1) #not enough samples
mean(y)
#just do as many samples as possible for computer 

