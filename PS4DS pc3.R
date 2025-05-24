#section titles and 
#answers to comprehension questions 
#will be listed in comments

###2 Writing Functions
Myfunction <- function()
  {cat("This is my first function in R!\n") 
  #\n prints a new line
   cat("No arguments are required to run this function")}
Myfunction()

mean.and.sd1 <- function(x) 
{aver <- mean(x) 
st.de <- sqrt(var(x))
c(aver,st.de)}

mean.and.sd2 <- function(x)
{aver <- mean(x)
 st.de <- sqrt(var(x))
 c(mean = aver, st.de=st.de)}

y1 <- mean.and.sd1(rnorm(200))
y2 <- mean.and.sd2(rnorm(200)) 
  
mean.and.sd3 <- function(x)
  {aver <- mean(x) 
  st.de <- sqrt(var(x)) 
  list(mean=aver, st.de=st.de)}
#built-in R functions used: 
  #mean
  #sqrt
  #var
  #function()
  #list
#five built-in R functions used? 
  #more if you count = and <- 
  #then six or seven. 

###2.1 Logical Operators 

x<- seq(1:10)
x < 5
x > 1 & x <=8 #returns per elt of vector
x > 1 && x<=8 #returned FALSE with warning 
#seems unsure whether to use length of vector 
#or first element of vector
#chosen to use just the first element of the vector
x > 3 || x< 10 #similarly dislikes this one. 
#the shorthand in this error is length(x) > 1 
# but logical(1)
x > 3 | x < 10


### 3 Control statements 

library(MASS)
data("anorexia")
attach(anorexia) 
summary(anorexia) #summarises dataset 
summary("anorexia") #just summarises character


if (Postwt[1]>Prewt[1])
  { print("Patient 1 gained weight after treatment") } else
  {print("Patient 1 did not gain weight after treatment")}

#Patient 1 did not gain weight after treatment. 

if (Postwt[15]>Prewt[15])
  + { print("Patient 15 gained weight after treatment") } else
    + {print("Patient 15 did not gain weight after treatment")}

#Patient 15 gained weight after treatment. 


######ASIDE 
#experimental more complex printing. 
#not yet successful.

index_print <- function(i)
{print("Patient" [i])}

treatment_report <- function(i)
{if (Postwt[i]>Prewt[i])
{print("Patient" i "gained weight after treatment!"
         else print("Patient" i "did not gain weight
                    after treatment."))}}



