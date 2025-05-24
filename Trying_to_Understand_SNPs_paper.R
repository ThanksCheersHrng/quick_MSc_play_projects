# Trying to understand SNPs 
# @Principal Components Analysis corrects for stratificiation 
# in genome-wide association studies 
# by Alkes L Price Et Al 

#example data given in figure 1 
Axis_of_Variation <- matrix(c(0.7,0.4,-0.1,-0.4,-0.5)) #generated from PCA 
Candidate_SNP <- matrix(c(2,2,1,1,0))
Candidate_Phenotype <- matrix(c(1,1,0,0,0))

score1 = sum(Axis_of_Variation*Candidate_SNP) #1.7
Axis_of_Variation*Candidate_SNP #doesn't give the list they're referring to. 0.o
Axis_of_Variation*Candidate_Phenotype #neither does this. 0.o 

# presume Axis_of_Variation*Candidate_SNP is the projected C_SNP onto A_o_V
# How do mathematically? Ans should be 1.0, 1.4, 1.1, 1.6, 0.8 
# a_vector = a_proj * b_vector #orthogonal projection 
# a_proj = magnitude(a_vector)*cos(angle between a and b) 
# or a_proj = a_vect dot b_vect
# so then the projection itself is the scalar a_proj times vector b_vector 

score1*Axis_of_Variation #that is not the answer I was expecting. 

#Where are the numbers in figure 1 coming from?? 

score1*Candidate_SNP #nope. Thank goodness bc that would make zero sense. 

#Let's play with the chisquared business for a bit. 
Mystery_Vector_From_Can_SNP <- matrix(c(1.0,1.4,1.1,1.6,0.8)) 
Mystery_Vector_From_Can_Phen <- matrix(c(0.3,0.6,0.1,0.4,0.5)) 

plot(Mystery_Vector_From_Can_Phen~Mystery_Vector_From_Can_SNP, pch = 19, col = "olivedrab")

MAT <- rbind(tabulate(Mystery_Vector_From_Can_Phen),tabulate(Mystery_Vector_From_Can_SNP))

chisq.test(MAT) # gives chisq = 4, df = 1, p-value = 0.0455 
#this is not the example of 0.07 from figure 1- I don't know what they're doing. 
 


# Questions/Observations... 

# I can't get the example from figure 1 
# to work- what are they leaving out? 

# Can we clarify that I have the right 
# conception of what they're doing with these 
# SNPS to generate data? 

# "10 data sets of 100,000 random SNPS"
# how long is the chain where they're SNPing? 
# also 100,000? 
# should I be picturing each of the 10^6 samples
# as a chain of 100,000 single nucleotide pairs, 
# and each sample has a *different* polymorphism? 

# Bunch of references to extra names 
# and methods- I assume we don't need to worry about this
# as it would be excessive specificity 
# and we'll be using our own methods? 
# e.g. "Armitage" trend chisq stat
# e.g. Balding-Nichols model 

# One advantage of their methodology in PCA 
# is that the top ten will include the top 1, 
# unlike in FA, where the top two factors
# doesn't contain the top one factor -- they are different models. 
# how does it work for ICA? 

# One critique I also have is the presumption that #
# population structure should only have one axis of variation
# surely in practice mixed race people are a counterexample to this? 

# By "leading to a loss in power", I assume they mean *not* that 
# there's somehow a loss of statistical power, 
# but rather that there's somehow a loss in the 
# ability to predict the population structure
# once they control for (remove) the differentiating SNP? 

# " We caution that this does not obviate the need to carefully match cases and
# controls when designing a disease study: in the current example, a
# more closely matched set of 500 cases and 500 controls would have
# achieved superior power to detect true associations."
# Now it sounds like they *are* using power in the statistical sense. 
# I think I need a diagram.

# I notice r > 1 in their settings. 
# In table 1 - "proportion of associations reported as significant" 
# Is that 1 - Power? 
# what are their definitions for 'random,' 'differentiated' and 'casual' SNPs? 
# I can guess but I'd rather be specific. 

# What does F_ST refer to? 

# False Positive -- e.g. 'you are of Ashkenazi descent so you 
# do have Gaucher disease' --- makes sense PCA would just plop 
# the disease in with the pop struc, but how did they 'demonstrate above' 
# that highly differentiated candidate SNPs are particularly 
# likely to produce false positive associations

# big idea for us: what's the point of doing ica when pca looks so good 99.9% 