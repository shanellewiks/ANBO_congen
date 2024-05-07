#### LFMM to identify candidtate loci associated with Bd in Boreal toads ####
#based on Brenna Forester's LandGen workshop for GWE#
#https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html

library(LandGenCourse)
library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)
library(adegenet)
library(vcfR)
library(here)
library(tidyverse)

####1) Import genomic data ####
#This data is the downsampled dataset from the genoscapeRtools workflow

loci<- read.vcfR(here("Downsampled", "cleaned_212I_3KP.vcf.recode.vcf")) #import vcf 

anbo.genind <- vcfR2genind(loci) #convert vcf to genind
anbo.genind

anbo.gl<- vcfR2genlight(loci) #convert to genlite
anbo_genpop

####2) Imputing ####
#data must be imputed to run lfmm
anbo_genpop<- genind2df(anbo.genind) #convert to df of sample and snps

dim(anbo_genpop)
#212 samples 30008 loci


#imputing with the most common genotype
sum(is.na(anbo_genpop))
#254912 NA's in the matrix

#imputing using the most common SNP at each genotype
gen.imp <- apply(anbo_genpop, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) #imputing
sum(is.na(gen.imp)) #no more Na's

#convert back to geinind
genind.imp <- df2genind(gen.imp, sep = ",")
genind.imp


####3) Importing the environmental data ####
#these are the predictor variables in an RDA

bd_bodycon<- read_csv(here("Bd_BodyCondition_clean.csv")) 
view(bd_bodycon)

#clean the data 
bd_bodycon_clean<- bd_bodycon %>% 
  mutate(Individual = as.character(Individual)) %>% #making sure sample id's are all characters and not factors 
  filter(!is.na(first_bd)) #filter missing env data
view(bd_bodycon_clean)

####4) subset the env and gen data ####

#confirm genotypes and envrionmental data are in the correct sequence 
identical(rownames(gen.imp), bd_bodycon_clean[,1])

#they are not so we need to sort this out 

#what samples are in the bd data that aren't in the rad data 
radsind<- rownames(gen.imp) #rad seq samplenames
bdind<- bd_bodycon_clean$Individual #body condition sample individuals

diff_bd<- setdiff(bdind,radsind)
view(diff_bd)
#There are 21 samples in the bd data that aren't in the rad data

#what samples are in the rad data that aren't in the bd data 
diff_rad<- setdiff(radsind,bdind)
view(diff_rad)
#There are 2 samples in the bd data that aren't in the rad data

#what samples are in common
common <- intersect(bdind, radsind)
view(common)
#This leaves us with 157 samples for the GEA

# We need to filter each of these data frames so that they have only the relavent samples that can be includedin the analysis 

# First filter samples to only include the ones in common
bd_bodycon_clean <- bd_bodycon_clean %>% 
  filter(Individual %in% common) 
  
view(bd_bodycon_clean)

# Filter samples from the imputed rad data to only include the ones in common

genind.imp.clean<- genind.imp[common]
#convert this to a matrix
gen.imp.clean<- genind2df(genind.imp.clean)
view(gen.imp.clean)

#lets check if things are in the correct order again
identical(rownames(gen.imp.clean), bd_bodycon_clean[,1])

#Still not in the same order so lets try to reorder the env data
radrownames<- rownames(gen.imp.clean) #rad seq samplenames
bd_rownames<- bd_bodycon_clean$Individual

identical(bd_rownames, radrownames)

#SUCCESS!

# subsetting env variables 
pred<-bd_bodycon_clean %>% 
  select(first_bd,prevalence) %>%  #subsetting 
  rename(fbd = first_bd, 
         prev = prevalence)

#Since we have only 2 predictors, I will run separate RDA's for each of these 

####5) LFMM: Latent Factor Mixed Model ####
# LFMM is a univariate test, so it builds a model for each SNP and each predictor variable
# This means 2* 30,008 = 60,016 tests in total

#LFMM requires an estimated number of populations, K (or run LFMM multiple times with different K values)
#Here, I run kmeans clustering to estimate the number of populations.
#This used Bayesian Information Criterion (BIC) to determine number of clusters
#The lowest BIC value corresponds to the estimated number of populations

cluster <- find.clusters(genind.imp.clean,max.n.clust = 4, n.pc=200)
view(cluster$Kstat)
plot(cluster$Kstat)
#Lowest BIC is for 3 populations. Therefore: 
K<- 3

#now we run the LFMM's

#5a)LFMM for year since first detection of Bd

first_bd.lfmm <- lfmm_ridge(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, K=K)

#Then we identify candidate loci using false discovery rate
#We calculate test statistics for the predictor 
first_bd.pv <- lfmm_test(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, lfmm=first_bd.lfmm, calibrate="gif")

names(first_bd.pv)

# I will look at the Genomic Inflation Factor (GIF). 
#This is an indicator of how well the model accounted for confounding variables. We want a value closer to 1.
# gif>1= too many candidate SNPs; gif<1= too stringent 

first_bd.pv$gif

#> 1.495071
#> 
#We will visualize the influence of GIF to p values. What you need is a flat histogram, with a peak at 0.
hist(first_bd.pv$pvalue[,1], main="Unadjusted p-values")        
hist(first_bd.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")


#this is not what we have so we need to adjust for a more conservative GIF value
#One option is to adjust the K values since this affects GIF. but for now we will adjust the GIF.
# Let's change the GIF and readjust the p-values:
zscore <- first_bd.pv$score[,1]   # zscores for first predictor, we only have one in our case...
(gif <-first_bd.pv$gif[1])       ## d.fault GIF for this predictor
#> 1.495071 original GIF


new.gif1<- 1.2 #randomly chose a more conservative GIF
adj.pv1 <-pchisq(zscore^2/new.gif1, df=1, lower = FALSE)

#plot the p values again
hist(first_bd.pv$pvalue[,1], main="Unadjusted p-values")        
hist(first_bd.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=2.8)")
hist(adj.pv1, main="REadjusted p-values (GIF=2.0)")

#none of these look great so I need to adjust K again
#This time I will run a PCA to determine K 

anbo.gl<- vcfR2genlight(loci)
view(anbo.gl)

anbo.gl.clean<- anbo.gl[anbo.gl@ind.names  %in% common] #subset the genlight object for only samples in included
anbo.gl.clean
anbo.pca <- glPca(anbo.gl.clean)

scatter(anbo.pca) #scatter plot for pca

myCol <- colorplot(anbo.pca$scores, anbo.pca$scores, transp=TRUE, cex=4) #Rough figure made from package 
abline(h=0,v=0, col="grey")
text(anbo.pca$scores[,1], 
     anbo.pca$scores[,2], 
     cex=0.3)

#the number is still a little weird
#I'm going to try this out for K values = 5, 6, 7 and then 4


#Testing K = 5

K <-5
first_bd.lfmm <- lfmm_ridge(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, K=K)
first_bd.pv <- lfmm_test(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, lfmm=first_bd.lfmm, calibrate="gif")

names(first_bd.pv)

# GIF
first_bd.pv$gif

#> 1.230937 #looks better than before
#> 
hist(first_bd.pv$pvalue[,1], main="Unadjusted p-values")        
hist(first_bd.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values") #dont love these histograms


#Testing K = 6

K <-6
first_bd.lfmm <- lfmm_ridge(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, K=K)
first_bd.pv <- lfmm_test(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, lfmm=first_bd.lfmm, calibrate="gif")

names(first_bd.pv)

#GIF
first_bd.pv$gif

#> 1.239063 #not a huge difference
#> 

hist(first_bd.pv$pvalue[,1], main="Unadjusted p-values")        
hist(first_bd.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values") #dont love these histograms

#Testing K = 7

K <-7
first_bd.lfmm <- lfmm_ridge(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, K=K)

first_bd.pv <- lfmm_test(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, lfmm=first_bd.lfmm, calibrate="gif")

names(first_bd.pv)

#GIF
first_bd.pv$gif

#> 1.242184 #this is worse ...
#> 
hist(first_bd.pv$pvalue[,1], main="Unadjusted p-values")        
hist(first_bd.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values") #dont love these histograms

#Testing K = 4 for the sake of it

K <-4
first_bd.lfmm <- lfmm_ridge(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, K=K)

first_bd.pv <- lfmm_test(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, lfmm=first_bd.lfmm, calibrate="gif")

names(first_bd.pv)

# GIF
first_bd.pv$gif

#> 1.383117 #this is also worse ...
#> 

hist(first_bd.pv$pvalue[,1], main="Unadjusted p-values")        
hist(first_bd.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values") #dont love these histograms

## So far K = 5 has been the best K value so we will stick to this ...

K <-5
first_bd.lfmm <- lfmm_ridge(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, K=K)

first_bd.pv <- lfmm_test(Y=genind.imp.clean, X=bd_bodycon_clean$first_bd, lfmm=first_bd.lfmm, calibrate="gif")

names(first_bd.pv)

first_bd.pv$gif

#> 1.230937 
#> 
#We will visualize the influence of GIF to p values. What you need is a flat histogram, with a peak at 0.
hist(first_bd.pv$pvalue[,1], main="Unadjusted p-values")        
hist(first_bd.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values") #dont love these histograms
 
# I will manually adjust this new gif value to see if this is bette

zscore <-first_bd.pv$score[,1]   # zscores for predictor
(gif <- first_bd.pv$gif[1])      

new.gif<- 1.0 #new gif chosen

#manually adjust the p value below
adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)

#plotting all the p values 
hist(first_bd.pv$pvalue[,1], main="Unadjusted p-values")        
hist(first_bd.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=2.8)")
hist(adj.pv1, main="REadjusted p-values (GIF=1.0)")


#I will stick to default GIF

#Now we adjust the p values to q values. 
#q values indicate the significance of each SNP while accounting for the sheer amount of SNPs being tested
first_bd_qv<-qvalue(first_bd.pv$calibrated.pvalue)$qvalues
length(which(first_bd_qv<0.1))

#Oh dear ... 0 SNPS's

(first_bd.FDR.1 <- colnames(genind.imp.clean)[which(first_bd_qv < 0.1)])

#Nothing ....

#Lets try this for the prevalence data instead
# We already know K so this makes things faster


#5b)LFMM for prevalence

#not as intuitive as I thought since there's no prevalence data for Wyoming so we need to exlude it now

bd_no_wy<-bd_bodycon_clean %>% 
  dplyr::filter(prevalence != "Unknown") %>%  #new bd dataframe
  mutate(prevalence = as.numeric(prevalence))
view(bd_no_wy)
keep<- bd_no_wy$Individual #vector w samples to keep 

genind.no.WY<- genind.imp.clean[keep] #filter genind 
genind.no.WY


# Test for K again using a PCA
anbo.gl.clean2<- anbo.gl[anbo.gl@ind.names  %in% keep] #subset genlight
anbo.gl.clean2
anbo.pca <- glPca(anbo.gl.clean2)

scatter(anbo.pca) #scatter plot for pca

myCol <- colorplot(anbo.pca$scores, anbo.pca$scores, transp=TRUE, cex=4) #Rough figure made from package 
abline(h=0,v=0, col="grey")
text(anbo.pca$scores[,1], 
     anbo.pca$scores[,2], 
     cex=0.3)

#still looks like K = 5

K<- 5

prev.lfmm <- lfmm_ridge(Y=genind.no.WY, X=bd_no_wy$prevalence, K=K)

#Calculate test statistics for the predictor 
prev.pv <- lfmm_test(Y=genind.no.WY, X=bd_no_wy$prevalence, lfmm=prev.lfmm, calibrate="gif")

names(prev.pv)

# I will look at the Genomic Inflation Factor (GIF). 
prev.pv$gif

#> 1.130096 looks lower 
#> 
#We will visualize the influence of GIF to p values. What you need is a flat histogram, with a peak at 0.
hist(prev.pv$pvalue[,1], main="Unadjusted p-values")        #looks better than calibrated
hist(prev.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

#transforming p values to q values
prev.qv <- qvalue(prev.pv$pvalue)$qvalues
prev.qv<-as.data.frame(prev.qv) %>% 
  rownames_to_column(var="loci")
view(prev.qv)
length(which(prev.qv<0.1)) #This means False discovery rate of 10%

#YAAAAS! There's 14 candidate loci thank goodness 
prev.qv<- prev.qv %>% 
  filter(V1 <0.1)
vals<- prev.qv$loci

#The loci are: 
vals

