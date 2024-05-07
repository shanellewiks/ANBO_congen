### RDA ###
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

#RDA is a particulary useful tool because it handles both response and predictor variables being multivariate
#The basic premise is that ir detects loci that covary with env variables

####1) Import genomic data ####
#This data is the downsampled dataset from the genoscapeRtools workflow

loci<- read.vcfR(here("Downsampled", "cleaned_212I_3KP.vcf.recode.vcf")) #import vcf 

anbo.genind <- vcfR2genind(loci) #convert vcf to genind
anbo.genind

anbo.gl<- vcfR2genlight(loci) #convert to genlite
anbo_genpop

#### 2) Imputing data ####

#Again, using the most common genotype at a locus to impute
anbo_genpop<- genind2df(anbo.genind) #convert to df 

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
bd_clean<- bd_bodycon %>% 
  mutate(Individual = as.character(Individual)) %>% #making sure sample id's are all characters and not factors 
  filter(!is.na(first_bd)) #filter missing env data
view(bd_clean)

####4) Subset data ####

#what samples are in the bd data that aren't in the rad data 
radsind<- rownames(gen.imp) #rad seq samplenames
bdind<- bd_clean$Individual #body condition sample individuals

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

# First filter samples to only include the ones in common
bd_clean <- bd_clean %>% 
  filter(Individual %in% common) 

view(bd_clean)

# Filter samples from the imputed rad data to only include the ones in common
genind.imp.clean<- genind.imp[common]
genind.imp.clean

#check if they are identical 
identical(bd_rownames, radrownames)
#>TRUE

#also forgot that I have to filter out Wyoming since there's no prevalence data for it
#I'll clean the data up a bit while I'm at it to give shorter names
bd_clean<-bd_clean %>% 
  filter(prevalence != "Unknown") %>% 
  mutate(prevalence=as.numeric(prevalence)) %>% 
  rename(fbd =first_bd,
         prev = prevalence)
view(bd_clean)
keep <- bd_clean$Individual

#filtering the genind data

gen.imp.clean<-genind.imp.clean[keep]
gen.imp.clean

#check if they are identical 
radsind<- rownames(gen.imp.clean) #rad seq samplenames
bdind<- bd_clean$Individual #body condition sample individuals

identical(bd_rownames, radrownames)
#>TRUE


####5) RDA ####

#First lets select only the environmental data 
pred<-bd_clean[,6:7]
colnames(pred)<-c("fbd", "prev")
pred

#run the RDA

anbo.rda <- rda(gen.imp.clean ~ ., data=pred, scale=T)
anbo.rda

#proportion of variance explained by env predictors is represented by the proportion column for the 
#constrained row. This is similar to a R squared value.
#Here it is, 1.399e-02 

#Adjusting the 'R squared' value
RsquareAdj(anbo.rda)
#0.0009286616
#therefore constrained ordination explains about 0.092% of the variation
#Apparently we shoudn't be alarmed since technically most SNP's are nuetral and shouldn't explain much variation

#The eigen values for the constrained axes reflect the variance explained by each canonical axis 
summary(anbo.rda)$concont

#visualize using a screeplot
screeplot(anbo.rda)

#running a formal test of statistical significance of each constrained axis
#anova.cca(anbo.rda, by="axis") #DONT RUN THIS IT TAKES AGES
#H0: no linear relationship between SNP data and env predictors 

#You can check Variance Inflation Factors for pedictor variables in the model to help find
#autocorrelated variables 
vif.cca(anbo.rda)
#Values below 10 or 5 are indicative of autocorrelation not being a problem

#RDA BIPLOT!


plot(anbo.rda, scaling=2)

#SNPS are in red in the center
#Individuals are the black circles 
#environmental predictors are in blue arrows

#Identify RDA candidate loci
#This is done using the loadings
load.rda<-summary(anbo.rda)$species[,1:3]

#visualize loadings with histograms 
#SNPs associated with env predictors are at the tail of the histograms
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")

#we can identify those SNPs with this script
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}


cand1 <- outliers(load.rda[,1], 3) 
cand2 <- outliers(load.rda[,2], 3)

anbo.rda.cand<-c(names(cand1), names(cand2))
length((anbo.rda.cand[duplicated(anbo.rda.cand)]))
#1 duplicate detection

anbo.rda.cand <- anbo.rda.cand[!duplicated(anbo.rda.cand)]
length(anbo.rda.cand)
#549 unique snp's??

#A nicer RDA biplot

# Set up the color scheme for plotting:
bgcol  <- ifelse(colnames(gen.imp.clean) %in% anbo.rda.cand, 'gray32', '#00000000')
snpcol <- ifelse(colnames(gen.imp.clean) %in% anbo.rda.cand, 'red', '#00000000')

## a.es 1 & 2 - zooming in to just the SNPs here...
plot(anbo.rda, type="n", xlim=c(-1,1), ylim=c(-1,1), main="ANBO RDA, axes 1 and 2")
points(anbo.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6')
points(anbo.rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol)
text(anbo.rda, scaling=3, display="bp", col="#0868ac", cex=1)

## a.es 2 & 3
plot(anbo.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="ANBO RDA, axes 2 and 3")
points(anbo.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(2,3))
points(anbo.rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3, choices=c(2,3))
text(anbo.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))


#We can compare this with the LFMM candidates
snps_both<-intersect(vals, anbo.rda.cand)
#14 snp's below
snps
snps_both
#basically all snp's identified in the prevalence LFMM were identified in the RDA