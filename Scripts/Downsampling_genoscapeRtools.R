### Down sampling for imputation ###
### May 1st 2024 ###
# This uses the genoscapeRtools package to find which individuals and positions have the least. 
# Data used here is 265 Boreal toads _Anaxyrus boreas_ from the Rocky Mountain populations (Trumbo et al. 2023)
# You can down sample based on these results
# code outputs are shown as 
#>

### Libraries ###
library(genoscapeRtools)

### Import data ###
# This data is obtained by running the -012 function in vcftools which creates 3 files, .012, .indv, .pos (and log)
# Direct the read_012 function to the -012 prefix
anbo <- read_012(prefix = "/Users/Shanelle/Dropbox/ANBO_congen/VCF_out/out")
#> 
#Read 68818 items
#Read 265 items

### Calculating missing data by individuals ###
indv <- miss_curves_indv(anbo)
indv$plot

### calculating missing data by locus ###
loci <- miss_curves_locus(anbo)
loci$plot

# GWAS is very demanding 
# It wants no missing data and huge sample sizes 
# Our sample sizes are already too low, so I wan't to keep as many individuals as possible
# But also dont want loci with a ton of missing data
# I will impute the rest

# Cleaning based on 212 individuals (80% of current sample size) and 30,000 position should mean we have no individuals with >13% missing data
clean <- miss_curves_indv(anbo, clean_pos = 30000, clean_indv = 212)

#plot the clean data
clean$plot

# I'll look at the resulting files: .012, .pos, .ind, with "cleaned_indv212_pos30000" prefix

dir(pattern = "cleaned_indv212_pos30000")
#>
#[1] "cleaned_indv212_pos30000.012"    "cleaned_indv212_pos30000.012.indv"
#[3] "cleaned_indv212_pos30000.012.pos" 

# We can read these back in to look at the distributions

anbo_clean <- read_012("cleaned_indv212_pos30000")

anbo_clean_miss <- anbo_clean == -1
missing_perc_in_indvs <- rowSums(anbo_clean_miss) / ncol(anbo_clean_miss)
missing_perc_in_loci <- colSums(anbo_clean_miss) / nrow(anbo_clean_miss)
par(mfrow=c(1,2))
hist(missing_perc_in_indvs, main = "Individuals")
hist(missing_perc_in_loci, main = "Positions")

# Missing percent in loci look good, but not so much in individuals. I'm just going to impute the rest.
# Now I will put these back on the server to get a vcf with only the individuals and loci that I want

# I used globus to import these files to Alpine
# I used vcftools to make a new downsamples vcf file. See code below

#vcftools --vcf populations.snps.vcf --out cleaned_212I_3KP.vcf --keep cleaned_indv212_pos30000/cleaned_indv212_pos30000.012.indv --positions cleaned_indv212_pos30000/cleaned_indv212_pos30000.012.pos --recode
#> 
#> 
#VCFtools - 0.1.16
#(C) Adam Auton and Anthony Marcketta 2009

#Parameters as interpreted:
#--vcf populations.snps.vcf
#--keep cleaned_indv212_pos30000/cleaned_indv212_pos30000.012.indv
#--out cleaned_212I_3KP.vcf
#--positions cleaned_indv212_pos30000/cleaned_indv212_pos30000.012.pos
#--recode

#Keeping individuals in 'keep' list
#After filtering, kept 212 out of 265 Individuals
#Outputting VCF file...
#After filtering, kept 30008 out of a possible 34409 Sites
#Run Time = 10.00 seconds




