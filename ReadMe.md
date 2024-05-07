# BZ562: Computational Approaches in Molecular Ecology Final Project
## Adaptive landscape genetics on the Boreal toad 

#### This data is obtained from [Trumbo et al. 2023](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.17175). This is a conservation landscape genomics study on the Rocky Mountain ecoregion (SRM) of Boreal Toads (_Anaxyrus boreas_). Populations of this species are threatened by chytrid fungus _Batrachochytrium dendrobatidis_ (Bd). The study used ddRADseq to sequence individuals to look at standard metrics of genetic diversity, environmental factors mediating geneflow, and loci associated with environment. 
#### Recently _Bd_ load data was obtained from these individuals. My initial goal for this class project was to conduct a GWAS and estimate heritability for the for _Bd_ load. Given some constrains, I haven't been able to do this, so instead I conducted GEA with the _Bd_ prevalence data and time of first _Bd_ detection in the paper. Trumbo et al. 2023 conducted an RDA for these _Bd_ variables. I attempted to replicate this RDA and in addition do an LFMM to see if there were overlapp in candidate adaptive loci between the two analyses (spoiler: all the RDA loci included the 14 LFMM loci). This is still beneficial for me since I will very likely be doing GEA analyses for my PhD research.

## My workflow
### Downsampling
#### For this study, sequences were filtered and aligned to a _A. boreas_ reference genome which was developed as part of this study. I started with the resulting .vcf file from this alignment (see Trumb et al. 2023 for details0. 
#### GEA analyses require no missing data genotypes at every locus.
#### Therefore, my first step to downsample as much as I could prior to imputing data.
#### To downsample, I used the [genoscapeRtools](https://github.com/eriqande/genoscapeRtools) package.
#### This included using vcftools to first create 012 files which are required for downsampling. 
#### These files are in the VCF_out folder.
#### I ended up downsampling such that I had no individuals with >13% missing data.I retained 212 individuals (80% of current samples).
#### The output includes a set of subsamples 012 files, which I imported back onto Alpine with Globus. These files are in the cleaned_indv212_pos3000 folder.
#### I then used vcftools --keep to subsample my original populations.snp.vcf file. 
#### The downsampled .vcf file is in the Downsampled folder.
### GEA in R
#### I conducted two types of GEA analysis in R studio; LFMM and RDA.
#### I followed the [script](https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html#WE_11) developed by Brenna Forester to conducuct the GEA's.
#### The individual scripts for the LFMM and RDA that I performed are in the Scripts folder.
#### Prior to conducting the analyses, I cleaned the _Bd_ data and genetic data such that only samples in both data sets were included. 
#### I used a simple imputation method in R that used the most common genotype at a locus to impute the missing genotypes.
### Results
#### The LFMM for first year since first _Bd_ detection yielded no candidate loci, while the LFMM for _Bd_ prevalence yielded 14 SNP's
#### The RDA yielded 549 SNP's, which included all SNP's from the LFMM.

###File structure
####ANBO_data: 
