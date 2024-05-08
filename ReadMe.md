# BZ562: Computational Approaches in Molecular Ecology Final Project
## Adaptive landscape genetics on the Boreal toad 

##### This data is obtained from [Trumbo et al. 2023](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.17175). This is a conservation landscape genomics study on the Rocky Mountain ecoregion (SRM) of Boreal Toads (_Anaxyrus boreas_). Populations of this species are threatened by chytrid fungus _Batrachochytrium dendrobatidis_ (Bd). The study used ddRADseq to sequence individuals to look at standard metrics of genetic diversity, environmental factors mediating geneflow, and loci associated with environment. 
##### Recently, _Bd_ load data was obtained from these individuals. My initial goal for this class project was to conduct a GWAS and estimate heritability for the for _Bd_ load. Given some constrains, I haven't been able to do this, so instead I conducted GEA with the _Bd_ prevalence data and time of first _Bd_ detection in the paper. Trumbo et al. 2023 conducted an RDA for these _Bd_ variables. I attempted to replicate this RDA and in addition do an LFMM to see if there were overlapp in candidate adaptive loci between the two analyses (spoiler: all the RDA loci included the 14 LFMM loci). This is still beneficial for me since I will very likely be doing GEA analyses for my PhD research.

## My workflow
##### For this study, sequences were filtered and aligned to a _A. boreas_ reference genome which was developed as part of this study. I started with the resulting .vcf file from this alignment (see Trumb et al. 2023 for details0. 
### Downsampling
##### GEA analyses require no missing data genotypes at every locus. Therefore, my first step to downsample as much as I could prior to imputing data.
* I used the [genoscapeRtools](https://github.com/eriqande/genoscapeRtools) package to downsample.
* First I used vcftools to create 012 files which are required for downsampling. 
* These files are in the VCF_out folder.
* I ended up downsampling such that I had no individuals with >13% missing data.
* I retained 212 individuals (80% of current samples).
* The output includes a set of subsamples 012 files, which I imported back onto Alpine with Globus. These files are in the cleaned_indv212_pos3000 folder.
* I then used vcftools --keep to subsample my original populations.snp.vcf file. 
* The downsampled .vcf file is in the Downsampled folder.

### GEA in R
##### I conducted two types of GEA analysis in R studio; LFMM and RDA.
##### I followed the [script](https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html#WE_11) developed by Brenna Forester to conducuct the GEA's.
##### The individual scripts for the LFMM and RDA that I performed are in the Scripts folder.
##### Prior to conducting the analyses, I cleaned the _Bd_ data and genetic data such that only samples in both data sets were included. 
##### I used a simple imputation method in R that used the most common genotype at a locus to impute the missing genotypes.
### Results
##### The LFMM for first year since first _Bd_ detection yielded no candidate loci, while the LFMM for _Bd_ prevalence yielded 14 SNP's
##### The RDA yielded 549 SNP's, which included all SNP's from the LFMM.

### File structure
#### [ANBO_data](https://github.com/shanellewiks/ANBO_congen/tree/main/ANBO_data):
* Bd_BodyCondition_clean.csv: Csv file with boreal toad individuals and their correponding Bd and condy condiion data
* Bd_BodyCondition_clean.xml: Excel file of the same data 
* mec17175-sup-0001-appendixs1.docx: Supplementary materials from Trumbo et al. 2023

#### [Scripts](https://github.com/shanellewiks/ANBO_congen/tree/main/Scripts):R scripts
* Downsampling_genoscapeRtools.R: genoscapeRtools script for downsampling
* LFMM.R: Script for running LFMM
* RDA.R: Script for running RDA
* Datacheck.R: Please ignore. Personal script for checking data.


### Challenges
##### Although I decided not to do a GWAS due to the lack of sufficient sample size, this was also in part being unable to impute my files and obtain output files that would let me run a GWAS.
##### My first attempt at imputing was using the [grur package](https://thierrygosselin.github.io/grur/articles/rad_genomics_computer_setup.html), which was developed to impute Radseq data. Unfortunately, this package and some of its dependencies are no longer supported in R. Even compiling packages locally proved impossible for me. 
##### I then tried to impute using [Beagle 4.1](https://faculty.washington.edu/browning/beagle/b4_1.html). Once again, this proved challenging since the loci in the vcf file were not sorted in order. After some attempts at sorting this (including using BCFtools sort and some awk use), I reached another dead end. 
##### My final attempt at imputing was using [LinkImputer](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3873-5). This seemed to work at the start, however, this gave me several errors and also proved to be a dead end. 
##### The simple imputation in R for GEA using the most common genotypes was helpful for those analyses. However, I couldn't find a way to export these filed into the input files required to run GWAS (.vsf's or .ped files from Plink).
##### Although I am disappointed that I could not run a GWAS for this project, I now have some experience doing GEA, which is an analysis I will actually be doing for my dissertation, so I have some practice now!
