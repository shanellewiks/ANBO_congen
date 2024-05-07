#### checking available data###
#load packages#
library(here)
library(tidyverse)
library(stringr)

#radseq sample names
radseq_samples<- read_csv(here("Radseq_sampleIDs.csv"))
view(radseq_samples)

#ANBO master spreadsheet
bd_bodycon<- read_csv(here("Bd_BodyCondition_clean.csv"))
view(bd_bodycon)

bd_bodycon_clean<- bd_bodycon %>% 
  dplyr::rename(Sample_ID = "Individual") 
view(bd_bodycon_clean)


check <- bd_bodycon_clean %>%   
  right_join(radseq_samples, by= "Sample_ID")
#231 of 265 radseqed individuals have bd samples
#107 of those have bd load = 0
#22 of those have non zero Bd values
#76 are tadpoles
#20 are metamorphs





