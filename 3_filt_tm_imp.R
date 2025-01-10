# SDS 20241231

# Filter TOPMed imputation to only the 967 NHW with metabolomics, as described
# by Erika's thesis.

# Setup ------------------------------------------------------------------------

setwd(paste0(Sys.getenv("RKJCOLLAB"), "/Collabs/metaboxcan"))

library(tidyverse)

# Load PLINK file after imputation of HumanOmniExpress array
genetic_nhw <- read_delim(
  "/Users/slacksa/temp_data/metaboxcan/imputed_clean/chr_all_concat.psam")

# Load metabolomics file
metab_mix <- read_csv(
  "raw/metab/COPDGene_P2_LT20miss_knnImp_NoOut_metabs_20211021.csv")

# Load exome chip PLINK file
genetic_ex_nhw <- read_delim(
  "raw/genetics/ExomeQCed_3.25.15_V1.4/NHWexomeQCed_3.25.15_V1.4.fam",
  col_names = F)

# Compare IDs ------------------------------------------------------------------

# Split FID and IID
genetic_nhw_mod <- genetic_nhw %>%
  dplyr::mutate(
    FID = gsub("_[[:alnum:]]+", "", `#IID`),
    IID = gsub("[[:alnum:]]+_", "", `#IID`))
identical(genetic_nhw_mod$FID, genetic_nhw_mod$IID)

# Check overlap
length(intersect(genetic_nhw_mod$IID, metab_mix$sid)) 
  # 957, 10 lower than Erika
  # But my starting N with metab was also lower - 1130 compared to her 1136
  # Our genetic Ns were different too - 6670 compared to her 6760

# See if exome array chip has more overlap - no!
length(intersect(genetic_ex_nhw$X2, metab_mix$sid))  # 940, even lower

# Make ID list for filtering ---------------------------------------------------

id_list <- intersect(genetic_nhw_mod$IID, metab_mix$sid)

genetic_nhw_filt <- genetic_nhw_mod %>%
  dplyr::filter(IID %in% id_list) %>%
  dplyr::select(`#IID`)

write_tsv(
  genetic_nhw_filt,
  "data/genetics/id_list_genetics_and_metab.txt",
  col_names = F)

# Use PLINK to filter ----------------------------------------------------------

# Filter
system(paste0(
  "plink2 --pfile /Users/slacksa/temp_data/metaboxcan/imputed_clean/chr_all_concat ",
  "--keep data/genetics/id_list_genetics_and_metab.txt ",
  "--make-pgen ",
  "--out /Users/slacksa/temp_data/metaboxcan/imputed_clean_filt/chr_all_concat_filt"))
