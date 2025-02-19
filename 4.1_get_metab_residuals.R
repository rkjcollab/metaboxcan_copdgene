# SDS 20250122

# Modified from the script MLR.R in PredictDB tool on BioDataCatalyst.

# Running with two different sets of covariates:
  # 1. age and sex (most closely matches paper draft)
  # 2. age, sex, smoking status, pack years, BMI, and FEV1pp (suggested
    # by Katerina)
    # TO NOTE: do not need to adjust for race because already subset to NHW

# Subset to the 957 individuals that calculating accuracy in before running
# model to get residuals here. To do this and keep all covars for option 2,
# reduce to 953 individuals. Will run option 1 with all 957 and with 953 to
# be comparable.

# TO NOTE: select the above options and N individuals in script below.

# R Script used for the Multiple Linear Regression
# Script will loop through all the columns of the transposed gene expression
# which correspond to each  gene and for each gene it runs linear regression
# on the covariates. Then it sets the residuals to the new expression for that
# gene.

# Omics file format: A tab delimited file with .tab file extension containing
# N + 1 rows and G + 1 columns, where N is the number of samples, and G is the
# number of features (genes, methylation sites, chromatin accessibility windows,
# etc.). The first row and column must contain sample IDs and feature IDs
# respectively. Feature values should be normalized across samples and variance
# stabilized.
# Covariates file format: A tab deliminated file with N+1 rows and K+1 columns,
# where N is the number of samples, and K is the desired covariates.
# Output prefix format: File name prefix for output files.

library(tidyverse)
setwd(paste0(Sys.getenv("RKJCOLLAB"), "/Collabs/metaboxcan"))

omic_file <- "raw/metab/COPDGene_P2_LT20miss_knnImp_NoOut_metabs_20211021.csv"
covariates_file <- "raw/pheno/COPDGene_P1P2P3_Flat_SM_NS_Sep24.txt"

# Also load list of individuals for accuracy to subset
ids <- read_tsv("data/genetics/id_list_genetics_and_metab_split.txt",
                col_names = c("fid", "iid"))

# Read in data
gene.exp <- read.csv(omic_file) %>%
  column_to_rownames("sid")

covar.data <- read_tsv(covariates_file) %>%
  dplyr::select(sid, Age_P2, gender, smoking_status_P2, ATS_PackYears_P2,
                BMI_P2, Pred_FEV1_P2) %>%
  column_to_rownames("sid")
# TO NOTE: this throws warnings, but doesn't affect the columns I need
    # _P2 indicates phase 2 visit, which is visit with metabolomics
  # Age = Age_P2, age at phase 2 visit?
  # Sex = gender, 1= male, 2 = female
  # Smoking status = smoking_status_P2
  # Pack years = ATS_PackYears_P2
  # BMI = BMI_P2
  # FEV1pp = FEV1pp_post_P2 (FEV1 % pred, post-BD, confirmed with Katerina)

# Subset individuals and check completeness of above variables
covar.data.filt <- covar.data %>%
  dplyr::filter(rownames(.) %in% ids$iid)
gene.exp.filt <- gene.exp %>%
  dplyr::filter(rownames(.) %in% ids$iid)

sum(is.na(covar.data.filt))  # 13
sum(is.na(covar.data.filt$Age_P2))  # 0
sum(is.na(covar.data.filt$gender))  # 0
sum(is.na(covar.data.filt$smoking_status_P2))  # 3
sum(is.na(covar.data.filt$ATS_PackYears_P2))  # 4
sum(is.na(covar.data.filt$BMI_P2))  # 3, TODO: could calculate?
sum(is.na(covar.data.filt$Pred_FEV1_P2))  # 3

covar.check <- covar.data.filt %>%
  dplyr::filter(is.na(ATS_PackYears_P2))
# All NAs are for same 4 people

# Run option 1
output_prefix <- "data/metab/COPDGene_P2_metabs_age_sex_covar"
covar.data.filt <- covar.data.filt %>%
  dplyr::select(Age_P2, gender)
sum(is.na(covar.data.filt))  # 0

# Run option 1, subset to individuals for option 2
# output_prefix <- "data/metab/COPDGene_P2_metabs_age_sex_covar_subset"
# covar.data.filt <- na.omit(covar.data.filt)  # 953
# covar.data.filt <- covar.data.filt %>%
#   dplyr::select(Age_P2, gender)

# Run option 2
# output_prefix <- "data/metab/COPDGene_P2_metabs_full_covar_subset"
# covar.data.filt <- na.omit(covar.data.filt)  # 953
# gene.exp.filt <- gene.exp.filt %>%
#   dplyr::filter(rownames(.) %in% rownames(covar.data.filt))

# Order omics and covariates
gene.exp.filt <- gene.exp.filt[rownames(covar.data.filt), ]
identical(rownames(covar.data.filt), rownames(gene.exp.filt))

# Make a copy of the gene.exp df and fill in with the residuals
expression <- gene.exp.filt

# Run MLR
for (i in 1:length(colnames(gene.exp.filt))) {
  fit <- lm(gene.exp.filt[,i] ~ as.matrix(covar.data.filt))
  expression[,i] <- fit$residuals
}

# Reformat for accuracy calculation
expression <- rownames_to_column(expression, "sid")

# Write results
write_tsv(expression, paste0(output_prefix, "_adj_resid.txt"))


