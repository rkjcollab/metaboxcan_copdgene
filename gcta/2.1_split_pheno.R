#!/usr/bin/env Rscript

# Script to split single metabolite file into one file for each metabolite
# right before GCTA analysis. Each metabolite is given as the pheno input.
# First two columns must be family ID (FID), then individual ID (IID).
# Heritability analysis is being run on the observed metabolite values, not the
# predicted metabolite values.

# Setup -----------------------------------------------------------------

# Load packages
library(tidyverse)
library(doParallel)
library(foreach)

# Catch arguments
args <- commandArgs(trailingOnly = TRUE)
metab_path <- args[1]
output_path <- args[2]

# Local testing
# setwd(paste0(Sys.getenv("RKJCOLLAB"), "/Collabs/metaboxcan"))
# metab_path <- "raw/metab/COPDGene_P2_LT20miss_knnImp_NoOut_metabs_20211021.csv"
# output_path <- "data/gcta"

# Load data --------------------------------------------------------------------

# Load metabolite data
metab <- read_csv(metab_path)

# Update to ID format needed by GCTA
metab <- metab %>%
  dplyr::rename(IID = sid) %>%
  dplyr::mutate(FID = IID) %>%
  dplyr::relocate(FID, .before = IID)

# Make pheno files -------------------------------------------------------------

cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Get list of metab names
phen_names <- metab %>%
  select(-FID, -IID) %>%
  colnames(.)

# Make pheno file for each metab name
foreach(x = phen_names, .packages = c("dplyr"), .verbose = T) %dopar% {
  # Get one metab name at a time
  phen_indv <- metab %>% select(FID, IID, all_of(x))
  
  # Write out each single metab file
  output_file <- paste0(output_path, "/", x, ".phen")
  write.table(
    phen_indv,
    file = output_file,
    sep = " ", row.names = F, col.names = F, quote = F
  )
}

stopCluster(cl)
