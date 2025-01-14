#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

print("This script merges GCTA output files into a single csv file.")
print("First trailing arg should be location of .hsq & .log files.")
print("Second trailing arg should be path/name for merged csv.")

library(tidyverse)

# Function to merge reml output
merge_reml <- function(path_in, path_out) {
  reml_files <- list.files(path = path_in, pattern = "*.hsq", full.names = TRUE)
  text_files <- list.files(path = path_in, pattern = "*.log")
  print(paste0("For file path: ", path_in))
  print(paste0("Number of .log files should match number of metabolites input: ", length(text_files)))
  print(paste0("Heritability converged for ", length(reml_files), " metabolites."))

  # Throw error if all files not distinct, or if NA file names
  if (n_distinct(reml_files) != length(reml_files)) stop("File names not distinct.")
  if (sum(is.na(reml_files)) != 0) stop("At least one file name is NA.")
  
  # Merge all metab results from one reml into one df
  # To note: the read_delim throws error, expected due to empty col 3 cells
  reml_raw <- read_delim(reml_files, id="name", show_col_types = FALSE)
  
  # Format for export
  # To note: assumes all inputs will have metab names followed by "_data..."
  # but also works as is if name is just gene name.
  reml_exp <- reml_raw %>%
    pivot_wider(names_from = Source, values_from = c(Variance, SE)) %>%
    select(-`SE_logL`, -`SE_logL0`, -`SE_LRT`, -`SE_df`, -`SE_Pval`, -`SE_n`) %>%
    rename_all(~ sub("Variance_", "", .)) %>%
    rowwise() %>%
    mutate(metab = sub("_data.*", "", tail(str_split(name, "/")[[1]], 1))) %>%
    relocate(metab) %>%
    select(-name)
  
  # Export
  write_csv(reml_exp, paste0(path_out, ".csv"))  # all data
}

# Call function on path to merge reml output - all data
merge_reml(args[1], args[2])
