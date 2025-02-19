# SDS 20240106

# Calculate prediction accuracy.

# Initially ran with metabolite input as provided. Now re-running with values
# adjusted for different sets of covariates:

# Running with two different sets of covariates:
  # 1. age and sex (most closely matches paper draft)
  # 2. age, sex, smoking status, pack years, BMI, and FEV1pp (suggested
    # by Katerina)

# Subset to the 957 individuals that calculating accuracy in before running
# model to get residuals. To do this and keep all covars for option 2, reduce
# to 953 individuals. Will run option 1 with all 957 and with 953 to be
# comparable.

# Setup ------------------------------------------------------------------------

setwd(paste0(Sys.getenv("RKJCOLLAB"), "/Collabs/metaboxcan"))

library(tidyverse)
library(data.table)
library(conflicted)

pred_gex <- read_delim(
  "data/metsim_prediction/filt/COPDGene_metsim_metaboxcan_predict.txt")

# Unadjusted
# metab <- read_csv(
#   "raw/metab/COPDGene_P2_LT20miss_knnImp_NoOut_metabs_20211021.csv")
# output_prefix <- ""

# Covar option 1
metab <- read_tsv(
  "data/metab/COPDGene_P2_metabs_age_sex_covar_adj_resid.txt")
output_prefix <- "age_sex_covar_adj_"

# Covar option 1, subset to same individuals as option 2
# metab <- read_tsv(
#   "data/metab/COPDGene_P2_metabs_age_sex_covar_subset_adj_resid.txt")
# output_prefix <- "age_sex_covar_subset_adj_"

# Covar option 2
# metab <- read_tsv(
#   "data/metab/COPDGene_P2_metabs_full_covar_subset_adj_resid.txt")
# output_prefix <- "full_covar_subset_adj_"


metab_info <- read_csv(
  "raw/metab/COPDGene_P2_MetaboliteInformation_20211021.csv")

pheno <- read_tsv(
  "raw/pheno/COPDGene_P1P2P3_Flat_SM_NS_Sep24.txt")
# TO NOTE: this throws warnings, but doesn't affect the column I need (gender)
# 1 = male, 2 = female

# Map and prep -----------------------------------------------------------------

# Mapping:
# metabolite data: metab_id --> CHEM_ID
# predicted expression: remove C from IDs, match to CHEM_ID

# Update metabolite data
old  <- metab_info$metab_id
new  <- as.character(metab_info$CHEM_ID)
metab_mod <- setnames(
  metab, old, new, skip_absent = T)

# Update predicted metabolite expression
pred_gex_mod <- pred_gex
colnames(pred_gex_mod) <- gsub("^C", "", colnames(pred_gex_mod))

pred_gex_mod <- pred_gex_mod %>%
  dplyr::mutate(ID = gsub("^.+_", "", IID)) %>%
  dplyr::relocate(ID, .after = IID) %>%
  dplyr::select(-FID, -IID)

# Merge data -------------------------------------------------------------------

metab_mod$pop = "obs.expr"
pred_gex_mod$pop = "pred.expr"

# Subset to only overlapping
metab_ids <- intersect(colnames(metab_mod), colnames(pred_gex_mod))  # 790
n_distinct(metab_ids)
  # 789 without pop
indv_ids <- intersect(metab_mod$sid, pred_gex_mod$ID)  # 957/953
n_distinct(indv_ids)

metab_mod_filt <- metab_mod %>%
  dplyr::filter(sid %in% indv_ids) %>%
  dplyr::select(sid, all_of(metab_ids)) %>% # 957/953 by 791 (with pop and id col)
  dplyr::rename(ID = sid)

pred_gex_mod_filt <- pred_gex_mod %>%
  dplyr::filter(ID %in% indv_ids) %>%
  dplyr::select(ID, all_of(metab_ids))  # 957/953 by 791

# Merge
merge <- rbind(metab_mod_filt, pred_gex_mod_filt) %>%
  dplyr::relocate(pop, .after = ID)

# Calculate accuracy -----------------------------------------------------------

# Calculating Spearman and Pearson correlation, but likely only using Pearson
# correlation based on experience in Predixcan project, see complete notes in
# calc_geuvadis_prediction_accuracy.Rmd:
  # Checking the Predixcan paper by Gamazon (PMID 6258848). From the methods:
  # "To assess performance, we used the square of the Pearson correlation, R2,
  # between predicted and observed expression levels."
  # UPDATE: Since now analyzing residuals, thought Spearman would no longer have
  # issues with ties, but still does.

# Function to get a data frame of each metab, R2, and correlation
# df <- merge
metabs <- grep("ID|pop", colnames(merge), invert = T, value = T)

get_pred_vs_obs <- function(df, metabs) {
  result_df <- NULL
  
  for(i in 1:length(metabs)){
    # i <- 1
    m <- df %>%
      dplyr::select(ID, pop, metabs[i]) %>%
      spread(pop, metabs[i]) %>%
      na.omit()
    metab <<- metabs[i]
    
    test.p <- cor.test(m$obs.expr, m$pred.expr, method="pearson",
                       alternative="two.sided")
    test.s <- cor.test(m$obs.expr, m$pred.expr, method="spearman",
                       alternative="two.sided")
    
    row <- tibble(metab = metabs[i],
                  pearson.cor = test.p$estimate,
                  pearson.r2 = test.p$estimate^2,
                  t.stat = test.p$statistic,
                  pearson.pval = test.p$p.value,
                  conf.int = paste0(round(test.p$conf.int[1],4),
                                    ", ", round(test.p$conf.int[2],4)),
                  df = test.p$parameter,
                  hyp = test.p$alternative,
                  spearman.cor = test.s$estimate,
                  spearman.r2 = test.s$estimate^2,
                  spearman.pval = test.s$p.value)
    
    result_df <- rbind(result_df, row)
  }
  return(result_df)
  
}

  
result <- get_pred_vs_obs(merge, metabs)
# As expected, Spearman test threw many "cannot compute exact p-value with ties"
# errors, likely due to identical predicted expression as seen in Predixcan
# project.


# Calculate accuracy stratified by sex -----------------------------------------

pheno_filt <- pheno %>% dplyr::select(sid, gender)
merge_mod <- left_join(
  merge,
  pheno_filt,
  by = c("ID" = "sid"))
sum(is.na(merge_mod$gender))  # 0 NA

merge_mod_m <- merge_mod %>%
  dplyr::filter(gender == 1)
n_distinct(merge_mod_m$ID)  # 489/487 males
merge_mod_f <- merge_mod %>%
  dplyr::filter(gender == 2)
n_distinct(merge_mod_f$ID)  # 468/466 females

result_m <- get_pred_vs_obs(merge_mod_m, metabs)
result_f <- get_pred_vs_obs(merge_mod_f, metabs)

# Write out --------------------------------------------------------------------

write_tsv(merge, paste0(
  "data/metsim_prediction/filt/COPDGene_", output_prefix, "pred_v_obs_metsim.txt"))

write_tsv(result, paste0(
  "data/metsim_prediction/filt/COPDGene_", output_prefix, "metsim_prediction_performance.txt"))

write_tsv(result_m, paste0(
  "data/metsim_prediction/filt/COPDGene_", output_prefix, "metsim_prediction_performance_males.txt"))

write_tsv(result_f, paste0(
  "data/metsim_prediction/filt/COPDGene_", output_prefix, "metsim_prediction_performance_females.txt"))
