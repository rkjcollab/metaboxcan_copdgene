#!/bin/bash

# Need a SNPs only VCF file for input into PredictDB tool on Seven Bridges.
# Also want to apply MAF filter of 0.01 - stricter than used for general TOPMed
# imputation (MAF 0).

cd "${RKJCOLLAB}/Collabs/metaboxcan"

plink2 --pfile data/genetics/tm_r3_imp/imputed_clean_maf0_rsq0.3/chr_all_concat \
    --snps-only --export vcf bgz --maf 0.01 \
    --out data/genetics/tm_r3_imp/imputed_clean_maf0.01_rsq0.3/chr_all_concat_snps_only
    