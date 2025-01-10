# SDS 20250109

# Create GCTA genetic relatedness matrix (GRM) from chip data in PLINK1.9.

plink_input="${RKJCOLLAB}/Collabs/metaboxcan/raw/genetics/CG10k_NHW_hg19_Oct2017/CG10k_NHW_hg19_Oct2017.bed"
gcta_dir="${RKJCOLLAB}/Collabs/metaboxcan/data/gcta"
maf_filter=0.05
grm_filter=0.05

# Create matrix with filter MAF < 0.05
plink_name=$(basename "${plink_input%.*}")
gcta-1.94.1 --bfile "$plink_input" \
    --make-grm --maf $maf_filter \
    --out ${gcta_dir}/${plink_name}_maf${maf_filter}_grm

# Filter matrix to a GRM cutoff of 0.05, removing related individuals
gcta-1.94.1 --grm ${gcta_dir}/${plink_name}_maf${maf_filter}_grm \
--grm-cutoff $grm_filter --make-grm \
--out ${gcta_dir}/${plink_name}_maf${maf_filter}_grm${grm_filter}