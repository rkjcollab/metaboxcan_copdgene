# SDS 20250109

# Create GCTA genetic relatedness matrix (GRM) from chip data in PLINK1.9.

plink_input="$1"
gcta_dir="$2"
maf_filter=$3
grm_filter=$4

# Create matrix with filter MAF < 0.05
plink_name=$(basename "$plink_input")
gcta-1.94.1 --bfile "$plink_input" \
    --make-grm --maf $maf_filter \
    --out ${gcta_dir}/${plink_name}_maf${maf_filter}_grm

# Filter matrix to a GRM cutoff of 0.05, removing related individuals
gcta-1.94.1 --grm ${gcta_dir}/${plink_name}_maf${maf_filter}_grm \
--grm-cutoff $grm_filter --make-grm \
--out ${gcta_dir}/${plink_name}_maf${maf_filter}_grm${grm_filter}
