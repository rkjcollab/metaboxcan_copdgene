# SDS 20250110

# Filter PLINK input file to only individuals with metabolites.

plink_input="${RKJCOLLAB}/Collabs/metaboxcan/raw/genetics/CG10k_NHW_hg19_Oct2017/CG10k_NHW_hg19_Oct2017"
    # formatted as FID, IID
plink_out_dir="${RKJCOLLAB}/Collabs/metaboxcan/data/genetics/chip"
id_list="${RKJCOLLAB}/Collabs/metaboxcan/data/genetics/id_list_genetics_and_metab.txt"
    # formatted as FID_IID

# Make split version of ID list
id_list_new="${id_list%.*}"
awk '{print $1}' "${id_list}" | cut -d'_' -f1 > "tmp_fid.txt"
awk '{print $1}' "${id_list}" | cut -d'_' -f2 > "tmp_iid.txt"
paste -d'\t' "tmp_fid.txt" "tmp_iid.txt" > "${id_list_new}_split.txt"

# Filter PLINK input
plink_name=$(basename "$plink_input")
plink2 --bfile "$plink_input" \
    --keep "${id_list_new}_split.txt" \
    --make-bed --out "${plink_out_dir}/${plink_name}_filt"

# Clean up
rm tmp_*