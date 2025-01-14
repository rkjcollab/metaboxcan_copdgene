# SDS 20250113

# Script runs GCTA analysis

if [ "$#" -eq 0 ]
then
    echo "Usage: ${0##*/} <metab_phen>"
    echo "Script runs GCTA analysis. All paths are set in bash script,"
    echo "only metabolite pheno files passed by batch script."
    exit
fi

metab_path="gcta/metab_phen"
grm="gcta/grm/CG10k_NHW_hg19_Oct2017_filt_maf0.05_grm0.05"
out_dir="gcta/${SLURM_JOB_NAME}"

m_name=$(basename "$1" .phen)

# Running constrained REML analysis
gcta-1.94.1 --reml \
	--grm "$grm" \
	--pheno "${metab_path}/${1}" \
	--out "${out_dir}/${m_name}"