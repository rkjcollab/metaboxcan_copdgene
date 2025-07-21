#!/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --ntasks=6
#SBATCH --job-name="tm_r3_imp"
#SBATCH -o "/scratch/alpine/sslack@xsede.org/temp_data/copdgene_tm_r3_imp/job_2_run_%x.out"
#SBATCH -e "/scratch/alpine/sslack@xsede.org/temp_data/copdgene_tm_r3_imp/job_2_run_%x.err"
#SBATCH --account=amc-general
#SBATCH --time=08:00:00
#SBATCH --mem=200G
#SBATCH --qos=normal

# This is currently called inside of the limactl amd64.
# Set paths and options in config.yml file. Paths are automatically 
# relative to snakefile location.

# Uncomment step want to run through
# step="submit_initial_input"  # pre-imp QC
# step="submit_fix_strands"  # submit imp
step="concat_convert_to_plink" # unzip, clean, and merge

# Set RKJCOLLAB since not present in lima
RKJCOLLAB="/Users/slacksa/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver"

# Set number cores or get from SLURM
n_cores=6  #"$SLURM_NTASKS"

# Local limactl
# --bind ${RKJCOLLAB}/Collabs/metaboxcan/data/genetics/tm_r3_imp:/output_data \
apptainer exec \
    --writable-tmpfs \
    --bind /Users/slacksa/repos/imputation_snakemake:/repo \
    --bind /Users/slacksa/repos/metaboxcan_copdgene:/proj_repo \
    --bind ${RKJCOLLAB}/Collabs/metaboxcan/data/genetics/chip:/input_data \
    --bind /Users/slacksa/temp_data/metaboxcan/copdgene_tm_r3_imp:/output_data \
    /Users/slacksa/repos/imputation_snakemake/envs/topmed_imputation.sif \
    snakemake --rerun-triggers mtime --snakefile /repo/Snakefile \
        --configfile /proj_repo/config_tm_r3_imp.yml \
        --cores "$n_cores" "$step"
