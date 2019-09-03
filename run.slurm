#!/bin/bash

#SBATCH --time=3-0
#SBATCH --mem=24G
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/stucci_temp/bed_%A_%a.out"
#SBATCH -A eeb
#SBATCH --array=1-4

#here are the vcf files for modern humans
#modern_vcf=/tigress/AKEY/akey_vol2/serenatu/analysis/rampa/filteredVCF_FINAL_OCT25/new.19.05.2017/dataset/merged_rps_1000g_png27

#vcf files for Altai Neandertal genome
altai_vcf=/tigress/AKEY/akey_vol1/home/bvernot/archaic_exome/data/neanderthal_altai_vcfs/2014.09.29/filtered_vcfs_altai_mpi_minimal_filters.noheader

#vcf files for the Denisovan genome
deni_vcf=/tigress/AKEY/akey_vol1/home/bvernot/archaic_exome/data/denisova_vcfs/2013.06.18/filtered_vcfs_den_mpi_minimal_filters

#here are the bed files which refer to the individuals in the vcf file. These are three columns that include #chr #start #end of a genomic region
bed_files=/tigress/AKEY/akey_vol1/home/serenatu/BACKUP_UW_Aug.2017/serenatu/STAR_2017/output/analysis/ALL.CALLS.AUG.2017

module load anaconda3
conda activate bed2vcf

set -euo pipefail

## NOTE ## Need to make output dirs before running
case $SLURM_ARRAY_TASK_ID in
1)  # UV altai
    ./bed2vcf.py \
        --bed_output \
        --canc_correction \
        --bed_files $bed_files/UV*.PNG.*_hap?.bed.merged.bed \
        --output_dir /tigress/tcomi/stucci_temp/altai_CAnc_inf/ \
        --modern_vcfs /tigress/tcomi/stucci_temp/processing/chr{chr}.vcf.gz \
        --archaic_vcfs $altai_vcf/chr{chr}.altai_neand_mpi_minimal_filtered_lowqual.vcf.gz
    ;;

2)  # UV den
    ./bed2vcf.py \
        --bed_output \
        --canc_correction \
        --bed_files $bed_files/UV*.PNG.*_hap?.bed.merged.bed \
        --output_dir /tigress/tcomi/stucci_temp/denisova_2013.06.18_CAnc_inf/ \
        --modern_vcfs /tigress/tcomi/stucci_temp/processing/chr{chr}.vcf.gz \
        --archaic_vcfs $deni_vcf/chr{chr}.den_filtered.vcf.gz \
    ;;

3)  # RPS CHB altai
    ./bed2vcf.py \
        --bed_output \
        --canc_correction \
        --bed_files $bed_files/RPS*.RPS*_hap?.bed.merged.bed \
        --bed_files $bed_files/NA*.CHB*_hap?.bed.merged.bed \
        --output_dir /tigress/tcomi/stucci_temp/altai_CAnc_inf/ \
        --modern_vcfs /tigress/tcomi/stucci_temp/processing2/chr{chr}.vcf.gz \
        --archaic_vcfs $altai_vcf/chr{chr}.altai_neand_mpi_minimal_filtered_lowqual.vcf.gz
    ;;

4)  # RPS CHB den
    ./bed2vcf.py \
        --bed_output \
        --canc_correction \
        --bed_files $bed_files/RPS*.RPS*_hap?.bed.merged.bed \
        --bed_files $bed_files/NA*.CHB*_hap?.bed.merged.bed \
        --output_dir /tigress/tcomi/stucci_temp/denisova_2013.06.18_CAnc_inf/ \
        --modern_vcfs /tigress/tcomi/stucci_temp/processing2/chr{chr}.vcf.gz \
        --archaic_vcfs $deni_vcf/chr{chr}.den_filtered.vcf.gz \
    ;;
esac
