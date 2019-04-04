#!/bin/bash

#SBATCH --time=0-8
#SBATCH --mem=32G
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/stucci_temp/bed_%A.out"

#here are the vcf files for modern humans
modern_vcf=/tigress/AKEY/akey_vol2/serenatu/analysis/rampa/filteredVCF_FINAL_OCT25/new.19.05.2017/dataset/merged_rps_1000g_png27

#vcf files for Altai Neandertal genome
altai_vcf=/tigress/AKEY/akey_vol1/home/bvernot/archaic_exome/data/neanderthal_altai_vcfs/2014.09.29/filtered_vcfs_altai_mpi_minimal_filters.noheader

#vcf files for the Denisovan genome
deni_vcf=/tigress/AKEY/akey_vol1/home/bvernot/archaic_exome/data/denisova_vcfs/2013.06.18/filtered_vcfs_den_mpi_minimal_filters

#here are the bed files which refer to the individuals in the vcf file. These are three columns that include #chr #start #end of a genomic region
bed_files=/tigress/AKEY/akey_vol1/home/serenatu/BACKUP_UW_Aug.2017/serenatu/STAR_2017/output/analysis/ALL.CALLS.AUG.2017
out_dir=/tigress/tcomi/stucci_temp/

module load anaconda3
conda activate bed2vcf

set -euo pipefail

# TODO always submit haplotypes to same function
./main.py \
    --bed_files $bed_files/UV5??.PNG.*_hap?.bed.merged.bed \
    --bed_output \
    --output_dir $out_dir \
    --modern_vcfs $modern_vcf/merged_rampa_1000g_png_phased_20.05.2017_chr*.vcf.gz \
    --archaic_vcfs $altai_vcf/chr[0-9].altai_neand_mpi_minimal_filtered_lowqual.vcf.gz \
    --archaic_vcfs $altai_vcf/chr[0-9][0-9].altai_neand_mpi_minimal_filtered_lowqual.vcf.gz \
    #--modern_vcfs $modern_vcf/merged_rampa_1000g_png_phased_20.05.2017_chr22.vcf.gz \
    #--archaic_vcfs $altai_vcf/chr22.altai_neand_mpi_minimal_filtered_lowqual.vcf.gz \
