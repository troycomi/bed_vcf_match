#!/bin/bash

#SBATCH --time=1-0
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/stucci_temp/processing/thin_%A_%a.out"
#SBATCH --array=1-22
#SBATCH --mem=1G
#SBATCH -A eeb

#here are the vcf files for modern humans
modern_vcf=/tigress/AKEY/akey_vol2/serenatu/analysis/rampa/filteredVCF_FINAL_OCT25/new.19.05.2017/dataset/merged_rps_1000g_png27

#vcf files for Altai Neandertal genome
altai_vcf=/tigress/AKEY/akey_vol1/home/bvernot/archaic_exome/data/neanderthal_altai_vcfs/2014.09.29/filtered_vcfs_altai_mpi_minimal_filters.noheader

#vcf files for the Denisovan genome
deni_vcf=/tigress/AKEY/akey_vol1/home/bvernot/archaic_exome/data/denisova_vcfs/2013.06.18/filtered_vcfs_den_mpi_minimal_filters

#here are the bed files which refer to the individuals in the vcf file. These are three columns that include #chr #start #end of a genomic region
bed_files=/tigress/AKEY/akey_vol1/home/serenatu/BACKUP_UW_Aug.2017/serenatu/STAR_2017/output/analysis/ALL.CALLS.AUG.2017
out_dir=/tigress/tcomi/stucci_temp/processing
vcftools=/tigress/AKEY/akey_vol1/software/Vcftools/vcftools-v0.1.15/bin/vcftools
bedtools=/tigress/AKEY/akey_vol1/software/bedtools2/bin/mergeBed
filter_bed=/tigress/AKEY/akey_vol2/serenatu/analysis/rampa/filteredVCF_FINAL_OCT25/filters_Jan26/FINAL_archaic_paper/Keep_by_Chr

module load anaconda3
conda activate bed2vcf

set -euo pipefail

chr=$SLURM_ARRAY_TASK_ID

zcat $modern_vcf/merged_rampa_1000g_png_phased_20.05.2017_chr${chr}.vcf.gz | \
python thin_vcf.py \
    --bed_file $filter_bed/keep_intersect_rps_arch_chr${chr}.bed.sorted.bed \
    --merge 10 \
    --individuals_file indivs.txt | \
gzip -c > $out_dir/chr${chr}.vcf.gz

exit

$bedtools \
    -i $filter_bed/keep_intersect_rps_arch_chr${chr}.bed.sorted.bed \
    -d 10 \
    > $out_dir/chr${chr}.temp.bed

$vcftools \
    --recode \
    --gzvcf $modern_vcf/merged_rampa_1000g_png_phased_20.05.2017_chr${chr}.vcf.gz \
    --bed $out_dir/chr${chr}.temp.bed \
    --keep indivs.txt \
    --stdout | gzip -c > $out_dir/chr${chr}.vcf.gz \

rm $out_dir/chr${chr}.temp.bed
