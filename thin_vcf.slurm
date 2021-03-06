#!/bin/bash

#SBATCH --time=0-1
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/stucci_temp/processing2/thin_%A_%a.out"
#SBATCH --array=1-22
#SBATCH --mem=500M
#SBATCH -A eeb

#here are the vcf files for modern humans
modern_vcf=/tigress/AKEY/akey_vol2/serenatu/analysis/rampa/filteredVCF_FINAL_OCT25/new.19.05.2017/dataset/merged_rps_1000g_png27

bed_files=/tigress/AKEY/akey_vol1/home/serenatu/BACKUP_UW_Aug.2017/serenatu/STAR_2017/output/analysis/ALL.CALLS.AUG.2017
out_dir=/tigress/tcomi/stucci_temp/processing2
vcftools=/tigress/AKEY/akey_vol1/software/Vcftools/vcftools-v0.1.15/bin/vcftools
bedtools=/tigress/AKEY/akey_vol1/software/bedtools2/bin/mergeBed
filter_bed=/tigress/AKEY/akey_vol2/serenatu/analysis/two.den/filters/sorted

module load anaconda3
conda activate bed2vcf

set -euo pipefail

chr=$SLURM_ARRAY_TASK_ID

zcat $modern_vcf/merged_rampa_1000g_png_phased_20.05.2017_chr${chr}.vcf.gz | \
python thin_vcf.py \
    --bed_file $filter_bed/rampa.final.sites.to.keep.chr${chr}.bed.sorted.bed \
    --individuals_file indivs2.txt | \
gzip -c > $out_dir/chr${chr}.vcf.gz

exit # The above should replace this

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
