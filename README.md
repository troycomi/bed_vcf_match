# BED to VCF match

> Correlate calls from bed files to vcf to determine percent match to archaic
genotypes.

## Installation
All dependencies are specified in the environment.yaml file.  A conda
environment can be created with
```
conda create --name bed2vcf --file environment.yaml
```
Once created, activate the environment with: `conda activate bed2vcf`

## Tests
All unit tests can be run within the conda environment by calling `pytest`

## Usage
There are two main functions in the top level directory

### thin\_vcf.py
This function performs filtering of vcfs based on individual names or bed files
to retain.  The script thin\_vcf.slurm shows a submission along with the
vcftools implementation.  Leveraging zcat and gzip to perform compression and
using the sorted nature of both files, this runs more quickly than the vcftools
version.

### bed2vcf.py
This performs the main analysis and is run with run.slurm.  Bed files are 
provided with the `--bed_files` flag and outputs are written to the
`--output_dir` with the same base filename as the inputs, appended with 
".matched".  Vcfs should be provided for each chromosome and labeled with a 
wildcard for python formatting.  E.g. if the files are named chrom\_1\_filter.vcf
the command line expects the input chrom\_\{chr\}\_filter.vcf and will expand
to each chromosome in \[1, 22\]

### Output Format
Currently only bed output is supported.  Each row corresponds to a row in the
input bed file.  The columns are:
1. Chromosome
2. Start position
3. End position
4. Number of entries in the modern VCF within this region.
5. Number of derived variants within the region for a haplotype
6. Number of derived variants in the modern VCF joined with archaic.
Intersection of sites, not necessarily matching variants. 
7. Number of derived variants in the archaic VCF within the region (number of alleles/2)
8. Number of derived variants in archaic which match a derived variant
in the modern vcf.

## License

MIT © [Troy Comi](https://github.com/troycomi)
