#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
bed2vcf

Calculate archaic match percent from provided bed file
'''

import argparse
from typing import List
from itertools import chain
from bed_vcf_match.read_vcf import import_vcf, import_archaic_vcf
from bed_vcf_match.analyze_bed import bed_structure
import gzip
import os


def main():
    args = read_args()

    # read in bed files, building output files and bed structure
    beds = []
    for bed_file in args.bed_files:
        beds.append(bed_structure(bed_file, args.output_dir))

    indivs = list(set([bed.individual for bed in beds]))
    print(f'found {len(indivs)} individuals')
    print(f'found {len(beds)} bed files')

    for chrm in range(1, 23):
        print(f'starting chromosome {chrm}')
        vcf = args.modern_vcfs[0].format(chr=chrm)
        print(os.path.split(vcf)[1], flush=True)
        if os.path.splitext(vcf)[1] == '.gz':
            reader = gzip.open(vcf, 'rt')
        else:
            reader = open(vcf)
        modern_db = import_vcf(reader,
                               check_phasing=True,
                               individuals=indivs)
        reader.close()

        vcf = args.archaic_vcfs[0].format(chr=chrm)
        print(os.path.split(vcf)[1], flush=True)
        if os.path.splitext(vcf)[1] == '.gz':
            reader = gzip.open(vcf, 'rt')
        else:
            reader = open(vcf)

        archaic_db = import_archaic_vcf(reader)
        reader.close()

        print('starting bed output...')
        for bed in beds:
            bed.process_chrom(chrm, modern_db, archaic_db)
        print(f'finished chromosome {chrm}')

    for bed in beds:
        bed.close()
    print('done!')


def read_args(args: List[str] = None) -> argparse.Namespace:
    '''
    read in command line arguments, returning namespace object
    '''
    parser = argparse.ArgumentParser(description='Match bed files with vcfs')

    parser.add_argument('--modern_vcfs',
                        default=[],
                        action='append',
                        nargs='*',
                        help='List of modern vcf files to add to database. '
                        'If files are provided with an existing database to '
                        'load, the vcf files will be added to the database.'
                        )

    parser.add_argument('--archaic_vcfs',
                        default=[],
                        action='append',
                        nargs='*',
                        help='List of archaic vcf files to add to database.'
                        'Archaic sites are merged with modern such that only '
                        'sites found in the existing database are retained.'
                        )

    parser.add_argument('--bed_files',
                        default=None,
                        action='append',
                        nargs='*',
                        help='List of bed files to query against database.'
                        )

    parser.add_argument('--vcf_output',
                        action='store_true',
                        help='If set, will output a vcf-like file for each '
                        'input bed file.'
                        )

    parser.add_argument('--bed_output',
                        action='store_true',
                        help='If set, will output a bed-like file for each '
                        'input bed file.'
                        )

    parser.add_argument('--output_dir',
                        default=None,
                        help='The output directory to store all files. '
                        'If not set but output is requested, output will be '
                        'created in the same directory as the input files.'
                        )

    args = parser.parse_args(args)

    # need to flatten the file args since they are lists of lists
    # using multiple args and append
    to_flatten = ['bed_files', 'modern_vcfs', 'archaic_vcfs']
    for flatten in to_flatten:
        arg = getattr(args, flatten)
        if arg is not None:
            setattr(args, flatten, list(chain.from_iterable(arg)))
    return args


if __name__ == '__main__':
    main()
