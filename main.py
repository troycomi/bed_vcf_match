#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
bed_vcf_match

Calculate archaic match percent from provided bed file
'''

import argparse
from typing import List
from itertools import chain
from bed_vcf_match.read_vcf import import_vcf
from bed_vcf_match.analyze_bed import read_bed, summarize_region
import pandas as pd
import gzip
import os


def main():
    args = read_args()

    if not any([args.save_db, args.vcf_output, args.bed_output]):
        return  # nothing to do

    modern_db = None

    # load old db
    if args.load_db:
        modern_db = pd.read_pickle(args.load_db)

    # read in modern vcfs
    print('reading modern vcfs', flush=True)
    for vcf in args.modern_vcfs:
        if os.path.splitext(vcf)[1] == '.gz':
            reader = gzip.open(vcf, 'rt')
        else:
            reader = open(vcf)

        modern_db = import_vcf(reader,
                               modern_db,
                               check_phasing=True)
        reader.close()

    # save new vcfs
    if args.save_db:
        modern_db.to_pickle(args.save_db)

    if not any([args.vcf_output, args.bed_output]):
        return  # nothing left to do

    # read in archaic files
    print('reading archaic vcfs', flush=True)
    archaic_dbs = []
    for vcf in args.archaic_vcfs:
        if os.path.splitext(vcf)[1] == '.gz':
            reader = gzip.open(vcf, 'rt')
        else:
            reader = open(vcf)

        archaic_dbs.append(import_vcf(reader,
                                      check_phasing=True))
        reader.close()

    print('processing bed files', flush=True)
    for bed_file in args.bed_files:
        in_file = os.path.split(bed_file)
        if args.out_dir is None:
            args.out_dir = in_file[0]
        # assume file is encoded as {indiv}.CODE.type_hap{haplotype}.bed...
        tokens = in_file[1].split('.')
        individual = tokens[0]
        haplotype = int(tokens[2][-1])
        outfile = os.path.join(args.out_dir,
                               in_file[1] + ".matched")

        # read in vcf entire vcf as pandas, on cache convert to current format
        # TODO cache the indiv, haplo merge with archaic
        with open(outfile) as bed_out:
            bed_lines = read_bed(bed_file)
            for line in bed_lines:
                bed_out.write(summarize_region(line.values[0],
                                               haplotype,
                                               individual,
                                               modern_db,
                                               *archaic_dbs))


def read_args(args: List[str] = None) -> argparse.Namespace:
    '''
    read in command line arguments, returning namespace object
    '''
    parser = argparse.ArgumentParser(description='Match bed files with vcfs')
    parser.add_argument('--save_db',
                        default=None,
                        help='If set to a valid filename, save the vcf '
                        'database as a pandas binary object. Valid extensions '
                        'include pkl, gz, bz2, zip or xz.'
                        )

    parser.add_argument('--load_db',
                        default=None,
                        help='If set to a valid filename, load the vcf '
                        'database from a previous execution. If set, no vcfs '
                        'are required for input.'
                        )

    parser.add_argument('--modern_vcfs',
                        default=None,
                        action='append',
                        nargs='*',
                        help='List of modern vcf files to add to database. '
                        'If files are provided with an existing database to '
                        'load, the vcf files will be added to the database.'
                        )

    parser.add_argument('--archaic_vcfs',
                        default=None,
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
