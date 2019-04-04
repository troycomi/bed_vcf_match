#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
thin_vcf

Remove unneeded rows from vcf. Equivalent function to merge bedfile and
vcftools recode with keep individuals
'''

import argparse
import numpy as np
import sys
from typing import List, TextIO


def main():
    args = read_args()

    # read in and merge bedfile
    if args.bed_file is not None:
        with open(args.bed_file) as bed_file:
            bed = read_bed(bed_file, args.merge)
    else:
        bed = None

    if args.individuals_file is not None:
        with open(args.individuals_file) as indiv_file:
            indivs = read_indivs(indiv_file)
    else:
        indivs = None

    parser = line_parser(bed, indivs)

    for line in sys.stdin:
        print(parser.process_line(line), end='')


def read_args(args: List[str] = None) -> argparse.Namespace:
    '''
    read in command line arguments, returning namespace object
    '''
    parser = argparse.ArgumentParser(description='Filter vcfs from stdin')
    parser.add_argument('--bed_file',
                        default=None,
                        help='If set, retain only the sites within bed '
                        'file regions.  Assumes indexing is identical and '
                        'inclusive with start and end. Must be sorted.')

    parser.add_argument('--merge',
                        default=None,
                        type=int,
                        help='If set to an integer, merge bed regions within '
                        '<merge> away from eachother, inclusive.')

    parser.add_argument('--individuals_file',
                        default=None,
                        help='If set, only individuals from the file are kept.'
                        ' One individual per line of input file.')

    args = parser.parse_args(args)
    return args


def read_bed(reader: TextIO, merge: int = None) -> np.array:
    '''
    read in the opened file as a numpy array.
    If merge is provided, adjacent sites are merged if
    start_i - end_i-1 <= merge.
    Return only start and end sites, chromosome is assumed constant and equal
    '''
    result = np.loadtxt(reader,
                        dtype=int,
                        delimiter='\t',
                        usecols=(1, 2))
    if merge is not None:
        # find locations to keep (diff between end and start < merge)
        tokeep = np.concatenate(
            ([True],
             np.diff(result.flatten())[1::2] > merge))

        swap_inds = np.where(tokeep[:-1] != tokeep[1:])[0]
        # nothing to merge/swap
        if swap_inds.size > 0:
            # this happens when last ind needs to swap
            if swap_inds.size % 2 == 1:
                swap_inds = np.append(swap_inds, result.shape[0] - 1)
            result[swap_inds[0::2], 1] = result[swap_inds[1::2], 1]
            result = result[tokeep, :]
    return result


def read_indivs(reader: TextIO) -> List[str]:
    '''
    read in each line of the file, return as list with newline stripped
    '''
    return [line.rstrip() for line in reader]


def match_indivs(header: List[str], indivs: List[str]) -> List[int]:
    header = np.array(header)
    indivs = np.array(indivs)

    matches = np.char.equal(header[:, None], indivs[None, :])
    return np.where(matches)[0]


class line_parser():
    def __init__(self, bed: np.array, indivs: List[str]):
        self.bed = bed
        self.indivs = indivs
        if self.bed is None:
            self.parser = self.no_op
        else:
            self.bed += 1
            self.parser = self.filter_bed
        self.retain = 9
        self.bed_ind = 0

    def process_line(self, line: str) -> str:
        if line[1] == '#':
            return line

        if line[0] == '#':
            if self.indivs is not None:
                header = line.split('\t')
                self.indivs = match_indivs(header, self.indivs)
                self.indivs = np.concatenate((np.array(range(self.retain)),
                                             self.indivs))
                if self.bed is None:
                    self.parser = self.filter_indiv
                else:
                    self.parser = self.filter_both

                return self.filter_indiv(line)
            else:
                return line

        return self.parser(line)

    def no_op(self, line: str):
        return line

    def filter_indiv(self, line: str):
        line = line.split('\t')
        return '\t'.join([line[i] for i in self.indivs]).rstrip() + '\n'

    def filter_bed(self, line: str):
        pos = int(line.split('\t')[1])
        if self.bed_ind >= self.bed.shape[0]:
            return ''
        while pos >= self.bed[self.bed_ind, 1]:
            self.bed_ind += 1
            if self.bed_ind >= self.bed.shape[0]:
                return ''
        # by the time we are here, pos < end
        if pos >= self.bed[self.bed_ind, 0]:
            return line
        else:
            return ''

    def filter_both(self, line: str):
        tokens = line.split('\t')
        pos = int(tokens[1])
        if self.bed_ind >= self.bed.shape[0]:
            return ''
        while pos >= self.bed[self.bed_ind, 1]:
            self.bed_ind += 1
            if self.bed_ind >= self.bed.shape[0]:
                return ''
        # by the time we are here, pos < end
        if pos >= self.bed[self.bed_ind, 0]:
            return '\t'.join([tokens[i] for i in self.indivs]).rstrip() + '\n'
        else:
            return ''


if __name__ == "__main__":
    main()
