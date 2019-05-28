'''
analyze_bed

Module for performing analyses with bedfiles
'''


import pandas as pd
import numpy as np
import os
from typing import TextIO, List, Dict, Tuple


CHROMOSOME = 0
START = 1
END = 2


class bed_structure():
    def __init__(self, filename, out_dir=None):
        with open(filename, 'r') as reader:
            self.bed = structure_bed(reader)

        # split filename into path and file
        self.filename = filename
        t = os.path.split(filename)

        # pull out individual, haplotype
        tokens = t[1].split('.')
        self.individual = tokens[0]
        self.haplotype = int(tokens[2][-1])

        # setup output file
        if out_dir is None:
            out_dir = t[0]
        outfile = os.path.join(out_dir, t[1] + ".matched")
        self.writer = open(outfile, 'w')

    def process_chrom(self, chromosome: int, modern_db, archaic_db):
        chrm = str(chromosome)
        if chrm not in self.bed:
            return
        for start, end in self.bed[chrm]:
            self.writer.write(summarize_region(
                [chromosome, start, end],
                self.haplotype,
                self.individual,
                modern_db,
                archaic_db))

    def close(self):
        self.writer.close()


def structure_bed(reader: TextIO) -> Dict[str, List[Tuple[int, int]]]:
    '''
    read in the bed file, returning a dictionary keyed by chromosome
    with a list of (start, end) tuples
    '''
    result = {}
    for line in reader:
        chrom, start, end = line.split()
        if chrom not in result:
            result[chrom] = []
        result[chrom].append((int(start), int(end)))

    return result


def read_bed(reader: TextIO) -> pd.io.parsers.TextFileReader:
    '''
    Read in bed file returning a pandas TextFileReader to iterator
    over the file.
    '''
    names = [''] * 3
    names[CHROMOSOME] = 'chrom'
    names[START] = 'start'
    names[END] = 'end'
    return pd.read_csv(reader, delimiter='\t', header=None,
                       names=names,
                       chunksize=1)


def summarize_region(bed_line: List[int],
                     haplotype: int,
                     individual: str,
                     modern_vcf: pd.DataFrame,
                     *archaic_vcfs: pd.DataFrame) -> str:
    '''
    Given a line from a bed file, the individual and haplotype,
    look up the corresponding region in the modern vcf database,
    match chromosome, site, and genotype with archaic vcfs.
    Return the bed line with added columns for the # of sites, number of
    modern and archaic variants, number of matches and fraction of sites
    matching for each archaic vcf provided.
    '''
    rows = filter_modern_db(modern_vcf,
                            bed_line[CHROMOSOME],
                            bed_line[START],
                            bed_line[END],
                            haplotype,
                            individual)

    # same for each archaic vcf
    line = '\t'.join([str(b) for b in bed_line])
    sites = len(rows)
    line += f'\t{sites}'  # number of sites
    line += f'\t{np.sum(rows["variant"])}'  # number of variants

    for archaic_vcf in archaic_vcfs:
        joined = join_vcf(rows, archaic_vcf)
        line += f'\t{np.sum(joined["archaic"])/2}'  # number archaic variants
        matches = np.sum(joined['archaic'] * joined['variant'])/2
        line += f'\t{matches}'
        if sites != 0:
            fraction = round(matches / sites, 6)
            line += f'\t{fraction}'
        else:
            line += '\tnan'

    return line + '\n'


def filter_modern_db(modern_vcf: pd.DataFrame,
                     chrom: int,
                     start: int,
                     end: int,
                     haplotype: int,
                     individual: str) -> pd.DataFrame:
    '''
    Get matching rows in modern database. Find start < position <= end
    '''

    result = modern_vcf.loc[
        (modern_vcf['chrom'] == chrom) &
        (modern_vcf['pos'] > start) &
        (modern_vcf['pos'] <= end),
        ['chrom', 'pos', 'ref', 'alt', individual]]
    result.rename({individual: 'variant'}, axis='columns', inplace=True)
    result.dropna(inplace=True)

    replace = {'0..': 0, '1..': 1}  # haplotype == 1
    if haplotype == 2:
        replace = {'..0': 0, '..1': 1}
    result.variant.replace(regex=replace, inplace=True)

    return result


def join_vcf(modern: pd.DataFrame, archaic: pd.DataFrame) -> pd.DataFrame:
    '''
    Merge the modern and archaic vcf dataframes, adding a new column of the
    archaic variant to the modern vcf.  If the chrom, pos, ref and alt do not
    match, a 0 is inserted.
    '''
    joined = pd.merge(modern,
                      archaic.rename(index=str,
                                     columns={'variant': 'archaic'}),
                      how='left',
                      on=['chrom', 'pos', 'ref', 'alt'])
    joined.fillna(0, inplace=True)
    # TODO add in canc filtering here if it is a column in result
    if 'CAnc' in joined.columns:
        match_alt = joined.CAnc == joined.alt
        joined.loc[match_alt, 'archaic'] = 2 - joined.loc[match_alt, 'archaic']
        joined.loc[match_alt, 'variant'] = 1 - joined.loc[match_alt, 'variant']
        match_neither = (joined.CAnc != joined.alt) & \
            (joined.CAnc != joined.ref)
        joined.loc[match_neither, ['archaic', 'variant']] = 0
        joined = joined.drop(columns='CAnc')
    return joined.astype({'archaic': int})
