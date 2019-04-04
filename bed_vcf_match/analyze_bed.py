'''
analyze_bed

Module for performing analyses with bedfiles
'''


import pandas as pd
import numpy as np
from typing import TextIO, List


CHROMOSOME = 0
START = 1
END = 2


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
                     *archaic_vcfs: pd.DataFrame) -> List[float]:
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
    return joined.astype({'archaic': int})
