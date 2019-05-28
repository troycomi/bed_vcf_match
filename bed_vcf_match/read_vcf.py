'''
read_vcf

Read a vcf file into pandas data frame
'''


import pandas as pd
from typing import TextIO, List
import numpy as np


def import_vcf(vcf_reader: TextIO,
               dataframe: pd.DataFrame = None,
               check_phasing: bool = False,
               individuals: List[str] = None) -> pd.DataFrame:
    '''
    Read in all lines of the provided, open vcf file and concatenate with
    provided pandas dataframe
    check_phasing: if true, raises value error for any unphased haplotype
    that is not ./.  If false, unphased haplotypes are converted to nan values,
    or collapsed.
    collapse_haplotype: if true, convert haplotypes to the sum of the
    genotypes.
    individuals: if specified, limit the imported data to only the provided
    individuals.  Individuals not found in the file raise value errors
    '''
    header_lines = 1  # 1-based indexing on error reporting
    for line in vcf_reader:
        header_lines += 1
        if line[1] == '#':  # comment string
            continue

        if line[0] == '#':  # header string
            header = line[1:].rstrip().split('\t')
            break

    indivs = header[9:]
    if individuals is not None:
        for indiv in individuals:
            if indiv not in indivs:
                raise ValueError(f'{indiv} not in file!')
        indivs = [indiv for indiv in set(individuals)]

    header = [h.lower() for h in header[:9]] + header[9:]
    usecols = [header[i] for i in [0, 1, 3, 4]] + indivs

    # precompute this as a dictionary for hopefully faster operations
    phase_decoder = {}
    for i in range(2):
        for j in range(2):
            phase_decoder[f'{i}|{j}'] = f'{i}|{j}'
            phase_decoder[f'{i}/{j}'] = np.nan
    phase_decoder['./.'] = 0

    new_frame = pd.read_csv(vcf_reader,
                            delimiter='\t',
                            header=None,
                            names=header,
                            usecols=usecols)

    new_frame = new_frame.loc[(new_frame.ref.str.len() == 1)
                              & (new_frame.alt.str.len() == 1)]

    new_frame[indivs] = new_frame[indivs].applymap(phase_decoder.get)

    if check_phasing:
        if new_frame.isna().any(axis=None):
            nanrow = new_frame[new_frame.isna().any(axis=1)].iloc[0]
            nanind = nanrow.loc[nanrow.isna()].index[0]
            raise ValueError('Unexpected unphased haplotype for '
                             f'{nanind} on position {nanrow.pos}')

    if dataframe is not None:
        return pd.concat([dataframe, new_frame], sort=False)

    if len(new_frame) == 0:
        return None

    return new_frame


def import_archaic_vcf(vcf_reader: TextIO,
                       dataframe: pd.DataFrame = None,
                       include_canc: bool = False) -> pd.DataFrame:
    '''
    Read in all lines of the provided, open vcf file and return a
    pandas dataframe.
    Expects a single individual, does not check phasing, and returns the
    number of alt sites only.  E.g. 0/1 -> 1, ./. -> 0, 1|1 -> 2
    '''
    usecols = ['chrom', 'pos', 'ref', 'alt', 'variant']
    if include_canc:
        usecols += ['infor']
    header = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual',
              'filter', 'infor', 'format', 'variant']
    result = pd.read_csv(vcf_reader,
                         delimiter='\t',
                         header=None,
                         names=header,
                         usecols=usecols,
                         comment='#')

    result = result.loc[(result.ref.str.len() == 1)
                        & (result.alt.str.len() == 1)]
    if include_canc:
        result = result.loc[result.infor.str.contains('CAnc')]

        # extract CAnc into separate column, drop info
        def extract_canc(info):
            for token in info.split(';'):
                tokens = token.split('=')
                if tokens[0] == 'CAnc':
                    return tokens[1]
        result.insert(len(result.columns),
                      "CAnc",
                      result.infor.map(extract_canc))
        result = result.drop(columns='infor')

    phase_decoder = {}
    for i in range(2):
        for j in range(2):
            phase_decoder[f'{i}|{j}'] = i + j
            phase_decoder[f'{i}/{j}'] = i + j
    phase_decoder['./.'] = 0

    if include_canc:
        phase_decoder['./.'] = -1  # to filter later

    def process_phase(x):
        return phase_decoder[x[:3]]

    result.variant = result.variant.map(process_phase)

    if include_canc:
        result = result.loc[result.variant != -1]

    if dataframe is not None:
        return pd.concat([dataframe, result], sort=False)

    return result
