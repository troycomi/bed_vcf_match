'''
read_vcf

Read a vcf file into pandas data frame
'''


import pandas as pd
from typing import TextIO, List


# contants for column numbers
CHROMOSOME = 0
POSITION = 1
REFERENCE = 3
ALTERNATIVE = 4
HAPLOTYPES = 9


def import_vcf(vcf_reader: TextIO,
               dataframe: pd.DataFrame = None,
               check_phasing: bool = False) -> pd.DataFrame:
    '''
    Read in all lines of the provided, open vcf file and merge into provided
    pandas dataframe
    '''
    header_lines = 1  # 1-based indexing on error reporting
    for line in vcf_reader:
        header_lines += 1
        if line[1] == '#':  # comment string
            continue

        if line[0] == '#':  # header string
            indivs = line.rstrip().split('\t')[9:]
            break

    print(f'read in {header_lines} of header')
    print(f'finished header with {len(indivs)} indivs', flush=True)

    try:
        if check_phasing:
            new_frame = pd.concat([parse_line(line, indivs, i + header_lines)
                                   for i, line in enumerate(vcf_reader)])
        else:
            new_frame = pd.concat([parse_line(line, indivs)
                                   for line in vcf_reader])

    except ValueError as e:
        if str(e) != "All objects passed were None":
            raise e
        else:
            new_frame = None

    if dataframe is not None:
        return pd.concat([dataframe, new_frame])

    return new_frame


def parse_line(line: str,
               individuals: List[str],
               all_phased: int = False) -> pd.DataFrame:
    '''
    read in line of vcf, returning a datafram with columns
    chrom, pos, ref, alt, individual, haplotype, variant
    all_phased will check that any unphased sites are ./.
    raising an exception if not
    '''
    tokens = line.rstrip().split('\t')
    chrom = int(tokens[CHROMOSOME])
    pos = int(tokens[POSITION])
    ref = tokens[REFERENCE]
    alt = tokens[ALTERNATIVE]
    haplotypes = tokens[HAPLOTYPES:]
    entries = []
    if all_phased and all_phased % 100 == 0:
        print(all_phased, flush=True)

    if len(ref) + len(alt) > 2:
        return None
    for i, entry in enumerate(haplotypes):
        if entry[1] == '|':  # only consider phased sites
            entries.append((individuals[i], 1, int(entry[0])))
            entries.append((individuals[i], 2, int(entry[2])))
        elif entry == './.':  # unless it's this
            entries.append((individuals[i], 1, 0))
            entries.append((individuals[i], 2, 0))
        elif all_phased is not False:
            raise ValueError(f'Unexpected unphased haplotype {entry} '
                             f'for {individuals[i]} on line {all_phased}!')

    result = pd.DataFrame.from_records(entries,
                                       columns=['individual', 'haplotype',
                                                'variant'])
    return result.assign(chrom=chrom,
                         pos=pos,
                         ref=ref,
                         alt=alt)


def import_archaic_vcf(vcf_reader: TextIO) -> pd.DataFrame:
    '''
    Read in all lines of the provided, open vcf file and return a
    pandas dataframe.
    Expects a single individual, does not check phasing, and returns the
    number of alt sites only.  E.g. 0/1 -> 1, ./. -> 0, 1|1 -> 2
    '''
    entries = []
    for line in vcf_reader:
        tokens = line.split('\t')

        chrom = int(tokens[CHROMOSOME])
        pos = int(tokens[POSITION])
        ref = tokens[REFERENCE]
        alt = tokens[ALTERNATIVE]

        if len(ref) + len(alt) > 2:
            continue

        genotype = tokens[HAPLOTYPES][:3].replace('.', '0')
        variant = int(genotype[0]) + int(genotype[2])

        entries.append(
            (chrom,
             pos,
             ref,
             alt,
             variant)  # variant
        )
    return pd.DataFrame.from_records(entries,
                                     columns=['chrom', 'pos', 'ref',
                                              'alt', 'variant'])
