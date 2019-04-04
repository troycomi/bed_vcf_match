from bed_vcf_match import read_vcf
from io import StringIO
import numpy as np
import pytest


def test_import_vcf():
    # add individual to target
    vcf = StringIO(
        '##comment\n'
        '##comment\n'
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV1\tUV2\n'
        '8\t10346\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n'
        '3\t1036\t.\tC\tG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
        '1\t1336\t.\tG\tT\t.\tPASS\t.\tGT\t1/0\t1|0\n'
    )

    with pytest.raises(ValueError) as e:
        df = read_vcf.import_vcf(vcf, individuals=['UV3'])
    assert 'UV3 not in file!' in str(e)

    vcf = StringIO(
        '##comment\n'
        '##comment\n'
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV1\tUV2\n'
        '8\t10346\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n'
        '3\t1036\t.\tC\tG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
        '1\t1336\t.\tG\tT\t.\tPASS\t.\tGT\t1/0\t1|0\n'
    )
    df = read_vcf.import_vcf(vcf, individuals=['UV2'])
    assert list(df.columns.values) ==\
        ['chrom', 'pos', 'ref', 'alt', 'UV2']
    assert list(df['chrom']) ==\
        [8, 3, 1]
    assert list(df['pos']) ==\
        [10346, 1036, 1336]
    assert list(df['ref']) ==\
        'A C G'.split()
    assert list(df['alt']) ==\
        'G G T'.split()
    assert list(df['UV2']) ==\
        ['0|0', np.nan, '1|0']

    vcf = StringIO(
        '##comment\n'
        '##comment\n'
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV1\tUV2\n'
        '8\t10346\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n'
        '3\t1036\t.\tC\tG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
        '1\t1336\t.\tG\tT\t.\tPASS\t.\tGT\t1/0\t1|0\n'
    )
    df = read_vcf.import_vcf(vcf, individuals=['UV2', 'UV1'])
    assert list(df.columns.values) ==\
        ['chrom', 'pos', 'ref', 'alt', 'UV1', 'UV2']
    assert list(df['chrom']) ==\
        [8, 3, 1]
    assert list(df['pos']) ==\
        [10346, 1036, 1336]
    assert list(df['ref']) ==\
        'A C G'.split()
    assert list(df['alt']) ==\
        'G G T'.split()
    assert list(df['UV1']) ==\
        ['0|1', '0|0', np.nan]
    assert list(df['UV2']) ==\
        ['0|0', np.nan, '1|0']

    # with new dataframe
    vcf = StringIO(
        '##comment\n'
        '##comment\n'
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV1\tUV2\n'
        '8\t10346\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n'
        '3\t1036\t.\tC\tG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
        '1\t1336\t.\tG\tT\t.\tPASS\t.\tGT\t1/0\t1|0\n'
    )

    df = read_vcf.import_vcf(vcf)
    assert list(df.columns.values) ==\
        ['chrom', 'pos', 'ref', 'alt', 'UV1', 'UV2']
    assert list(df['chrom']) ==\
        [8, 3, 1]
    assert list(df['pos']) ==\
        [10346, 1036, 1336]
    assert list(df['ref']) ==\
        'A C G'.split()
    assert list(df['alt']) ==\
        'G G T'.split()
    assert list(df['UV1']) ==\
        ['0|1', '0|0', np.nan]
    assert list(df['UV2']) ==\
        ['0|0', np.nan, '1|0']

    # with existing dataframe
    vcf = StringIO(
        '##comment\n'
        '##comment\n'
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV3\tUV2\n'
        '10\t10346\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n'
        '1\t1036\t.\tC\tG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
    )

    df = read_vcf.import_vcf(vcf, df)
    assert list(df.columns.values) ==\
        ['chrom', 'pos', 'ref', 'alt', 'UV1', 'UV2', 'UV3']
    assert list(df['chrom']) ==\
        [8, 3, 1, 10, 1]
    assert list(df['pos']) ==\
        [10346, 1036, 1336, 10346, 1036]
    assert list(df['ref']) ==\
        'A C G A C'.split()
    assert list(df['alt']) ==\
        'G G T G G'.split()
    assert list(df['UV1']) ==\
        ['0|1', '0|0', np.nan, np.nan, np.nan]
    assert list(df['UV2']) ==\
        ['0|0', np.nan, '1|0', '0|0', np.nan]
    assert list(df['UV3']) ==\
        [np.nan, np.nan, np.nan, '0|1', '0|0']

    # skip multi allelic sites
    vcf = StringIO(
        '##comment\n'
        '##comment\n'
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV3\tUV2\n'
        '10\t10346\t.\tAA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n'
        '1\t1036\t.\tC\tGG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
        '1\t1036\t.\tC\tG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
    )

    df = read_vcf.import_vcf(vcf)
    assert list(df.columns.values) ==\
        ['chrom', 'pos', 'ref', 'alt', 'UV3', 'UV2']
    assert list(df['chrom']) ==\
        [1]
    assert list(df['pos']) ==\
        [1036]
    assert list(df['ref']) ==\
        ['C']
    assert list(df['alt']) ==\
        ['G']
    assert list(df['UV3']) ==\
        ['0|0']
    assert np.isnan(df['UV2']).all()

    # all multi-allelic
    vcf = StringIO(
        '##comment\n'
        '##comment\n'
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV3\tUV2\n'
        '10\t10346\t.\tAA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n'
        '1\t1036\t.\tC\tGG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
    )

    df = read_vcf.import_vcf(vcf)
    assert df is None

    # check phasing
    vcf = StringIO(
        '##comment\n'
        '##comment\n'
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV3\tUV2\n'
        '10\t10346\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n'
        '1\t1036\t.\tC\tG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
    )

    with pytest.raises(ValueError) as e:
        df = read_vcf.import_vcf(vcf, df, check_phasing=True)
    assert 'Unexpected unphased haplotype for UV2 on position 1036' in str(e)

    vcf = StringIO(
        '#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfor\tformat\tUV3\tUV2\n'
        '10\t10346\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0|0\n'
        '1\t1036\t.\tC\tG\t.\tPASS\t.\tGT\t0|0\t1/0\n'
    )

    with pytest.raises(ValueError) as e:
        df = read_vcf.import_vcf(vcf, df, check_phasing=True)
    assert 'Unexpected unphased haplotype for UV3 on position 10346' in str(e)


def test_import_archaic_vcf():
    # actual line
    vcf = StringIO('14\t19073582\t.\tG\tA\t3576.58\t.\t'
                   'AC=2;AF=1.00;AN=2;BaseQRankSum=-1.395;DP=121;Dels=0.00;'
                   'FS=0.000;HRun=0;HaplotypeScore=1.9955;MQ=30.12;MQ0=8;'
                   'MQRankSum=-0.322;QD=29.56;ReadPosRankSum=0.904\t'
                   'GT:DP:GQ:PL:A:C:G:T:IR\t'
                   '1/1:121:99:3577,296,0:79,40:0,0:1,0:0,0:0')
    df = read_vcf.import_archaic_vcf(vcf)
    assert list(df['variant']) ==\
        [2]
    assert list(df['chrom']) ==\
        [14]
    assert list(df['pos']) ==\
        [19073582]
    assert list(df['ref']) ==\
        ['G']
    assert list(df['alt']) ==\
        ['A']

    # more lines, simplify extra stuff
    vcf = StringIO(
        '14\t19073582\t.\tG\tA\t.\t.\t.\t.\t1/1:\n'
        '1\t1073582\t.\tA\tT\t.\t.\t.\t.\t0|1:\n'
        '4\t1973582\t.\tT\tC\t.\t.\t.\t.\t1/0:\n'
        '2\t1903582\t.\tC\tG\t.\t.\t.\t.\t0/0:\n'
        '3\t1907582\t.\tG\t.\t.\t.\t.\t.\t./.:\n'
        '3\t1907582\t.\tGG\tA\t.\t.\t.\t.\t./.:\n'
    )
    df = read_vcf.import_archaic_vcf(vcf)
    assert list(df['variant']) ==\
        [2, 1, 1, 0, 0]
    assert list(df['chrom']) ==\
        [14, 1, 4, 2, 3]
    assert list(df['pos']) ==\
        [19073582, 1073582, 1973582, 1903582, 1907582]
    assert list(df['ref']) ==\
        'G A T C G'.split()
    assert list(df['alt']) ==\
        'A T C G .'.split()
