from bed_vcf_match import analyze_bed
from io import StringIO
import pandas as pd
from pandas.util.testing import assert_frame_equal


def test_read_bed():
    bed_file = ('1\t2860418\t2909250\n'
                '2\t2913234\t2957133\n'
                '3\t3382004\t3460679\n'
                '1\t5017441\t5094205\n')
    bed_input = StringIO(bed_file)

    expected = [[int(t) for t in l.split()] for l in bed_file.split('\n')]

    for i, line in enumerate(analyze_bed.read_bed(bed_input)):
        assert list(line.values[0]) == expected[i]


def test_filter_modern_db():
    modern = StringIO(
        'chrom,pos,ref,alt,UV1,UV2\n'
        '1,100,A,T,0|1,0|0\n'
        '2,200,C,G,1|1,nan\n'
        '3,300,T,A,nan,0|1\n'
        '3,305,T,A,nan,0|0\n'
        '1,105,A,T,nan,1|0\n'
        '1,110,C,G,nan,1|1\n'
    )
    database = pd.read_csv(modern)
    # select nans
    rows = analyze_bed.filter_modern_db(database, 1, 100, 110, 1, 'UV1')
    assert len(rows) == 1

    # basic, test matches
    rows = analyze_bed.filter_modern_db(database, 1, 100, 110, 1, 'UV2')
    assert len(rows) == 3
    assert list(rows['ref']) == 'A A C'.split()
    assert list(rows['alt']) == 'T T G'.split()
    assert list(rows['variant']) == [0, 1, 1]
    assert list(rows['chrom']) == [1] * 3
    assert list(rows['pos']) == [100, 105, 110]

    # check <= for end then start
    rows = analyze_bed.filter_modern_db(database, 1, 100, 109, 1, 'UV2')
    assert len(rows) == 2
    assert list(rows['pos']) == [100, 105]
    rows = analyze_bed.filter_modern_db(database, 1, 101, 110, 1, 'UV2')
    assert len(rows) == 2
    assert list(rows['pos']) == [105, 110]

    # when start == end
    rows = analyze_bed.filter_modern_db(database, 1, 100, 100, 2, 'UV2')
    assert len(rows) == 1
    assert list(rows['pos']) == [100]

    # another individual
    rows = analyze_bed.filter_modern_db(database, 2, 200, 210, 1, 'UV1')
    assert list(rows['chrom']) == [2]
    assert list(rows['pos']) == [200]
    assert len(rows) == 1

    # non existant chrom
    rows = analyze_bed.filter_modern_db(database, 8, 100, 110, 2, 'UV1')
    assert len(rows) == 0


def test_join_vcfs():
    modern = StringIO(
        'chrom,pos,individual,haplotype,variant,ref,alt\n'
        '1,100,UV1,1,0,A,T\n'
        '1,100,UV1,2,1,A,T\n'
        '1,100,UV2,1,0,A,T\n'
        '1,100,UV2,2,0,A,T\n'
        '2,200,UV1,1,1,C,G\n'
        '2,200,UV1,2,1,C,G\n'
        '3,300,UV2,1,0,T,A\n'
        '3,300,UV2,2,1,T,A\n'
        '3,305,UV2,1,0,G,C\n'
        '3,305,UV2,2,0,G,C\n'
        '1,105,UV2,1,1,A,T\n'
        '1,105,UV2,2,0,A,T\n'
        '1,110,UV2,1,1,C,G\n'
        '1,110,UV2,2,1,C,G\n'
    )
    modern = pd.read_csv(modern)

    archaic = StringIO(
        'chrom,pos,ref,alt,variant\n'
        '1,100,A,T,1\n'  # match with 1
        '2,200,C,G,2\n'  # match with 2
        '3,300,T,G,2\n'  # match except alt
        '3,305,T,C,2\n'  # match except ref
        '1,105,A,T,0\n'  # match with 0
    )  # 1,110 is omitted
    archaic = pd.read_csv(archaic)

    expected = StringIO(
        'chrom,pos,individual,haplotype,variant,ref,alt,archaic\n'
        '1,100,UV1,1,0,A,T,1\n'
        '1,100,UV1,2,1,A,T,1\n'
        '1,100,UV2,1,0,A,T,1\n'
        '1,100,UV2,2,0,A,T,1\n'
        '2,200,UV1,1,1,C,G,2\n'
        '2,200,UV1,2,1,C,G,2\n'
        '3,300,UV2,1,0,T,A,0\n'
        '3,300,UV2,2,1,T,A,0\n'
        '3,305,UV2,1,0,G,C,0\n'
        '3,305,UV2,2,0,G,C,0\n'
        '1,105,UV2,1,1,A,T,0\n'
        '1,105,UV2,2,0,A,T,0\n'
        '1,110,UV2,1,1,C,G,0\n'
        '1,110,UV2,2,1,C,G,0\n'
    )
    expected = pd.read_csv(expected)

    joined = analyze_bed.join_vcf(modern, archaic)
    assert_frame_equal(joined, expected)


def test_summarize_region():
    modern = StringIO(
        'chrom,pos,ref,alt,UV1,UV2\n'
        '1,100,A,T,0|1,0|0\n'
        '1,105,A,T,nan,1|0\n'
        '1,110,C,G,nan,1|1\n'
        '1,115,A,T,nan,1|0\n'
        '1,120,C,G,nan,1|1\n'
    )
    modern = pd.read_csv(modern)

    # no archaic, just match number of sites and variants
    summary = analyze_bed.summarize_region([1, 100, 120], 1, 'UV2', modern)
    assert summary == '1\t100\t120\t5\t4\n'
    summary = analyze_bed.summarize_region([1, 100, 120], 2, 'UV2', modern)
    assert summary == '1\t100\t120\t5\t2\n'
    # no matches
    summary = analyze_bed.summarize_region([2, 100, 120], 2, 'UV2', modern)
    assert summary == '2\t100\t120\t0\t0\n'

    # one archaic, no matches
    archaic1 = StringIO(
        'chrom,pos,ref,alt,variant\n'
        '1,100,A,G,1\n'
        '1,105,A,C,2\n'
    )
    archaic1 = pd.read_csv(archaic1)
    summary = analyze_bed.summarize_region([1, 100, 120], 1, 'UV2',
                                           modern, archaic1)
    assert summary == '1\t100\t120\t5\t4\t0.0\t0.0\t0.0\n'
    # no matches
    summary = analyze_bed.summarize_region([2, 100, 120], 1, 'UV2',
                                           modern, archaic1)
    assert summary == '2\t100\t120\t0\t0\t0.0\t0.0\tnan\n'

    # another archaic, with matches
    archaic2 = StringIO(
        'chrom,pos,ref,alt,variant\n'
        '1,100,A,T,1\n'
        '1,105,A,T,2\n'
        '1,110,C,G,0\n'
    )
    archaic2 = pd.read_csv(archaic2)
    summary = analyze_bed.summarize_region([1, 100, 120], 1, 'UV2',
                                           modern, archaic1, archaic2)
    assert summary == '1\t100\t120\t5\t4\t0.0\t0.0\t0.0\t1.5\t1.0\t0.2\n'
