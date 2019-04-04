import thin_vcf as main
import numpy as np
from numpy.testing import assert_array_equal as aae
from io import StringIO


def arg_helper(args, nondefault={}):
    defaults = {
        'bed_file': None,
        'merge': None,
        'individuals_file': None,
    }
    for k, v in nondefault.items():
        defaults[k] = v

    for k, v in args.items():
        assert v == defaults[k]


def test_read_args():
    # test defaults
    args = main.read_args([])
    arg_helper(args.__dict__)

    # set each
    args = main.read_args('--bed_file test'.split())
    arg_helper(args.__dict__, {'bed_file': 'test'})

    args = main.read_args('--merge 1'.split())
    arg_helper(args.__dict__, {'merge': 1})

    args = main.read_args('--individuals_file test'.split())
    arg_helper(args.__dict__, {'individuals_file': 'test'})


def test_read_bed_basic():
    bed = StringIO(
        '1\t100\t105\n'
    )
    output = main.read_bed(bed)
    aae(output, np.array([100, 105]))

    bed = StringIO(
        '1\t100\t105\n'
        '1\t104\t115\n'
    )
    output = main.read_bed(bed)
    aae(output, np.array([[100, 105], [104, 115]]))

    bed = StringIO(
        '1\t100\t105\n'
        '1\t104\t115\n'
        '1\t116\t120\n'
        '1\t125\t140\n'
        '1\t160\t180\n'
    )
    output = main.read_bed(bed)
    aae(output, np.array([[100, 105],
                          [104, 115],
                          [116, 120],
                          [125, 140],
                          [160, 180],
                          ]))

def test_read_bed_simple_merge():
    bed = StringIO(
        '1\t100\t105\n'
        '1\t104\t115\n'
        '1\t116\t120\n'
        '1\t125\t140\n'
        '1\t160\t180\n'
    )
    output = main.read_bed(bed, 0)
    aae(output, np.array([[100, 115],
                          [116, 120],
                          [125, 140],
                          [160, 180],
                          ]))

    bed = StringIO(
        '1\t100\t105\n'
        '1\t104\t115\n'
        '1\t116\t120\n'
        '1\t125\t140\n'
        '1\t160\t180\n'
    )
    output = main.read_bed(bed, 1)
    aae(output, np.array([[100, 120],
                          [125, 140],
                          [160, 180],
                          ]))

    bed = StringIO(
        '1\t100\t105\n'
        '1\t104\t115\n'
        '1\t116\t120\n'
        '1\t125\t140\n'
        '1\t160\t180\n'
    )
    output = main.read_bed(bed, 5)
    aae(output, np.array([[100, 140],
                          [160, 180],
                          ]))


def test_read_bed_more_merge():
    # single line
    bed = StringIO(
        '1\t100\t105\n'
    )
    output = main.read_bed(bed, 0)
    aae(output, np.array([100, 105]))

    # no merge
    bed = StringIO(
        '1\t100\t105\n'
        '1\t110\t115\n'
        '1\t116\t120\n'
    )
    output = main.read_bed(bed, 0)
    aae(output, np.array([[100, 105],
                          [110, 115],
                          [116, 120],
                          ]))

    # all merge
    bed = StringIO(
        '1\t100\t105\n'
        '1\t104\t115\n'
        '1\t114\t120\n'
    )
    output = main.read_bed(bed, 0)
    aae(output, np.array([[100, 120],
                          ]))

    # merge first
    bed = StringIO(
        '1\t100\t105\n'
        '1\t104\t115\n'
        '1\t116\t120\n'
    )
    output = main.read_bed(bed, 0)
    aae(output, np.array([[100, 115],
                          [116, 120],
                          ]))

    # merge last
    bed = StringIO(
        '1\t100\t105\n'
        '1\t106\t115\n'
        '1\t114\t120\n'
    )
    output = main.read_bed(bed, 0)
    aae(output, np.array([[100, 105],
                          [106, 120],
                          ]))


def test_read_indivs():
    indivs = StringIO(
        ''
    )
    assert main.read_indivs(indivs) == []

    indivs = StringIO(
        'I1'
    )
    assert main.read_indivs(indivs) == ['I1']

    indivs = StringIO(
        'I1\n'
        'I2  \n'
        'I3'
    )
    assert main.read_indivs(indivs) == 'I1 I2 I3'.split()


def test_match_indivs():
    header = 'bla bla i1 i2 i3 i4 i5 i6 i7 i8 i9 i1000'.split()
    indivs = ['i1']
    matches = main.match_indivs(header, indivs)
    aae(matches, [2])

    header = 'bla bla i1 i2 i3 i4 i5 i6 i7 i8 i9 i1000'.split()
    indivs = ['i5', 'i1']
    matches = main.match_indivs(header, indivs)
    aae(matches, [2, 6])

    header = 'bla bla i1 i2 i3 i4 i5 i6 i7 i8 i9 i1000'.split()
    indivs = ['i15']
    matches = main.match_indivs(header, indivs)
    aae(matches, [])

    header = 'bla bla i1 i2 i3 i4 i5 i6 i7 i8 i9 i1000'.split()
    indivs = ['i5', 'i7', 'i10']
    matches = main.match_indivs(header, indivs)
    aae(matches, [6, 8])


def test_line_parser_no_op():
    p = main.line_parser(None, None)
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('# stuff') == '# stuff'
    assert p.process_line('stuff') == 'stuff'
    assert p.process_line('stuff\n') == 'stuff\n'
    assert p.process_line('stuff\tstuff\n') == 'stuff\tstuff\n'


def test_line_parser_bed():
    p = main.line_parser(np.array([[3, 4], [5, 6]]), None)
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('# stuff') == '# stuff'
    assert p.process_line(
        '3 3 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''

    assert p.process_line(
        '3 4 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        '3 4 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')

    assert p.process_line(
        '3 5 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''

    assert p.process_line(
        '3 6 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        '3 6 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')

    assert p.process_line(
        '3 7 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''


def test_line_parser_bed_2():
    p = main.line_parser(np.array([[10, 20], [30, 50]]), None)
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('# stuff') == '# stuff'
    assert p.process_line(
        '3 1 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''
    assert p.process_line(
        '3 10 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''
    assert p.process_line(
        '3 15 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        '3 15 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')
    assert p.process_line(
        '3 20 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        '3 20 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')
    assert p.process_line(
        '3 21 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''
    assert p.process_line(
        '3 35 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        '3 35 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')
    assert p.process_line(
        '3 51 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''
    assert p.process_line(
        '3 61 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''


def test_line_parser_indiv():
    p = main.line_parser(None, ['u1', 'u3'])
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line(
        '#c p a b c d e f g u1 u2 u3\n'.replace(' ', '\t')) == \
        '#c p a b c d e f g u1 u3\n'.replace(' ', '\t')
    assert p.process_line(
        '3 1 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        '3 1 a b c d e f g ./. 1|1\n'.replace(' ', '\t')

    p = main.line_parser(None, ['u1'])
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line(
        '#c p a b c d e f g u1 u2 u3\n'.replace(' ', '\t')) == \
        '#c p a b c d e f g u1\n'.replace(' ', '\t')
    assert p.process_line(
        '3 1 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        '3 1 a b c d e f g ./.\n'.replace(' ', '\t')


def test_line_parser_both():
    p = main.line_parser(np.array([[10, 20], [30, 50]]), ['u1', 'u3'])
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line('## stuff') == '## stuff'
    assert p.process_line(
        '#c p a b c d e f g u1 u2 u3\n'.replace(' ', '\t')) == \
        '#c p a b c d e f g u1 u3\n'.replace(' ', '\t')
    assert p.process_line(
        '3 11 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        '3 11 a b c d e f g ./. 1|1\n'.replace(' ', '\t')
    assert p.process_line(
        '3 61 a b c d e f g ./. 0|0 1|1\n'.replace(' ', '\t')) == \
        ''
