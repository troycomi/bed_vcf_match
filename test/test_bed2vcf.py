import bed2vcf as main


def arg_helper(args, nondefault={}):
    defaults = {
        'archaic_vcfs': [],
        'modern_vcfs': [],
        'bed_files': None,
        'output_dir': None,
        'vcf_output': False,
        'bed_output': False,
    }
    for k, v in nondefault.items():
        defaults[k] = v

    for k, v in args.items():
        assert v == defaults[k]


def test_read_args():
    # test defaults
    args = main.read_args([])
    arg_helper(args.__dict__)

    # set outputs
    args = main.read_args('--bed_output'.split())
    arg_helper(args.__dict__, {'bed_output': True})

    args = main.read_args('--bed_output --vcf_output'.split())
    arg_helper(args.__dict__,
               {'bed_output': True,
                'vcf_output': True})

    args = main.read_args('--vcf_output'.split())
    arg_helper(args.__dict__, {'vcf_output': True})

    # vcfs and bed files
    arg_name = ['bed_files', 'modern_vcfs', 'archaic_vcfs']
    arg_values = ['file1', 'file2 file3']

    # single calls
    for l in arg_name:
        for v in arg_values:
            args = main.read_args(f'--{l} {v}'.split())
            arg_helper(args.__dict__, {l: v.split()})

    # repeated calls
    for l in arg_name:
        args = main.read_args(
            f'--{l} {arg_values[0]} --{l} {arg_values[1]}'.split())
        arg_helper(args.__dict__, {l: ' '.join(arg_values).split()})
