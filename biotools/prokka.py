from biotools import accessoryfunctions


def kwargs_to_string(kwargs):
    """
    Given a set of kwargs, turns them into a string which can then be passed to a command.
    :param kwargs: kwargs from a function call.
    :return: outstr: A string, which is '' if no kwargs were given, and the kwargs in string format otherwise.
    """
    outstr = ''
    for arg in kwargs:
        outstr += ' --{} {}'.format(arg, kwargs[arg])
    return outstr


def prokka(input_fasta, output_dir, output_name, **kwargs):
    options = kwargs_to_string(kwargs)
    cmd = 'prokka --outdir {} --prefix {} {} {}'.format(output_dir, output_name, options, input_fasta)
    out, err = accessoryfunctions.run_subprocess(cmd)
    return out, err
