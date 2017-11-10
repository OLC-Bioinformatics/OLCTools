import gzip
import subprocess


def run_subprocess(command):
    """
    command is the command to run, as a string.
    runs a subprocess, returns stdout and stderr from the subprocess as strings.
    """
    x = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = x.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if x.returncode != 0:
        print('STDERR from called program: {}'.format(err))
        print('STDOUT from called program: {}'.format(out))
        raise subprocess.CalledProcessError(x.returncode, command)
    return out, err


def uncompress_gzip(infile, outfile='NA'):
    if not infile.endswith('.gz'):
        raise TypeError('Input file does not appear to be gzipped! Gzipped files should end with .gz!')
    if outfile == 'NA':
        outfile = infile.replace('.gz', '')
    with gzip.open(infile, 'rb') as f:
        file_content = f.read()
    with open(outfile, 'wb') as f:
        f.write(file_content)
    return outfile


def file_len(fname):
    i = -1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
