# Fasta file manipulation.


def clean_names(in_fasta, out_fasta='NA', method='split', delimiter='NA', truncate_value=10):
    """
    Takes a fasta file, and cleans up headers that are overly long. Default behaviour is to modify input fasta in place,
    splitting title on whitespace. Can also truncate to a certain number of characters.
    :param in_fasta: Path to input fasta file. Must be uncompressed.
    :param out_fasta: Path to output fasta file. If not specified, input is modified in place.
    :param method: Option between split on delimiter of your choice and truncate to number of characters.
    :param delimiter: Delimiter to split on. Defaults to splitting on whitespace, but can be changed to anything.
    :param truncate_value: Number of characters to truncate to if using method='truncate'
    :return:
    """
    with open(in_fasta) as fasta_file:
        lines = fasta_file.readlines()
    if out_fasta == 'NA':
        out_file = in_fasta
    else:
        out_file = out_fasta
    with open(out_file, 'w') as f:
        for line in lines:
            if '>' in line:
                if method == 'split':
                    if delimiter == 'NA':
                        f.write(line.split()[0] + '\n')
                    else:
                        f.write(line.split(delimiter)[0] + '\n')
                elif method == 'truncate':
                    if len(line) > truncate_value:
                        f.write(line[0:truncate_value] + '\n')
                    else:
                        f.write(line)
                else:
                    raise ValueError('Valid values of method are split or truncate. Please specify one of those options.')
            else:
                f.write(line)
