# Jellyfish Overview

Jellyfish is a popular kmer counter, which can be found [here](https://github.com/gmarcais/Jellyfish/releases).


NOTES: 

- The jellyfish executable must be accessible from your $PATH for these wrappers to function.
- Any options you want to specify for Jellyfish commands that are not implemented as keyword arguments should be entered as options='', where the options string contains the options exactly as they would be specified on the command line.
- Wrappers will return STDOUT and STDERR.
- All wrappers will automatically look for reverse reads if they're present, so you can be lazy and only specify your forward reads if the reads are paired. For this to work, forward reads and reverse reads must be in the same folder and follow Illumina's R1/R2 naming convention. If you have paired reads that don't follow these assumptions, use the keyword argument `reverse_in='path/to/paired_reads'`.

### Jellyfish Count

```python
from biotools import jellyfish
out, err = jellyfish.count(forward_in, reverse_in='NA', kmer_size=31, count_file='mer_counts.jf', hash_size='100M')
```

This command will count kmers of size 31 (changeable with `kmer_size`) and place the output in `mer_counts.jf`, which can be changed with the option `count_file`. Input reads can be uncompressed or gzip compressed, and as of Jellyfish 2.0 `hash_size` doesn't really matter - it will be increased automatically if the initial specification isn't enough.

### Jellyfish Dump

```python
from biotools import jellyfish
out, err = jellyfish.dump(mer_file, output_file='counts.fasta')
```

This command will dump the counts created by `jellyfish count` into a human-readable format, default FASTA.
