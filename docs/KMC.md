# KMC Overview

KMC is a kmer counter, which can be found [here](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download).
In addition to being able to count kmers, version 3.0 has some neat extra utilities, such as being able to find intersections between two sets of kmers.


NOTES: 

- The kmc executables must be accessible from your $PATH for these wrappers to function.
- Any extra options you want to specify for KMC commands can be specified using `parameter='argument'`.
- Wrappers will return STDOUT and STDERR.
- All wrappers will automatically look for reverse reads if they're present, so you can be lazy and only specify your forward reads if the reads are paired. For this to work, forward reads and reverse reads must be in the same folder and follow Illumina's R1/R2 naming convention. If you have paired reads that don't follow these assumptions, use the keyword argument `reverse_in='path/to/paired_reads'`.

### KMC (Count Kmers)

```python
from biotools import kmc
out, err = kmc.kmc(forward_in, database_name, min_occurrences=1, reverse_in='NA', k=31, cleanup=True, tmpdir='tmp')
```

This command will count kmers of size 31 (changeable with `k`) and place the output in `database_name`. Low frequency kmers can be screened out by changing `min_occurrences`. By default, this command will create a temporary directory needed by KMC called `tmp`, which will be deleted upon completion of kmer counting. If you want to keep the temporary directory, change the `cleanup` argument to `False`, and if you would prefer a different temporary directory, specify it with the `tmpdir` parameter.

IMPORTANT NOTE: KMC assumes that input will be in FASTQ format. If your input file is in FASTA format, add `fm=''` to your function call.

### KMC Dump

```python
from biotools import kmc
out, err = kmc.dump(database, output, min_occurences=1, max_occurences=250)
```

This command will dump the kmers contained in `database` to a human-readable format, ignoring kmers with fewer instances than `min_occurrences` and kmers with more instances than `max_occurrences`.


### KMC Intersect

```python
from biotools import kmc
out, err = kmc.intersect(database_1, database_2, results)
```

This command will create a database (`results`) that contains kmers that are present in both `database_1` and `database_2`, where both these databases were created using `kmc.kmc`. The `results` database can then be dumped using `kmc.dump` in order to inspect it.


### KMC Union

```python
from biotools import kmc
out, err = kmc.union(database_1, database_2, results)
```

This command will create a database (`results`) that contains all kmers that are present in either `database_1` or `database_2`, where both these databases were created using `kmc.kmc`. The `results` database can then be dumped using `kmc.dump` in order to inspect it.


### KMC Subtract

```python
from biotools import kmc
out, err = kmc.subtract(database_1, database_2, results, exclude_below=1)
```

This command will create a database (`results`) that contains kmers that are present in `database_1` but not in `database_2`, 
where both these databases were created using `kmc.kmc`. The `results` database can then be dumped using `kmc.dump` in order to inspect it.

By default, all kmers in database\_2 will be subtracted from database\_1. This can be changed with the exclude\_below argument. Setting it to 3, for example, will only subtract kmers that have counts of 3 or more - kmers with counts of 1 or 2 will not be subtracted.
