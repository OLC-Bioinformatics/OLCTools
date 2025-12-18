# Mash Overview

Mash is a set of tools for creating and using MinHash sketches, a neat way of turning a genome into a small signature which can easily be compared with other signatures. 
Documentation on Mash can be found [here](http://mash.readthedocs.io/en/latest/), the Genome Biology paper describing Mash can be found [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x), and Mash can be downloaded [here](https://github.com/marbl/Mash/releases).

NOTES: 

- The Mash executable must be accessible from your $PATH for these wrappers to function.
- All wrappers return out and err, the STDOUT and STDERR from the program being called.
- Any number of arguments for sequence files can be passed to `mash.dist`, `mash.sketch`, and `mash.screen`. Unix wildcards (*, ?) can be part of these arguments.
- The default setting for all wrappers is to only use one thread. This can be changed with the `threads` keyword argument. For example, to use 8 threads, set `threads=8`.
- Any other parameters you want to change for Mash are also possible to change, using the same parameter='argument' format.
- Wrappers can also have the command string returned. To do this, add the argument `returncmd=True`. This will cause the command string to be returned as the third return value.

The following wrappers for Mash have been written:

### Mash Sketch

```python
from biotools import mash
out, err = mash.sketch(file_to_sketch_1, file_to_sketch_2, output_sketch='sketch.msh', threads=1)
```

This command will allow you to sketch one or many files, and save them to a sketch file for for quick analysis later on. The default value for the output sketch is `sketch.msh`, but this can be changed with the keyword argument `output_sketch`. Wildcards can be used to sketch lots of files at once. For example, to sketch all FASTA files in the directory `example`, one could call `mash.sketch('example/*.fasta')`.


### Mash Dist

```python
from biotools import mash
out, err = mash.dist(query_file_1, query_file_2, output_file='distances.tab', threads=1)
```

Mash dist will find distances between either FASTA files, FASTQ files, or previously sketched files. By default, these distances will be output to `distances.tab`, which can then be read by `mash.read_mash_output`. 


### Mash Screen

```python
from biotools import mash
out, err = mash.screen(query_file_1, query_file_2, output_file='screen.tab', threads=1)
```

Mash screen will find how well one query is contained within another. The first argument must be a sketch file, and subsequent arguments can be FASTA files, FASTQ files, or other sketches. This command will only work if mash>=2.0 is installed. Output of the mash screen command can be read by `mash.read_mash_screen`.


### Read Mash Output

```python
from biotools import mash
mash_results = mash.read_mash_output(result_file)
```

This command takes a result file from `mash.dist` as input, and will return a list of results, with each index corresponding to one line of the result file. Each item in the list has the following attributes: reference, query, distance, pvalue, and matching_hash.

### Read Mash Screen

```python
from biotools import mash
mash_results = mash.read_mash_screen(screen_result)
```

This command takes a result file from `mash.screen` as input, and will return a list of results, with each index corresponding to one line of the result file. Each item in the list has the following attributes: identity, shared\_hashes, median\_multiplicity, pvalue, and query_id.



