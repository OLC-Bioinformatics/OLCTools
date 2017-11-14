# BBTools Overview

BBTools is an excellent suite of programs created by Brian Bushnell, which you can find [here](https://jgi.doe.gov/data-and-tools/bbtools/).

NOTES: 

- The shell scripts used for calling the BBTools programs must be accessible from your $PATH for these wrappers to function.
- All wrappers return out and err, the STDOUT and STDERR from the program being called.
- All wrappers will automatically look for reverse reads if they're present, so you can be lazy and only specify your forward reads if the reads are paired. For this to work, forward reads and reverse
reads must be in the same folder and follow Illumina's R1/R2 naming convention. If you have paired reads that don't follow these assumptions, use the keyword argument `reverse_in='path/to/paired_reads'`.
- The default setting for all programs is to use the number of cores on your computer, as BBTools tend to be able to take advantage of multiple processors fairly well. To change this, add the argument `threads='8'` for 8 cores, or change it to whatever else you want to use.
- Default behaviour for BBTools is to not overwrite output files if they're already there. If they are there, BBTools will crash. To change this behaviour, set `overwrite='t'`.
- Any other parameters you want to change for BBTools are also possible to change, using the same `parameter='argument'` format.
- Wrappers can also have the command string returned. To do this, add the argument `returncmd=True`. This will cause the 
command string to be returned as the third return value.

The following wrappers for BBTools have been written:

### BBMap

```python
from biotools import bbtools
out, err = bbtools.bbmap(reference, forward_in, out_bam, returncmd=False, reverse_in='NA')
```

Here, reference should be the path to a FASTA (or multi-FASTA) file for reads to be aligned against, forward\_in is a the path to a set of FASTQ reads (either compressed or uncompressed), and out\_bam is
an output file to write the alignment to (use a .bam ending to get BAM output, .sam to get SAM output).


### BBDuk_Trim

```python
from biotools import bbtools
out, err = bbtools.bbduk_trim(forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA')
```

Does quality trimming of reads. Default settings (those used in the OLC Assembly Pipeline) are to trim bases below quality 20, have a minimum read length of 50, and remove any adapter sequences.


### BBDuk_Filter

```python
from biotools import bbtools
out, err = bbtools.bbduk_filter(reference, forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA')
```

Filters out any reads that match to reference from the input reads passed with forward\_in, and writes the clean reads to forward\_out.

### BBDuk_Bait

```python
from biotools import bbtools
out, err = bbtools.bbduk_bait(reference, forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA')
```

Baits out any reads that match to reference from the input reads passed with forward\_in, and writes them to forward\_out. 

### BBMerge

```python
from biotools import bbtools
out, err = bbtools.bbduk_bbmerge(forward_in, merged_reads, returncmd=False, reverse_in='NA')
```

BBMerge will join forward and reverse reads that have overlapping regions due to small insert size and write these reads to the location specified with merged\_reads.

### Dedupe

```python
from biotools import bbtools
out, err = bbtools.dedupe(input_file, output_file, returncmd=False)
```

Given an input file, Dedupe will remove any duplicate sequences and write the un-duplicated sequences to output\_file.

### BBNorm

```python
from biotools import bbtools
out, err = bbtools.bbnorm(forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA')
```

BBNorm will do read normalization to ensure that read depths don't get to be too ridiculous, which can cause issues for _de novo_ assemblies. The default kmer depth that BBNorm targets is 100, which can 
be changed using target='depth'.


### Tadpole

```python
from biotools import bbtools
out, err = bbtools.tadpole(forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA', mode='correct')
```

Tadpole is the error-corrector and/or assembler of BBTools. The default mode here is error correction, but other modes can be specified with mode='alternate_mode'. Available alternate modes are 'contig' (make contigs), 'extend' (extend sequences to be longer), and 'insert' (measures insert sizes).
