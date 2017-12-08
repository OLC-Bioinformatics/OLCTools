import pytest
import os
from biotools import bbtools


def test_bbmap_command_paired():
    out, err, cmd = bbtools.bbmap(reference='tests/dummy_fasta/test.fasta',
                                  forward_in='tests/dummy_fastq/test_R1.fastq',
                                  out_bam='tests/out.bam', returncmd=True)
    assert cmd == 'bbmap.sh ref=tests/dummy_fasta/test.fasta in=tests/dummy_fastq/test_R1.fastq' \
                  ' in2=tests/dummy_fastq/test_R2.fastq out=tests/out.bam nodisk'
    os.remove('tests/out.bam')


def test_bbmap_command_single():
    out, err, cmd = bbtools.bbmap(reference='tests/dummy_fasta/test.fasta',
                                  forward_in='tests/dummy_fastq/single.fastq',
                                  out_bam='tests/out.bam', returncmd=True)
    assert cmd == 'bbmap.sh ref=tests/dummy_fasta/test.fasta in=tests/dummy_fastq/single.fastq ' \
                  'out=tests/out.bam nodisk'
    os.remove('tests/out.bam')


def test_bbduk_bait_command_paired():
    out, err, cmd = bbtools.bbduk_bait(forward_in='tests/dummy_fastq/test_R1.fastq',
                                       forward_out='tests/out_R1.fastq',
                                       reference='tests/dummy_fasta/test.fasta',
                                       returncmd=True)
    assert cmd == 'bbduk.sh in=tests/dummy_fastq/test_R1.fastq ' \
                  'in2=tests/dummy_fastq/test_R2.fastq outm=tests/out_R1.fastq outm2=tests/out_R2.fastq' \
                  ' ref=tests/dummy_fasta/test.fasta'
    os.remove('tests/out_R1.fastq')
    os.remove('tests/out_R2.fastq')


def test_bbduk_bait_command_single():
    out, err, cmd = bbtools.bbduk_bait(forward_in='tests/dummy_fastq/single.fastq',
                                       forward_out='tests/out.fastq',
                                       reference='tests/dummy_fasta/test.fasta',
                                       returncmd=True)
    assert cmd == 'bbduk.sh in=tests/dummy_fastq/single.fastq ' \
                  'outm=tests/out.fastq' \
                  ' ref=tests/dummy_fasta/test.fasta'
    os.remove('tests/out.fastq')


def test_bbduk_filter_command_paired():
    out, err, cmd = bbtools.bbduk_filter(forward_in='tests/dummy_fastq/test_R1.fastq',
                                         forward_out='tests/out_R1.fastq',
                                         reference='tests/dummy_fasta/test.fasta',
                                         returncmd=True)
    assert cmd == 'bbduk.sh in=tests/dummy_fastq/test_R1.fastq ' \
                  'in2=tests/dummy_fastq/test_R2.fastq out=tests/out_R1.fastq out2=tests/out_R2.fastq' \
                  ' ref=tests/dummy_fasta/test.fasta'
    os.remove('tests/out_R1.fastq')
    os.remove('tests/out_R2.fastq')


def test_bbduk_filter_command_single():
    out, err, cmd = bbtools.bbduk_filter(forward_in='tests/dummy_fastq/single.fastq',
                                         forward_out='tests/out.fastq',
                                         reference='tests/dummy_fasta/test.fasta',
                                         returncmd=True)
    assert cmd == 'bbduk.sh in=tests/dummy_fastq/single.fastq ' \
                  'out=tests/out.fastq' \
                  ' ref=tests/dummy_fasta/test.fasta'
    os.remove('tests/out.fastq')


def test_bbmap_command_paired_kwargs():
    out, err, cmd = bbtools.bbmap(reference='tests/dummy_fasta/test.fasta',
                                  forward_in='tests/dummy_fastq/test_R1.fastq',
                                  out_bam='tests/out.bam', returncmd=True, threads='3', ordered='t')
    assert 'bbmap.sh ref=tests/dummy_fasta/test.fasta in=tests/dummy_fastq/test_R1.fastq' \
           ' in2=tests/dummy_fastq/test_R2.fastq out=tests/out.bam nodisk' in cmd and 'ordered=t' in cmd and 'threads=3' in cmd
    os.remove('tests/out.bam')


def test_bbduk_filter_exception():
    with pytest.raises(ValueError):
        out, err, cmd = bbtools.bbduk_filter(forward_in='tests/dummy_fastq/test_R1.fastq',
                                             forward_out='tests/out_aa.fastq',
                                             reference='tests/dummy_fasta/test.fasta',
                                             returncmd=True)


# Would have a test for bbnorm here, but it seems to like to use 50GB of memory, which will cause some
# major issues with Travis-CI. Will leave it out.

