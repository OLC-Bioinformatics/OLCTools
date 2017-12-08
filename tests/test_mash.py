from biotools import mash
import pytest
import os


def test_sketch_no_input_files():
    with pytest.raises(ValueError):
        mash.sketch()


def test_dist_no_input_files():
    with pytest.raises(ValueError):
        mash.dist()


def test_mash_dist_call():
    out, err, cmd = mash.dist('tests/dummy_fastq/*fastq', output_file='tests/distances.tab', returncmd=True)
    assert cmd == 'mash dist tests/dummy_fastq/*fastq  -p 1  > tests/distances.tab'
    os.remove('tests/distances.tab')


def test_mash_dist_call_multithreaded():
    out, err, cmd = mash.dist('tests/dummy_fastq/*fastq', output_file='tests/distances.tab', returncmd=True, threads=4)
    assert cmd == 'mash dist tests/dummy_fastq/*fastq  -p 4  > tests/distances.tab'
    os.remove('tests/distances.tab')


def test_mash_dist_call_kwargs():
    out, err, cmd = mash.dist('tests/dummy_fastq/*fastq', output_file='tests/distances.tab', returncmd=True, s='34')
    assert cmd == 'mash dist tests/dummy_fastq/*fastq  -p 1  -s 34 > tests/distances.tab'
    os.remove('tests/distances.tab')


def test_mash_sketch_call():
    out, err, cmd = mash.sketch('tests/dummy_fastq/*fastq', output_sketch='tests/sketch.msh', returncmd=True)
    assert cmd == 'mash sketch tests/dummy_fastq/*fastq -o tests/sketch.msh -p 1 '
    os.remove('tests/sketch.msh')


def test_mash_sketch_call_multithreaded():
    out, err, cmd = mash.sketch('tests/dummy_fastq/*fastq', output_sketch='tests/sketch.msh', returncmd=True, threads=4)
    assert cmd == 'mash sketch tests/dummy_fastq/*fastq -o tests/sketch.msh -p 4 '
    os.remove('tests/sketch.msh')


def test_read_mash_dist():
    out, err, cmd = mash.dist('tests/dummy_fastq/*fastq', output_file='tests/distances.tab', returncmd=True)
    results = mash.read_mash_output('tests/distances.tab')
    assert results[1].reference == 'tests/dummy_fastq/single.fastq' \
        and results[1].query == 'tests/dummy_fastq/test_R2.fastq' \
        and results[1].distance == 0.00763536
    os.remove('tests/distances.tab')
