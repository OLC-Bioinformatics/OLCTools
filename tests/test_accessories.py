from olctools.accessoryFunctions.accessoryFunctions import combinetargets, dependency_check, filer, \
    find_paired_reads, find_unpaired_reads, make_path
import pytest
import shutil
import glob
import os


def test_dependency_check_good():
    assert dependency_check('ls') is True


def test_dependency_check_bad():
    assert dependency_check('aaaaaaddsfsdss') is False


def test_paired_read_find():
    assert find_paired_reads('tests/dummy_fastq') == [['tests/dummy_fastq/test_R1.fastq',
                                                       'tests/dummy_fastq/test_R2.fastq']]


def test_unpaired_read_find():
    assert find_unpaired_reads('tests/dummy_fastq') == ['tests/dummy_fastq/single.fastq']


def test_paired_read_notfound():
    assert find_paired_reads('tests/dummy_fasta') == []


def test_make_path():
    make_path('tests/new_folder')
    assert os.path.isdir('tests/new_folder') is True
    shutil.rmtree('tests/new_folder')


def test_make_path_but_file_exists():
    with pytest.raises(OSError):
        make_path('tests/test_imports.py')


def test_make_path_nested():
    make_path('tests/new_folder/new_folder')
    assert os.path.isdir('tests/new_folder/new_folder')
    shutil.rmtree('tests/new_folder')


def test_filer_fastq():
    filelist = ['tests/dummy_fastq/test_R1.fastq', 'tests/dummy_fastq/test_R2.fastq',
                'tests/dummy_fastq/test2_S1_L001_R1_001.fastq.gz', 'tests/dummy_fastq/test2_S1_L001_R2_001.fastq.gz',
                'tests/dummy_fastq/test3_S55_L001_R1_001.fastq.gz', 'tests/dummy_fastq/test4_1.fastq',
                'tests/dummy_fastq/test4_2.fastq']
    assert filer(filelist) == {'tests/dummy_fastq/test', 'tests/dummy_fastq/test2',
                               'tests/dummy_fastq/test3', 'tests/dummy_fastq/test4'}


def test_combine_fastas():  # TODO: Actually inspect contents of created file.
    combinetargets(targets=glob.glob('tests/dummy_fasta/*fasta'), targetpath='tests')
    assert os.path.isfile('tests/combinedtargets.fasta')
    assert os.path.getsize('tests/combinedtargets.fasta')
    os.remove('tests/combinedtargets.fasta')


def test_combine_fastqs():
    combinetargets(targets=glob.glob('tests/dummy_fasta/*fastq'), targetpath='tests')
    assert os.path.isfile('tests/combinedtargets.fasta') is False

