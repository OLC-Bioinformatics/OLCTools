import shutil
import os
import pytest
from accessoryFunctions.accessoryFunctions import *


def test_dependency_check_good():
    assert dependency_check('ls') is True


def test_dependency_check_bad():
    assert dependency_check('aaaaaaddsfsdss') is False


def test_paired_read_find():
    assert find_paired_reads('tests/dummy_fastq') == [['tests/dummy_fastq/test_R1.fastq',
                                                       'tests/dummy_fastq/test_R2.fastq']]


def test_paired_read_notfound():
    assert find_paired_reads('tests/dummy_fasta') == []


def test_make_path():
    make_path('tests/new_folder')
    assert os.path.isdir('tests/new_folder') is True
    shutil.rmtree('tests/new_folder')


def test_make_path_butfile():
    with pytest.raises(OSError):
        make_path('tests/test_mash.py')


def test_make_path_nested():
    make_path('tests/new_folder/new_folder')
    assert os.path.isdir('tests/new_folder/new_folder')
    shutil.rmtree('tests/new_folder')


def test_filer_fastq():
    filelist = ['tests/dummy_fastq/test_R1.fastq', 'tests/dummy_fastq/test_R2.fastq']
    assert filer(filelist) == {'tests/dummy_fastq/test'}
