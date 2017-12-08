from biotools import kmc
import os


def test_kmc_call_paired():
    out, err, cmd = kmc.kmc(forward_in='tests/dummy_fastq/test_R1.fastq',
                            database_name='tests/kmc_db',
                            returncmd=True)
    assert cmd == 'kmc -k31 -ci1  @tmp/filelist.txt tests/kmc_db tmp'
    os.remove('tests/kmc_db.kmc_pre')
    os.remove('tests/kmc_db.kmc_suf')


def test_kmc_call_single():
    out, err, cmd = kmc.kmc(forward_in='tests/dummy_fastq/single.fastq',
                            database_name='tests/kmc_db',
                            returncmd=True)
    assert cmd == 'kmc -k31 -ci1  tests/dummy_fastq/single.fastq tests/kmc_db tmp'
    os.remove('tests/kmc_db.kmc_pre')
    os.remove('tests/kmc_db.kmc_suf')


def test_kmc_call_paired_k22():
    out, err, cmd = kmc.kmc(forward_in='tests/dummy_fastq/test_R1.fastq',
                            database_name='tests/kmc_db',
                            returncmd=True, k=22)
    assert cmd == 'kmc -k22 -ci1  @tmp/filelist.txt tests/kmc_db tmp'
    os.remove('tests/kmc_db.kmc_pre')


def test_kmc_call_kwargs():
    out, err, cmd = kmc.kmc(forward_in='tests/dummy_fastq/test_R1.fastq',
                            database_name='tests/kmc_db',
                            returncmd=True, m='2')
    assert cmd == 'kmc -k31 -ci1  -m2 @tmp/filelist.txt tests/kmc_db tmp'
    os.remove('tests/kmc_db.kmc_pre')
    os.remove('tests/kmc_db.kmc_suf')


def test_kmc_dump_call():
    out, err, cmd = kmc.dump(database='tests/kmc_dbs/db_1', output='tests/kmc_dump', returncmd=True)
    assert cmd == 'kmc_tools dump -ci1 -cx250 tests/kmc_dbs/db_1 tests/kmc_dump'
    os.remove('tests/kmc_dump')


def test_kmc_union_call():
    out, err, cmd = kmc.union(database_1='tests/kmc_dbs/db_1', database_2='tests/kmc_dbs/db_2',
                              results='tests/kmc_db', returncmd=True)
    assert cmd == 'kmc_tools union tests/kmc_dbs/db_1 tests/kmc_dbs/db_2 tests/kmc_db'
    os.remove('tests/kmc_db.kmc_pre')
    os.remove('tests/kmc_db.kmc_suf')


def test_kmc_subtract_call():
    out, err, cmd = kmc.subtract(database_1='tests/kmc_dbs/db_1', database_2='tests/kmc_dbs/db_2',
                                 results='tests/kmc_db', returncmd=True)
    assert cmd == 'kmc_tools kmers_subtract tests/kmc_dbs/db_1 tests/kmc_dbs/db_2 -ci1 tests/kmc_db'
    os.remove('tests/kmc_db.kmc_pre')
    os.remove('tests/kmc_db.kmc_suf')


def test_kmc_intersect_call():
    out, err, cmd = kmc.intersect(database_1='tests/kmc_dbs/db_1', database_2='tests/kmc_dbs/db_2',
                                  results='tests/kmc_db', returncmd=True)
    assert cmd == 'kmc_tools intersect tests/kmc_dbs/db_1 tests/kmc_dbs/db_2 tests/kmc_db'
    os.remove('tests/kmc_db.kmc_pre')
    os.remove('tests/kmc_db.kmc_suf')
