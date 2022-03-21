#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="OLCTools",
    version="1.2.3",
    packages=find_packages(),
    include_package_data=True,
    url="https://github.com/OLC-Bioinformatics/OLCTools",
    install_requires=['biopython',
                      'interop',
                      'sipprverse',
                      'pandas',
                      'numpy',
                      'geneseekr',
                      'psutil',
                      'urllib3',
                      'confindr>=0.5',
                      'pytest',
                      'pysam',
                      'xlsxwriter',
                      'seaborn',
                      'rauth',
                      'requests',
                      'selenium',
                      'chromedriver-autoinstaller',
                      'click']
)
