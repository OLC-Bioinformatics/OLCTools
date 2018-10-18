#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="OLCTools",
    version="0.5.4",
    packages=find_packages(),
    include_package_data=True,
    author="Andrew Low",
    author_email="andrew.low@canada.ca",
    url="https://github.com/lowandrew/OLCTools",
    install_requires=['biopython',
                      'interop',
                      'sipprverse',
                      'pandas',
                      'numpy',
                      'geneseekr',
                      'psutil',
                      'urllib3',
                      'confindr',
                      'pytest',
                      'pysam',
                      'xlsxwriter',
                      'seaborn',
		      'rauth',
		      'click']
)
