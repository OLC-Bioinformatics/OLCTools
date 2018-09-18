#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="OLCTools",
    version="0.4.7",
    packages=find_packages(),
    include_package_data=True,
    author="Andrew Low",
    author_email="andrew.low@inspection.gc.ca",
    url="https://github.com/lowandrew/OLCTools",
    install_requires=['biopython', 'interop']
)
