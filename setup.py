#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="OLCTools",
    version="0.3.96",
    packages=find_packages(),
    author="Andrew Low",
    author_email="andrew.low@inspection.gc.ca",
    url="https://github.com/lowandrew/OLCTools",
    install_requires=['biopython', 'interop']
)
