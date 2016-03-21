#! /usr/bin/env python

import sys

extra = {}
try:
    from setuptools import setup
    if sys.version_info < (2, 7):
        extra["install_requires"] = ["argparse"]
    if sys.version_info >= (3,):
        extra["use_2to3"] = True
except ImportError:
    from distutils.core import setup
    if sys.version_info < (2, 7):
        extra["dependencies"] = ["argparse"]

setup(
    name="alchemy",
    packages=[
        "alchemy",
        "alchemy.oboparser",
        "alchemy.data"],
    version="0.1",
    description="Drug overlap",
    license="GPLv2",
    long_description=open('README.md').read(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    keywords="bioinformatics, chemical compounds, drug, ontology, enrichment",
    author="Andre Rendeiro",
    author_email="arendeiro@cemm.oeaw.ac.at",
    url="https://github.com/epigen/alchemy/",
    install_requires=[
        "numpy",
        "pandas>=0.14",
        "scipy",
        "statsmodels"
    ],
    entry_points={
        "console_scripts": [
            'alchemy = alchemy.alchemy:main'
        ],
    },
    package_data={
        'alchemy': [
            "../README.md",
            "../LICENSE",
            "../MANIFEST.in",
            'data/full_mapping.existing.csv',
            'data/chebi.obo']
    },
    include_package_data=True
)
