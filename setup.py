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
    name="drugannotation",
    packages=["drugannotation"],
    version="0.1",
    description="Drug overlap",
    author="Andre Rendeiro",
    author_email="arendeiro@cemm.oeaw.ac.at",
    url="https://github.com/epigen/drugannotation/",
    package_data={'drugannotation': ['annotations/*.csv']}
)
