#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="shel",
    package_dir={"": "src/python"},
    packages=find_packages(where="src/python"),
)
