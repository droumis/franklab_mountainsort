#!/usr/bin/env python3
from setuptools import find_packages, setup

INSTALL_REQUIRES = ['franklab_msdrift', 'franklab_mstaggedcuration']
TESTS_REQUIRE = []

setup(
    name='franklab_mountainsort',
    version='0.1.1.dev0',
    license='',
    description=(''),
    author='',
    author_email='',
    url='https://github.com/LorenFrankLab/franklab_mountainsort',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    tests_require=TESTS_REQUIRE,
)
