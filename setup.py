#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools import find_packages


# load README.rst
def readme():
    with open('README.rst') as file:
        return file.read()


setup(
    name='csdemux',
    version='0.0.3',
    description='python3 wrapper for demultiplexing CS sequencing results',
    long_description=readme(),
    url='https://github.com/tomharrop/csdemux',
    author='Tom Harrop',
    author_email='twharrop@gmail.com',
    license='GPL-3',
    packages=find_packages(),
    install_requires=[
        'pandas>=1.1.0',
        'snakemake>=5.20.1'
    ],
    entry_points={
        'console_scripts': [
            'csdemux = csdemux.__main__:main'
            ],
    },
    scripts={
        'csdemux/src/calculate_hamming_distance.R',
        'csdemux/src/combine_step_logs.R',
        'csdemux/src/plot_adaptor_content.R',
        'csdemux/src/plot_barcode_content.R',
        'csdemux/src/plot_barcode_distance.R',
    },
    package_data={
        'csdemux': [
            'Snakefile',
            'README.rst'
        ],
    },
    zip_safe=False)
