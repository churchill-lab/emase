#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Copyright (c) 2015 Kwangbom Choi, The Jackson Laboratory
This software was developed by Kwangbom "KB" Choi in Gary Churchill's Lab.
This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this software. If not, see <http://www.gnu.org/licenses/>.
"""

import os
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import emase


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

# requirements = [
#     'numpy==1.8.2',
#     'numexpr==2.3.1',
#     'scipy==0.13.3',
#     'pysam>=0.6',
#     'cython>=0.13',
#     'tables==3.1.0',
#     'biopython>=1.63'
# ]

requirements = []
on_rtd = os.environ.get('READTHEDOCS', None)
if not on_rtd:
    requirements.append('numpy==1.8.2')
    requirements.append('numexpr==2.3.1')
    requirements.append('scipy==0.13.3')
    requirements.append('pysam>=0.6')
    requirements.append('cython>=0.13')
    requirements.append('tables==3.1.0')
    requirements.append('biopython>=1.63')

test_requirements = [
    'pytest'
]

setup(
    name='emase',
    version=emase.__version__,
    description="EMASE: Expectation-Maximization algorithm for Allele Specific Expression",
    long_description=readme + '\n\n' + history,
    author="Kwangbom \"KB\" Choi, Ph. D., The Jackson Laboratory",
    author_email='kb.choi@jax.org',
    url='https://github.com/churchill-lab/emase',
    packages=[
        'emase',
    ],
    package_dir={'emase':
                 'emase'},
    include_package_data=True,
    scripts=[
        'scripts/prepare-emase',
        'scripts/create-hybrid',
        'scripts/bam-to-emase',
        'scripts/combine-emase-files',
        'scripts/run-emase',
        'scripts/count-alignments',
        'scripts/get-common-alignments',
        'scripts/pull-out-unique-reads',
        'scripts/count-shared-multireads-pairwise',
        'scripts/simulate-reads'
    ],
    install_requires=requirements,
    license="GPLv3",
    zip_safe=False,
    keywords='emase',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
