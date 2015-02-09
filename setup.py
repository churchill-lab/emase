#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import emase


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    'numpy>=1.7',
    'scipy>=0.12',
    'pysam>=0.6',
    'tables>=3.0'
]

test_requirements = [
    'pytest'
]

setup(
    name='emase',
    version=emase.__version__,
    description="EMASE: Expectation-Maximization algorithm for Allele Specific Expression",
    long_description=readme + '\n\n' + history,
    author="Kwangbom \"KB\" Choi",
    author_email='kb.choi@jax.org',
    url='https://github.com/jax-cgd/emase',
    packages=[
        'emase',
    ],
    package_dir={'emase':
                 'emase'},
    include_package_data=True,
    scripts=[
        'scripts/prepare-emase',
        'scripts/bam-to-emase',
        'scripts/combine-emase-files',
        'scripts/run-emase',
        'scripts/count-alignments',
        'scripts/count-shared-multireads-pairwise'
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
