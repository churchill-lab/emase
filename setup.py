#!/usr/bin/env python
import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read().replace(".. :changelog:", "")

def get_emase_version():
    sys.path.insert(0, "emase")
    import version
    return version.__version__

requirements = []
on_rtd = os.environ.get("READTHEDOCS", None)
if not on_rtd:
    with open("requirements.txt") as requirements_file:
        requirements_lines = requirements_file.readlines()
        for line in requirements_lines:
            requirements.append(line)

test_requirements = ["pytest"]

setup(
    name="emase",
    version=get_emase_version(),
    description="EMASE: Expectation-Maximization algorithm for Allele Specific Expression",
    long_description=readme + "\n\n" + history,
    author='Kwangbom "KB" Choi, Ph. D., The Jackson Laboratory',
    author_email="kb.choi@jax.org",
    url="https://github.com/churchill-lab/emase",
    packages=[
        "emase",
    ],
    package_dir={"emase": "emase"},
    entry_points={
        "console_scripts": [
            "emase = emase.commands:app",
        ]
    },
    include_package_data=True,
    scripts=[
        "scripts/prepare-emase",
        "scripts/create-hybrid",
        "scripts/bam-to-emase",
        "scripts/combine-emase-files",
        "scripts/run-emase",
        "scripts/count-alignments",
        "scripts/get-common-alignments",
        "scripts/pull-out-unique-reads",
        "scripts/count-shared-multireads-pairwise",
        "scripts/simulate-reads",
    ],
    install_requires=requirements,
    license="GPLv3",
    zip_safe=False,
    keywords="emase",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.10",
    ],
    test_suite="tests",
    tests_require=test_requirements,
)
