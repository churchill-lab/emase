========================================================================
EMASE: Expectation-Maximization algorithm for Allele Specific Expression
========================================================================

.. image:: https://badge.fury.io/py/emase.png
    :target: http://badge.fury.io/py/emase

.. image:: https://anaconda.org/kbchoi/emase/badges/version.svg
    :target: https://anaconda.org/kbchoi/emase

.. image:: https://travis-ci.org/churchill-lab/emase.png?branch=master
    :target: https://travis-ci.org/churchill-lab/emase

.. image:: https://readthedocs.org/projects/emase/badge/?version=latest
    :target: http://emase.readthedocs.org/en/latest/?badge=latest
    :alt: Documentation Status


Overview
--------

**EMASE** is a software program written in Python to quantify allele-specific
expression and gene expression simultaneously from RNA-seq data. EMASE takes in
the diploid transcriptome alignment BAM file and GTF file as inputs and
estimates expression abundance for each isoforms and each alleles using
Expectation Maxmization algorithm.

* Free software: GNU General Public License v3 (GPLv3)
* Documentation: https://emase.readthedocs.org.


Reference
---------

Allele-specific gene expression in F1 Hybrids from model organisms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we have F1 hybrids with parental genetic variants information, one can use
Seqnature to build strain specific genomes and extract diploid transcriptome.
RNA-seq alignment bam file obtained by aligbning RNA-seq reads to the diploid
transcriptome is used as input for EMASE.

* [Hierarchical Analysis of Multi-mapping RNA-Seq Reads Improves the Accuracy of Allele-specific Expression](http://www.biorxiv.org/content/early/2017/07/22/166900) Narayanan Raghupathy, Kwangbom Choi, Matthew J Vincent, Glen L Beane, Keith Sheppard, Steven C Munger, Ron Korstanje, Fernando Pardo Manuel de Villena, Gary A Churchill. Manuscript in review.


Applications
------------

Personalized ASE analysis
^^^^^^^^^^^^^^^^^^^^^^^^^

EMASE can be used to do personalized RNA-seq analysis in human. For this,  we use the phased genetic
variation (SNP and Indel) information to construct personalized diploid genome and align reads to the diploid transcriptome.

* [RNA-Seq Alignment to Individualized Genomes Improves Transcript Abundance Estimates in Multiparent Populations](http://www.genetics.org/content/198/1/59.short) Steven C. Munger, Narayanan Raghupathy, Kwangbom Choi, Allen K. Simons, Daniel M. Gatti, Douglas A. Hinerfeld, Karen L. Svenson, Mark P. Keller, Alan D. Attie, Matthew A. Hibbs, Joel H. Graber, Elissa J. Chesler and Gary A. Churchill. Genetics. 2014 Sep;198(1):59-73. doi: 10.1534/genetics.114.165886.

Allele-specific Binding using Chip-Seq in F1 Hybrids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although we explained the use of EMASE to quantify Allele-Specific Expression
from RNA-seq data, the tool can be used with other types of sequencing data. We
have successfully used EMASE to quantify allele-specific binding from ChIP-seq
data. While useing ChIP-seq, one needs to use diploid binding target sequences
instead of diploid transcriptome for alignment target sequences.

* [PRDM9 Drives Evolutionary Erosion of Hotspots in Mus musculus through Haplotype-Specific Initiation of Meiotic Recombination](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004916) Christopher L. Baker, Shimpei Kajita, Michael Walker, Ruth L. Saxl, Narayanan Raghupathy, Kwangbom Choi, Petko M. Petkov, Kenneth Paigen PLOS Genetics: published 08 Jan 2015 | info:doi/10.1371/journal.pgen.1004916
