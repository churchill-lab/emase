===============================
Overview
===============================

.. image:: https://badge.fury.io/py/emase.png
    :target: http://badge.fury.io/py/emase

.. image:: https://travis-ci.org/jax-cgd/emase.png?branch=master
        :target: https://travis-ci.org/jax-cgd/emase

.. image:: https://pypip.in/d/emase/badge.png
        :target: https://pypi.python.org/pypi/emase


EMASE: Expectation-Maximization algorithm for Allele Specific Expression
------------------------------------------------------------------------
Narayanan Raghupathy, Kwangbom Choi, Steve Munger, and Gary Churchill

* Free software: MIT license
* Documentation: https://emase.readthedocs.org.

What is EMASE?
~~~~~~~~~~~~~~

EMASE is a software program written in Python to quantify allele-specific
expression and gene expression simultaneously from RNA-seq data. EMASE takes in
the diploid transcriptome alignment BAM file and GTF file as inputs and
estimates expression abundance for each isoforms and each alleles using
Expectation Maxmization algorithm.

Why Use EMASE?
~~~~~~~~~~~~~~

Current RNA-seq analysis pipeline employ two steps to quantify gene expression
and allele-specific expression (ASE); gene expression is estimated from all
read alignments, while ASE is assessed separately by using only reads that
overlap known SNP locations.

Large-scale genome sequencing efforts have characterized millions of genetic
variants across in human and model organisms. However development of tools that
can effectively utilize this individual/strain-specific variation to inform
quantitation of gene expression abundance have lagged behind.

In F1 hybrids from model organisms, EMASE allows us to utilize parental
strain-specific genetic variation in RNA-seq analysis to quantify gene
expression and allele-specific expression (ASE) simultaneously

In humans, EMASE allows us to utilize the individual's genetic variation in
doing personalized RNA-seq analysis and quantify gene expression and
allele-specific expression (ASE) simultaneously

Briefly, EMASE: EM for allele-specific expression, uses individualized diploid
genomes/transcriptomes adjusted for known genetic variations and quantifies
allele-specific gene expression and total gene expressionsimultaneously. The EM
algorithm employed in EMASE models multi-reads at the level of gene, isoform,
and allele and apportions them probabilistically.

Earlier, we developed Seqnature to utilize known genetic variations, bot SNPs
and Indels, to build indivdualized genomes and adjust annotations. One can use
Seqnature to create individualized diploid transcritpme and align RNA-seq reads
simultaneously to the diploid transcriptome and get alignment file in BAM
format. This diploid BAM file can be used as input to EMASE.

Applications
~~~~~~~~~~~~

Allele-specific gene expression in F1 Hybrids from model organisms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we have F1 hybrids with parental genetic variants information, one can use
Seqnature to build strain specific genomes and extract diploid transcriptome.
RNA-seq alignment bam file obtained by aligbning RNA-seq reads to the diploid
transcriptome is used as input for EMASE.

Personalized ASE analysis in Human
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use EMASE for personalized RNA-seq analysis in human, we need phased genetic
variation information.

Allele-specific Binding using Chip-Seq in F1 Hybrids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although we explained the use of EMASE to quantify Allele-Specific Expression
from RNA-seq data, the tool can be used with other types of sequencing data. We
have successfully used EMASE to quantify allele-specific binding from ChIP-seq
data. While useing ChIP-seq, one needs to use diploid binding target sequences
instead of diploid transcriptome for alignment target sequences.

Mining Diploid alignments and alignment probabilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Features
~~~~~~~~

* TODO
