========
Examples
========

Estimation of allele-specific expression from a diploid sample
--------------------------------------------------------------

Building Diploid Transcriptome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Seqnature incorporates known polymorphisms and short indels from
genetically diverse and heterozygous model organisms into reference
genomes, and can construct individualized haploid or diploid
transcriptomes suitable for read alignment by common aligners.

Aligning reads to the Diploid Transcriptome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EMASE has been extensively used with the aligner bowtie1. RNA-seq reads
need to be aligned simultaneously to a diploid transcriptome. Use bowtie
aligner with the following option::

    bowtie -a -best -strata -v 3 -m 100

Convert bam file to emase sparse 3-dim matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run EMASE
~~~~~~~~~

Results
~~~~~~~

*Files
emase.isoforms.effective_read_counts
emase.isoforms.tpm
emase.genes.effective_read_counts
emase.genes.tpm
