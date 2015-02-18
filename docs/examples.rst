==============
Real Use Cases
==============

Estimation of allele-specific expression from a diploid sample
--------------------------------------------------------------

Building Diploid Transcriptome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Seqnature incorporates known polymorphisms and short indels from genetically
diverse and heterozygous model organisms into reference genomes, and can
construct individualized haploid or diploid transcriptomes suitable for read
alignment by common aligners.

Aligning reads to the Diploid Transcriptome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EMASE has been extensively used with the aligner bowtie1. RNA-seq reads need to
be aligned simultaneously to a diploid transcriptome. Use bowtie aligner with
the following option::

    bowtie -a -best -strata -v 3 -m 100

Convert bam file to emase sparse 3-dim matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run EMASE
~~~~~~~~~

Estimation of allele-specific binding from ChIP-seq sample
----------------------------------------------------------

We assume you have a set of individualized genome and annotation files created by
Seqnature package. We also assume you have a bed file that specifies genomic regions
of your interest. First, you need to convert your bed file into a simple gtf format::

    Narayanan, add your command line here.

...

Finally, run::

    prepare-emase -G S1.fa,S2.fa -g S1.gtf,S2.gtf -s S1,S2 -o S1xS2

This will store the following files in the folder 'S1xS2'
* emase.pooled.transcriptome.fa
* emase.pooled.transcriptome.info
* bowtie1 index files

Now you can align your RNA-seq reads against the pooled bowtie index of target region::

    bowtie -q -a --best --strata --sam -v 3 S1xS2/bowtie.transcriptome S1xS2.fastq | samtools view -bS -F 4 - > bowtie.transcriptome.bam

It is now ready to run emase. We assume the read length is 100bp.::

    run-emase -i bowtie.transcriptome.bam -L S1xS2/emase.pooled.transcriptome.info -M 4 -c

Estimation of allele-specific expression from human RNA-seq sample
------------------------------------------------------------------

We assume you have a vcf file that contains phased variant calls. In general,
parent information is not available, so we are going to distinguish two alleles
with the names of L(eft) and R(ight).

First build individualized genomes using Seqnature package (link) like this.

...