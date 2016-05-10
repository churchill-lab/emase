=====
Usage
=====

To use EMASE in a new project
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In python scripts, we can load 'emase' as a module::

    import emase

Or::

    from emase import AlignmentMatrixFactory as AMF
    from emase import AlignmentPropertyMatrix as APM
    from emase import EMfactory

To run EMASE on command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Note:** We will assume you installed EMASE in its own conda virtual environment. First of all, you have to "activate" the virtual environment by doing the following::

    source activate emase

The first step of the pipeline is to process reference genome::

    prepare-emase -G ${REF_GENOME} -g ${REF_GTF} -o ${REF_DIR} -m --no-bowtie-index

'prepare-emase' generates the following files for the reference genome::

    ${REF_DIR}/emase.transcriptome.fa
    ${REF_DIR}/emase.transcriptome.info         <== Used as ${TID_FILE} in the following steps
    ${REF_DIR}/emase.gene2transcripts.tsv       <== Used as ${GROUP_FILE} in the following steps

Then build a pooled transcriptome and prepare required files for EMASE::

    create-hybrid -G ${GENOME1},${GENOME2} -g ${GTF1},${GTF2} \
                  -s ${SUFFIX1},${SUFFIX2} -o ${EMASE_DIR}

Now the following files will be available::

    ${EMASE_DIR}/emase.pooled.transcriptome.fa
    ${EMASE_DIR}/emase.pooled.transcriptome.info  <== Used as ${TINFO_FILE} in the next steps
    ${EMASE_DIR}/bowtie.transcriptome.1.ebwt
    ${EMASE_DIR}/bowtie.transcriptome.2.ebwt
    ${EMASE_DIR}/bowtie.transcriptome.3.ebwt
    ${EMASE_DIR}/bowtie.transcriptome.4.ebwt
    ${EMASE_DIR}/bowtie.transcriptome.rev.1.ebwt
    ${EMASE_DIR}/bowtie.transcriptome.rev.2.ebwt

RNA-seq reads should be aligned against the pooled transcriptome::

    bowtie -q -a --best --strata --sam -v 3 ${EMASE_DIR}/bowtie.transcriptome ${FASTQ} \
        | samtools view -bS - > ${BAM_FILE}

Before running EMASE, we need to convert the bam file into the emase format::

    bam-to-emase -a ${BAM_FILE} -i ${TID_FILE} -s ${SUFFICE1},${SUFFIX2} -o ${EMASE_FILE}

For paired-end data, perform upto this step with R1 and R2 end independently, and get their common alignments::

    get-common-alignments -i ${EMASE_FILE_R1},${EMASE_FILE_R2} -o ${EMASE_FILE}

Finally, to run EMASE::

    run-emase -i ${EMASE_FILE} -g ${GROUP_FILE} -L ${TINFO_FILE} -M ${MODEL} -o ${OUTBASE} \
              -r ${READLEN} -p ${PSEUDOCOUNT} -m ${MAX_ITERS} -t ${TOLERANCE}

'run-emase' outputs the following files::

    ${OUTBASE}.isoforms.expected_read_counts
    ${OUTBASE}.isoforms.tpm
    ${OUTBASE}.genes.expected_read_counts
    ${OUTBASE}.genes.tpm

