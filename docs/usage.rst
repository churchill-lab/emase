=====
Usage
=====

To use EMASE in a new project
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In python scripts, we can load 'emase' as a module::

    import emase


To run EMASE on command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need to first process reference genome::

    prepare-emase -G ${REF_GENOME} -g ${REF_GTF} -o ${REF_DIR} -m

'prepare-emase' generates the following files for EMASE::

    ${REF_DIR}/emase.transcriptome.fa
    ${REF_DIR}/emase.transcriptome.info         <== Used as ${TID_FILE} in the next steps
    ${REF_DIR}/emase.gene2transcript.tsv        <== Used as ${GROUP_FILE} in the next steps
    ${REF_DIR}/bowtie.transcriptome.1.ebwt
    ${REF_DIR}/bowtie.transcriptome.2.ebwt
    ${REF_DIR}/bowtie.transcriptome.3.ebwt
    ${REF_DIR}/bowtie.transcriptome.4.ebwt
    ${REF_DIR}/bowtie.transcriptome.rev.1.ebwt
    ${REF_DIR}/bowtie.transcriptome.rev.2.ebwt

Then build a pooled transcriptome and prepare required files for EMASE::

    prepare-emase -G ${GENOME1},${GENOME2} -g ${GTF1},${GTF2} -s ${SUFFIX1},${SUFFIX2} -o ${EMASE_DIR}

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
        | samtools view -bS -F 4 - > ${BAM_FILE}

Before running EMASE, we need to convert the bam file into the emase format::

    bam-to-emase -a ${BAM_FILE} -i ${TID_FILE} -s ${SUFFICE1},${SUFFIX2} -o ${EMASE_FILE}

Finally, to run EMASE::

    run-emase -i ${EMASE_FILE} -g ${GROUP_FILE} -L ${TINFO_FILE} -M ${MODEL} -o ${OUTBASE} \
              -r ${READLEN} -p ${PSEUDOCOUNT} -m ${MAX_ITERS} -t ${TOLERANCE}

'run-emase' outputs the following files::

    ${OUTBASE}.isoforms.effective_read_counts
    ${OUTBASE}.isoforms.tpm
    ${OUTBASE}.genes.effective_read_counts
    ${OUTBASE}.genes.tpm
