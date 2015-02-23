=====
Usage
=====

To use EMASE in a project
~~~~~~~~~~~~~~~~~~~~~~~~~

In python scripts, we can load 'emase' as a module::

    import emase


To run EMASE
~~~~~~~~~~~~

We need to first build a pooled transcriptome and prepare required files for EMASE::

    prepare-emase -G ${GENOME1},${GENOME2} -g ${GTF1},${GTF2} -o ${EMASE_DIR} -s ${SUFFIX1},${SUFFIX2} -m

'prepare-emase' also generates the following files for EMASE::

    ${EMASE_DIR}/emase.pooled.transcriptome.fa
    ${EMASE_DIR}/emase.pooled.transcriptome.info  <== Used as ${TINFO_FILE} in the next steps
    ${EMASE_DIR}/emase.gene2transcript.tsv        <== Used as ${GROUP_FILE} in the next steps
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
