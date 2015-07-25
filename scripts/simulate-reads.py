#!/usr/bin/env python
import pysam
import numpy as np
from scipy.sparse import lil_matrix
from numpy.random import multinomial, randint, poisson, choice


__author__ = 'Kwangbom "KB" Choi, Ph. D.'


help_message = '''
    Usage:
        simulate-reads -F <target_fasta> -p <param_file> -m <simulation_model> -N <total_count> [ -r <read_len> -e <err_rate> ]

    Input:
        <target_fasta> : A target fasta file
        <hap_list>    : Names of haplotypes to be used instead (comma delimited, in the order of genomes)
        <out_file>    : Output file name (default: './emase.pooled.targets.fa')

    Parameters:
        --help, -h : shows this help message
        --create-bowtie-index : builds bowtie1 index

'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def introduce_mutation(seq, pos, new_base):
    return seq[:pos] + new_base + seq[pos+1:]


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hF:s:o:", ["help", "create-bowtie-index"])
        except getopt.error, msg:
            raise Usage(msg)

        # Default values of vars
        read_len = 100
        lamb = 0.1
        mut = {'A':('T', 'G', 'C'), \
               'T':('G', 'C', 'A'), \
               'G':('C', 'A', 'T'), \
               'C':('A', 'T', 'G')}

        # option processing (change this later with optparse)
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option == "-F":
                fastalist = value.split(',')
                num_haps = len(fastalist)
            if option == "-s":
                haplist = value.split(',')
            if option == "--create-bowtie-index":
                build_bowtie_index = True
            if option == "-o":
                outfile = value

        # Check if the required options are given


        #
        # Main body
        #

        # Generate read counts
        emase_count_transcript = np.zeros((len(tname), 2))
        with open(param_file) as fh:
            fh.next()
            for curline in fh:
                item = curline.rstrip().split('\t')
                emase_count_transcript[tid[item[0]], :] = map(float, item[1:3])
        emase_count_transcript.shape, emase_count_transcript.sum()

        # Simulate reads

        f = pysam.FastaFile(target)
        trange = np.zeros(emase_count_transcript.shape)
        with open('/data/kbchoi/data/mm10/R75-REL1410/B6xCAST/emase.pooled.transcriptome.info') as fh:
            for curline in fh:
                item = curline.rstrip().split('\t')
                t, h = item[0].split('_')
                trange[tid[t], hid[h]] = np.int(item[1])
        trange = trange - read_len + 1

        with open(outfile, 'w') as fhout:
            with open('emase.M1.simulated.transcripts.read_counts') as fh:
                curline = fh.next()
                item = curline.rstrip().split('\t')
                hname = item[1:]
                hid = dict(zip(hname, np.arange(len(hname))))
                for curline in fh:
                    item = curline.rstrip().split('\t')
                    t = item[0]
                    hcount = map(int, item[1:])
                    for h in xrange(num_haps):
                        th = '%s_%s' % (t, hname[h])
                        trange_max = trange[tid[t], h]
                        for k in xrange(hcount[h]):
                            start = randint(0, trange_max)
                            end = start + read_len
                            seq = f.fetch(th, start, end)
                            nerr = poisson(lamb)
                            errstr = ''
                            if nerr > 0:
                                errloc = sorted(choice(read_len, nerr, replace=False))
                                errstr = ':'
                                for epos in errloc:
                                    from_base = seq[epos]
                                    to_base = mut[from_base][randint(3)]
                                    seq = introduce_mutation(seq, epos, to_base)
                                    errstr += '[%d]%s->%s;' % (epos+1, from_base, to_base)
                            fhout.write(">%s:%d-%d%s\n%s\n" % (th, start+1, end, errstr, seq))

        #
        # End of Main body
        #

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        return 2


if __name__ == "__main__":
    sys.exit(main())
