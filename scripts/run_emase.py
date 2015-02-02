#!/usr/bin/env python
import sys
import getopt
from emase.EMfactory import EMfactory
from emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM

__author__ = 'Kwangbom "KB" Choi, Ph. D.'

help_message = '''

    Usage:
        run_emase.py -i <h5_file> -g <grp_file> -d <len_file> -o <out_file> \\
                     -b <min_uniq_reads> -p <pseudocount> -r <read_length> -m <max_iters> -t <tolerance>
    Input:
        <h5_file>        : Alignments stored in a PyTables HDF5 format
        <grp_file>       : Gene-to-transcript map (ENSMUSGxxx followed by a list of ENSMUSTyyy's)
        <len_file>       : File that contains transcript lengths
        <out_file>       : The name of output result file. (default: 'emase.isoforms')
        <min_uniq_reads> : The number of unique reads to be considered for convergence checking
        <pseudocount>    : Pseudocount for allele specificity (default: 0.0)
        <max_iters>      : The number of maximum iterations for EM (default: 1000)
        <tolerance>      : Tolerance for the termination of EM. (default: 0.01)
    Parameters:
        -h, --help: shows this help message
        -w: reports the posterior probability for each read

'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hi:g:d:b:p:r:m:t:o:w", ["help"])
        except getopt.error, msg:
            raise Usage(msg)

        # Default values of vars
        h5file = None
        grpfile = None
        lenfile = None
        outfile = 'emase.isoforms'
        min_uniq_reads = 0
        read_length = 100
        pseudocount = 0.2
        max_iters = 999
        tolerance = 0.01
        report_gene_counts = False
        report_posterior = False

        # option processing (change this later with optparse)
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option == '-i':
                h5file = value
            if option == '-g':
                grpfile = value
                report_gene_counts = True
            if option == '-d':
                lenfile = value
            if option == '-b':
                min_uniq_reads = int(value)
            if option == '-p':
                pseudocount = float(value)
            if option == '-r':
                read_length = int(value)
            if option == '-m':
                max_iters = int(value)
            if option == '-t':
                tolerance = float(value)
            if option == '-o':
                outfile = value
            if option == '-w':
                report_posterior = True

        # Check if the required options are given
        if h5file is None:  # If alignment file is not given
            raise Usage(help_message)

        #
        # Main body
        #

        alignments = APM(h5file=h5file, grpfile=grpfile)
        alignments.reset()
        em_factory = EMfactory(alignments, lenfile=lenfile, read_length=read_length)
        em_factory.prepare(pseudocount)
        em_factory.run(tol=tolerance, max_iters=max_iters, min_uniq_reads=min_uniq_reads, verbose=True)
        em_factory.report_depths(filename="%s.tpm" % outfile, tpm=True)
        em_factory.report_effective_read_counts(filename="%s.effective_read_counts" % outfile)
        if report_posterior:
            em_factory.export_posterior_probability(filename="%s.posterior.h5" % outfile)
        if report_gene_counts:
            grp_outfile = outfile.replace('isoforms', 'genes')
            em_factory.report_depths(filename="%s.tpm" % grp_outfile, tpm=True, grp_wise=True)
            em_factory.report_effective_read_counts(filename="%s.effective_read_counts" % grp_outfile, grp_wise=True)

        #
        # End of main body
        #

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        return 2


if __name__ == "__main__":
    sys.exit(main())
