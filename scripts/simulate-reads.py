#!/usr/bin/env python
import os
import sys
import pysam
import numpy as np
from scipy.sparse import lil_matrix
from numpy.random import multinomial, randint, poisson, choice
import getopt


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
            opts, args = getopt.getopt(argv[1:], "hF:g:p:m:N:r:e:", \
                                       ["help", "grouping-file", "parameter-file", "total-count", "error-rate"])
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
                target_file = value
            if option in ("-g", "--grouping-file"):
                grouping_file = value
            if option in ("-p", "--parameter-file"):
                param_file = value
            if option == "-m":
                model = int(value)
            if option in ("-N", "--total-count"):
                model = int(value)
            if option == "-r":
                model = int(value)
            if option in ("-e", "--error-rate"):
                model = int(value)
            if option == "-o":
                outfile = value

        # Check if the required options are given
        if model not in (1, 2, 3, 4):
            print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + 'Simulation model should be either 1, 2, 3 or 4.'
            return 2


        #
        # Main body
        #

        gname  = list()
        groups = list()
        with open(grouping_file) as fh:
            for curline in fh:
                item = curline.rstrip().split("\t")
                gname.append(item[0])
                tid_list = [ tid[t] for t in item[1:] ]
                groups.append(tid_list)
        gname = np.array(gname)
        gid = dict(zip(gname, np.arange(len(gname))))
        grp_conv_mat = lil_matrix((len(tname), len(gname)))
        for i in xrange(len(gname)):
            grp_conv_mat[groups[i], i] = 1.0
        grp_conv_mat = grp_conv_mat.tocsc()

        # Generate read counts

        emase_count_transcript = np.zeros((len(tname), 2))
        with open(param_file) as fh:
            curline = fh.next()
            item = curline.rstrip().split('\t')
            hname = item[1:]
            num_haps = len(hname)
            hid = dict(zip(hname, np.arange(num_haps)))
            for curline in fh:
                item = curline.rstrip().split('\t')
                emase_count_transcript[tid[item[0]], :] = map(float, item[1:3])
        emase_count_transcript.shape, emase_count_transcript.sum()

        if model == 1:
            gacount = grp_conv_mat.transpose() * emase_count_transcript
            gcount = gacount.sum(axis=1)
            theta = gcount / gcount.sum()
            gexpr = gcount > 0
            gcount_sim = multinomial(N, theta)[0]
            phi = np.zeros(gacount.shape)
            phi[gexpr, :] = gacount[gexpr, :] / gcount[gexpr, np.newaxis]
            phi[np.logical_not(gexpr), :] = np.ones(num_haps) / num_haps
            gacount_sim = np.zeros(gacount.shape)
            for g in xrange(num_genes):
                gacount_sim[g] = multinomial(gcount_sim[g], phi[g])
            tacount_sim = np.zeros(emase_count_transcript.shape)
            for g in xrange(num_genes):
                tindex = groups[g]
                num_isoforms = len(tindex)
                for h in xrange(num_haps):
                    delta = emase_count_transcript[tindex, h]
                    delta_sum = delta.sum()
                    if delta_sum > 0:
                        delta /= delta_sum
                    else:
                        delta = np.ones(num_isoforms) / num_isoforms
                    tacount_sim[tindex, h] = multinomial(gacount_sim[g, h], delta)
        elif model == 2:
            pass
        elif model == 3:
            pass
        else:  # if model == 4:
            pass


        # Generate reads

        f = pysam.FastaFile(target_file)
        target_info_file = os.path.splitext(target_file)[0] + '.info'
        trange = np.zeros(emase_count_transcript.shape)
        with open(target_info_file) as fh:
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
