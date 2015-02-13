#!/usr/bin/env python
import numpy as np
import time
from scipy.sparse import eye, lil_matrix
from .AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM

__author__ = 'Kwangbom "KB" Choi, Ph. D.'

def enum(**enums):
    return type('Enum', (), enums)

class EMfactory:

    Model = enum(GENE_ISOFORM_ALLELE_MODEL=0, GENE_ALLELE_ISOFORM_MODEL=1, GENE_COMMUNITY_MODEL=2, RSEM_MODEL=3)

    def __init__(self, alignments, lenfile=None, read_length=100):
        self.alignments = alignments
        self.allelic_expression = None
        self.grp_conv_mat = None
        if self.alignments.num_groups > 0:
            self.grp_conv_mat = lil_matrix((self.alignments.num_loci, self.alignments.num_groups))
            for i in xrange(self.alignments.num_groups):
                self.grp_conv_mat[self.alignments.groups[i], i] = 1.0
            self.grp_conv_mat = self.grp_conv_mat.tocsc()
        self.t2t_mat = eye(self.alignments.num_loci, self.alignments.num_loci)
        self.t2t_mat = self.t2t_mat.tolil()
        for tid_list in self.alignments.groups:
            for ii in xrange(len(tid_list)):
                for jj in xrange(ii):
                    i = tid_list[ii]
                    j = tid_list[jj]
                    self.t2t_mat[i, j] = 1
                    self.t2t_mat[j, i] = 1
        self.t2t_mat = self.t2t_mat.tocsc()
        self.target_lengths = None
        if lenfile is not None:
            hid = dict(zip(self.alignments.hname, np.arange(len(self.alignments.hname))))
            self.target_lengths = np.zeros((self.alignments.num_loci, self.alignments.num_haplotypes))
            with open(lenfile) as fh:
                for curline in fh:
                    item = curline.rstrip().split("\t")
                    locus, hap = item[0].split("_")
                    self.target_lengths[self.alignments.lid[locus], hid[hap]] = max(float(item[1]), 1.0)
            self.target_lengths = self.target_lengths.transpose() / read_length  # lengths in terms of read counts
            if not np.all(self.target_lengths > 0.0):
                raise RuntimeError('There exist transcripts missing length information.')

    def prepare(self, pseudocount=0.0):
        """
        Initialize the probability of read origin according to the alignment profile
        :param pseudocount: Uniform prior for allele specificity estimation
        :return: Nothing (as it performs an in-place operations)
        """
        # self.alignments.initialize()
        self.alignments.normalize_reads(axis=APM.Axis.READ)
        self.allelic_expression = self.alignments.sum(axis=APM.Axis.READ)
        if self.target_lengths is not None: # allelic_expression will be at depth-level
            self.allelic_expression = np.divide(self.allelic_expression, self.target_lengths)
        if pseudocount > 0.0:  # pseudocount is at depth-level
            orig_allelic_expression_sum = self.allelic_expression.sum()
            nzloci = np.nonzero(self.allelic_expression)[1]
            self.allelic_expression[:, nzloci] += pseudocount
            self.allelic_expression *= (orig_allelic_expression_sum / self.allelic_expression.sum())  # original depth scale

    def get_allelic_expression(self, at_group_level=False):
        if at_group_level:
            return self.allelic_expression * self.grp_conv_mat
        else:
            return self.allelic_expression.copy()

    def update_probability_at_read_level(self, model=1):
        self.alignments.reset()
        if model == self.Model.GENE_ISOFORM_ALLELE_MODEL:
            self.alignments.multiply(self.allelic_expression, axis=APM.Axis.READ)
            self.alignments.normalize_reads(axis=APM.Axis.LOCUS)  # Locus-wise normalization
            self.alignments.multiply(self.allelic_expression.sum(axis=0), axis=APM.Axis.HAPLOTYPE)
            self.alignments.normalize_reads(axis=APM.Axis.GROUP, grouping_mat=self.t2t_mat)
            self.alignments.multiply((self.allelic_expression * self.t2t_mat).sum(axis=0), axis=APM.Axis.HAPLOTYPE)
            self.alignments.normalize_reads(axis=APM.Axis.READ)
        elif model == self.Model.GENE_ALLELE_ISOFORM_MODEL:
            pass
        elif model == self.Model.GENE_COMMUNITY_MODEL:
            pass
        elif model == self.Model.RSEM_MODEL:
            self.alignments.multiply(self.allelic_expression, axis=APM.Axis.READ)
            self.alignments.normalize_reads(axis=APM.Axis.READ)

    def update_allelic_expression(self, model=1):
        """
        A single EM step
        :param model:
        :return:
        """
        err_states = np.seterr(all='warn')
        err_states = np.seterr(**err_states)
        self.update_probability_at_read_level(model)
        self.allelic_expression = self.alignments.sum(axis=APM.Axis.READ)
        if self.target_lengths is not None:
            self.allelic_expression = np.divide(self.allelic_expression, self.target_lengths)

    def run(self, model=1, tol=0.01, max_iters=999, min_uniq_reads=1, verbose=True):
        """
        Run EM iterations
        :param model:
        :param tol:
        :param max_iters:
        :param min_uniq_reads:
        :param verbose:
        :return:
        """
        if verbose:
            print
            print "Iter No  Time (hh:mm:ss)  Total error (depth)  Max error (%)  Locus of max error  Allele expression change"
            print "-------  ---------------  -------------------  -------------  ------------------  ------------------------"
        num_iters = 0
        err_max = 1.0
        time0 = time.time()
        num_uniq_reads = self.alignments.count_unique_reads()
        errchk_locs = np.where(num_uniq_reads > min_uniq_reads - np.nextafter(0, 1))
        while err_max > tol and num_iters < max_iters:
            prev_allelic_expression = self.get_allelic_expression()
            self.update_allelic_expression(model=model)
            curr_allelic_expression = self.get_allelic_expression()
            err = np.abs(curr_allelic_expression - prev_allelic_expression)
            err_sum = err.sum()
            err_locs = np.nonzero(prev_allelic_expression)
            err[err_locs] /= np.maximum(prev_allelic_expression[err_locs], 1.0)
            err_masked = err[errchk_locs]
            err_max = err_masked.max()
            err_maxlocus = np.unravel_index(err_masked.argmax(), err_masked.shape)
            err_max_hid = errchk_locs[0][err_maxlocus]
            err_max_lid = errchk_locs[1][err_maxlocus]
            num_iters += 1
            if verbose:
                time1 = time.time()
                delmin, s = divmod(int(time1 - time0), 60)
                h, m = divmod(delmin, 60)
                print " %5d      %4d:%02d:%02d     %16.2f     %9.2f%%     %s  %s: %.2f ==> %.2f" % \
                      (num_iters, h, m, s, err_sum, err_max * 100,
                       self.alignments.lname[err_max_lid], self.alignments.hname[err_max_hid],
                       prev_allelic_expression[err_max_hid, err_max_lid],
                       curr_allelic_expression[err_max_hid, err_max_lid])

    def report_effective_read_counts(self, filename, grp_wise=False, reorder='as-is'):
        """
        Get counts
        :param filename:
        :param grp_wise:
        :param reorder:
        :return:
        """
        effective_read_counts = self.alignments.sum(axis=APM.Axis.READ)
        if grp_wise:
            lname = self.alignments.gname
            effective_read_counts = effective_read_counts * self.grp_conv_mat
        else:
            lname = self.alignments.lname
        total_read_counts = effective_read_counts.sum(axis=0)
        if reorder == 'decreasing':
            report_order = np.argsort(total_read_counts.flatten())
            report_order = report_order[::-1]
        elif reorder == 'increasing':
            report_order = np.argsort(total_read_counts.flatten())
        elif reorder == 'as-is':
            report_order = np.arange(len(lname))  # report in the original locus order
        cntdata = np.vstack((effective_read_counts, total_read_counts))
        fhout = open(filename, 'w')
        fhout.write("locus\t" + "\t".join(self.alignments.hname) + "\ttotal\n")
        for locus_id in report_order:
            fhout.write("\t".join([lname[locus_id]] + map(str, cntdata[:, locus_id].ravel())) + "\n")
        fhout.close()

    def report_depths(self, filename, tpm=True, grp_wise=False, reorder='as-is'):
        if grp_wise:
            lname = self.alignments.gname
            depths = self.allelic_expression * self.grp_conv_mat
        else:
            lname = self.alignments.lname
            depths = self.allelic_expression
        if tpm:
            depths = depths * (1000000.0 / depths.sum())
        total_depths = depths.sum(axis=0)
        if reorder == 'decreasing':
            report_order = np.argsort(total_depths.flatten())
            report_order = report_order[::-1]
        elif reorder == 'increasing':
            report_order = np.argsort(total_depths.flatten())
        elif reorder == 'as-is':
            report_order = np.arange(len(lname))  # report in the original locus order
        cntdata = np.vstack((depths, total_depths))
        fhout = open(filename, 'w')
        fhout.write("locus\t" + "\t".join(self.alignments.hname) + "\ttotal\n")
        for locus_id in report_order:
            fhout.write("\t".join([lname[locus_id]] + map(str, cntdata[:, locus_id].ravel())) + "\n")
        fhout.close()

    def export_posterior_probability(self, filename, title="Posterior Probability"):
        self.alignments.save(h5file=filename, title=title)


if __name__ == "__main__":
    pass  # TODO: Put some simple test
