#!/usr/bin/env python

import os
import pysam
import tables
from scipy.sparse import coo_matrix
import numpy as np
import struct


class AlignmentMatrixFactory():

    def __init__(self, alnfile): #, haplotype_names=None, locus_names=None, read_names=None):
        self.alnfile = alnfile
        self.hname = None
        self.lname = None
        self.rname = None
        self.tmpfiles = None
        # self.num_alignments = None

    def prepare(self, haplotypes, loci, delim='_', outdir=None):
        if len(haplotypes) > 0:  # Suffices given
            self.hname = haplotypes
        else:  # Suffix not given
            self.hname = ['h0']
        self.lname = loci
        self.rname = set()
        fh = pysam.Samfile(self.alnfile, 'rb')
        for aln in fh.fetch(until_eof=True):
            self.rname.add(aln.qname)
        self.rname = np.array(sorted(list(self.rname)))
        num_loci = len(self.lname)
        num_reads = len(self.rname)
        lid = dict(zip(self.lname, np.arange(num_loci)))
        rid = dict(zip(self.rname, np.arange(num_reads)))
        self.tmpfiles = dict.fromkeys(self.hname)
        if outdir is None:
            outdir = os.path.dirname(self.alnfile)
        fhout = dict.fromkeys(self.hname)
        for hap in self.hname:
            outfile = os.path.join(outdir, "%s_%d.bin" % (hap, os.getpid()))
            self.tmpfiles[hap] = outfile
            fhout[hap] = open(outfile, "wb")
        fh = pysam.Samfile(self.alnfile, 'rb')
        if len(haplotypes) > 0:  # Suffices given
            for aln in fh.fetch(until_eof=True):
                if aln.flag != 4 and aln.flag != 8:
                    locus, hap = fh.getrname(aln.tid).split(delim)
                    fhout[hap].write(struct.pack('>I', rid[aln.qname]))
                    fhout[hap].write(struct.pack('>I', lid[locus]))
        else:  # Suffix not given
            hap = self.hname[0]
            for aln in fh.fetch(until_eof=True):
                if aln.flag != 4 and aln.flag != 8:
                    locus = fh.getrname(aln.tid)
                    fhout[hap].write(struct.pack('>I', rid[aln.qname]))
                    fhout[hap].write(struct.pack('>I', lid[locus]))
        for hap in self.hname:
            fhout[hap].close()

    def produce(self, h5file, title='Alignments', index_dtype='uint32', data_dtype=float, complib='zlib', incidence_only=True):
        h5fh = tables.open_file(h5file, 'w', title=title)
        fil  = tables.Filters(complevel=1, complib=complib)
        h5fh.set_node_attr(h5fh.root, 'incidence_only', incidence_only)
        h5fh.set_node_attr(h5fh.root, 'mtype', 'csc_matrix')
        h5fh.set_node_attr(h5fh.root, 'shape', (len(self.lname), len(self.hname), len(self.rname)))
        h5fh.set_node_attr(h5fh.root, 'hname', self.hname)
        h5fh.create_carray(h5fh.root, 'lname', obj=self.lname, title='Locus Names', filters=fil)
        h5fh.create_carray(h5fh.root, 'rname', obj=self.rname, title='Read Names', filters=fil)
        for hid in xrange(len(self.hname)):
            hap = self.hname[hid]
            infile = self.tmpfiles[hap]
            dmat = np.fromfile(open(infile, 'rb'), dtype='>I')
            dmat = dmat.reshape((len(dmat)/2, 2)).T
            if dmat.shape[0] > 2:
                dvec = dmat[2]
            else:
                dvec = np.ones(dmat.shape[1])
            spmat = coo_matrix((dvec, dmat[:2]), shape=(len(self.rname), len(self.lname)))
            spmat = spmat.tocsc()
            hgroup = h5fh.create_group(h5fh.root, 'h%d' % hid, 'Sparse matrix components for Haplotype %d' % hid)
            i1 = h5fh.create_carray(hgroup, 'indptr', obj=spmat.indptr.astype(index_dtype), filters=fil)
            i2 = h5fh.create_carray(hgroup, 'indices', obj=spmat.indices.astype(index_dtype), filters=fil)
            if not incidence_only:
                d = h5fh.create_carray(hgroup, 'data', obj=spmat.data.astype(data_dtype), filters=fil)
        h5fh.flush()
        h5fh.close()

    def cleanup(self):
        for tmpfile in self.tmpfiles.itervalues():
            os.remove(tmpfile)


if __name__ == "__main__":
    pass # TODO: Put some simple test
