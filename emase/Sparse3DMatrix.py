#!/usr/bin/env python

import copy
import tables
import numpy as np
from scipy.sparse import lil_matrix, csc_matrix, csr_matrix, coo_matrix, hstack, vstack
from numbers import Number


class Sparse3DMatrix:
    """
    3-dim sparse matrix designed for "pooled" RNA-seq alignments
    """
    def __init__(self, other=None, h5file=None, datanode='/', shape=None, dtype=float):
        self.shape = (0, 0, 0)
        self.ndim = 3
        self.data = list()
        self.finalized = False # The data matrix has been populated with alignment data and converted to csc format
        # If copying from 'other' Sparse3DMatrix
        if other is not None:
            if other.finalized:
                self.shape = other.shape
                self.data  = copy.deepcopy(other.data)
                self.finalized = True
            else:
                raise RuntimeError('The original matrix must be finalized.')
        # If reading in from emase format file
        elif h5file is not None:
            h5fh = tables.open_file(h5file, 'r')
            self.shape = h5fh.get_node_attr(datanode, 'shape')
            for hid in xrange(self.shape[1]):
                self.data.append(self._reconstruct_spmat(h5fh, hid, datanode, dtype))
            h5fh.close()
            self.finalize() # This also converts the data matrices into csc format
        # If creating from scratch
        elif shape is not None:
            if len(shape) != 3:
                raise RuntimeError('The shape must be a tuple of three positive integers.')
            elif np.any(shape < 1):
                raise RuntimeError('The shape must be a tuple of three positive integers.')
            else:
                self.shape = shape
                for hid in xrange(self.shape[1]):
                    self.data.append(lil_matrix((self.shape[2], self.shape[0]), dtype=dtype))

    def _reconstruct_spmat(self, h5fh, hid, datanode, dtype):
        try:
            mtype = h5fh.get_node_attr(datanode, 'mtype')
            incidence_only = h5fh.get_node_attr(datanode, 'incidence_only')
        except AttributeError:
            mtype = 'coo_matrix'
            incidence_only = False
        hapnode = h5fh.get_node(datanode + '/h%d' % hid)
        if mtype == 'csc_matrix':
            indptr = h5fh.get_node(hapnode, 'indptr').read()
            indptr = indptr.astype(int)
            indices = h5fh.get_node(hapnode, 'indices').read()
            indices = indices.astype(int)
            if not incidence_only:
                data = h5fh.get_node(hapnode, 'data').read()
                data = data.astype(dtype)
            else:
                data = np.ones(len(indices), dtype=dtype)
            spmat = csc_matrix((data, indices, indptr), shape=(self.shape[2], self.shape[0]))
        elif mtype is 'coo_matrix':
            coor = h5fh.get_node(hapnode, 'coor').read()
            data = h5fh.get_node(hapnode, 'data').read()
            data = data.astype(dtype)
            spmat = coo_matrix((data, coor), shape=(self.shape[2], self.shape[0]))
        else:
            raise RuntimeError('Only csc or coo matrices are supported.')
        return spmat

    def copy(self):
        if self.finalized:
            dmat = self.__class__()
            dmat.shape = self.shape
            dmat.data = copy.deepcopy(self.data)
            dmat.finalized = True
            return dmat
        else:
            raise RuntimeError('The original matrix must be finalized.')

    #
    # Binary operators
    #

    def __add__(self, other):
        if self.finalized:
            dmat = self.__class__()
            dmat.shape = self.shape
            if isinstance(other, Sparse3DMatrix):
                if other.finalized:
                    for hid in xrange(self.shape[1]):
                        dmat.data.append(self.data[hid] + other.data[hid])
                else:
                    raise RuntimeError('Both matrices must be finalized.')
            elif isinstance(other, (csc_matrix, csr_matrix, coo_matrix, lil_matrix)):
                other_csc = other.tocsc()
                for hid in xrange(self.shape[1]):
                    dmat.data.append(self.data[hid] + other_csc)
            else:
                raise TypeError('This operator is not supported between the given types.')
            dmat.finalized = True
            return dmat
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def __sub__(self, other):
        if self.finalized:
            dmat = self.__class__()
            dmat.shape = self.shape
            if isinstance(other, Sparse3DMatrix):
                if other.finalized:
                    for hid in xrange(self.shape[1]):
                        dmat.data.append(self.data[hid] - other.data[hid])
                else:
                    raise RuntimeError('Both matrices must be finalized.')
            elif isinstance(other, (csc_matrix, csr_matrix, coo_matrix, lil_matrix)):
                other_csc = other.tocsc()
                for hid in xrange(self.shape[1]):
                    dmat.data.append(self.data[hid] - other_csc)
            else:
                raise TypeError('This operator is not supported between the given types.')
            dmat.finalized = True
            return dmat
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def __mul__(self, other):
        if self.finalized:
            dmat = self.__class__()
            dmat.shape = self.shape
            if isinstance(other, Sparse3DMatrix):  # element-wise multiplication between same kind
                if other.finalized:
                    for hid in xrange(self.shape[1]):
                        dmat.data.append(self.data[hid].multiply(other.data[hid]))
                else:
                    raise RuntimeError('Both matrices must be finalized.')
            elif isinstance(other, (np.ndarray, csc_matrix, csr_matrix)):  # matrix-matrix multiplication
                for hid in xrange(self.shape[1]):
                    dmat.data.append(self.data[hid] * other)
                dmat.shape = (other.shape[1], self.shape[1], self.shape[2])
            elif isinstance(other, (coo_matrix, lil_matrix)):              # matrix-matrix multiplication
                other_csc = other.tocsc()
                for hid in xrange(self.shape[1]):
                    dmat.data.append(self.data[hid] * other_csc)
                dmat.shape = (other_csc.shape[1], self.shape[1], self.shape[2])
            elif isinstance(other, Number):  # rescaling of matrix
                for hid in xrange(self.shape[1]):
                    dmat.data.append(self.data[hid] * other)
            else:
                raise TypeError('This operator is not supported between the given types.')
            dmat.finalized = True
            return dmat
        else:
            raise RuntimeError('The original matrix must be finalized.')

    #
    # Methods for populating the data matrix
    #

    def set_value(self, lid, hid, rid, value):
        if np.all(self.shape > 0):
            self.data[hid][rid, lid] = value
        else:
            raise RuntimeError('The Sparse3DMatrix has only been declared.')

    def add_value(self, lid, hid, rid, value):
        if np.all(self.shape > 0):
            self.data[hid][rid, lid] += value
        else:
            raise RuntimeError('The Sparse3DMatrix has only been declared.')

    def reset(self):
        if self.finalized:
            for hid in xrange(self.shape[1]):
                # Todo: inherit the dtype from orig
                self.data[hid].data = np.ones(self.data[hid].nnz, dtype=self.data[hid].dtype)
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def finalize(self): # Run this when all the populating job is done
        if not self.finalized:
            for hid in xrange(self.shape[1]):
                self.data[hid] = self.data[hid].tocsc()
            self.finalized = True

    #
    # Summarization methods
    #

    def sum(self, axis=2):
        if self.finalized:
            if axis == 0:
                sum_mat = [] # sum along axis=0
                for hid in xrange(self.shape[1]):
                    sum_mat.append(self.data[hid].sum(axis=1).A)
                sum_mat = np.hstack(sum_mat)
            elif axis == 1:  # sum along axis=1
                sum_mat = self.data[0]
                for hid in xrange(1, self.shape[1]):
                    sum_mat = sum_mat + self.data[hid] # Unlike others, this sum_mat is still sparse matrix
            elif axis == 2:  # sum along axis=2
                sum_mat = []
                for hid in xrange(self.shape[1]):
                    sum_mat.append(self.data[hid].sum(axis=0).A)
                sum_mat = np.vstack(sum_mat) # TODO: Fix the dimension by using hstack?
            else:
                raise RuntimeError('The axis should be 0, 1, or 2.')
            return sum_mat
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def get_cross_section(self, index, axis=0):
        if self.finalized:
            if axis == 0:
                cols = []
                for hid in xrange(len(self.data)):
                    cols.append(self.data[hid][:, index])
                return hstack(cols).tocsc()
            elif axis == 1:
                return self.data[axis]
            elif axis == 2:
                rows = []
                for hid in xrange(len(self.data)):
                    rows.append(self.data[hid][index, :])
                return vstack(rows).tocsc()
            else:
                raise RuntimeError('The axis should be 0, 1, or 2.')
        else:
            raise RuntimeError('The original matrix must be finalized.')

    #
    # In-place operations
    #

    def add(self, addend_mat, axis=1):
        """
        In-place addition

        :param addend_mat: A matrix to be added on the Sparse3DMatrix object
        :param axis: The dimension along the addend_mat is added
        :return: Nothing (as it performs in-place operations)
        """
        if self.finalized:
            if axis == 0:
                raise NotImplementedError('The method is not yet implemented for the axis.')
            elif axis == 1:
                for hid in xrange(self.shape[1]):
                    self.data[hid] = self.data[hid] + addend_mat
            elif axis == 2:
                raise NotImplementedError('The method is not yet implemented for the axis.')
            else:
                raise RuntimeError('The axis should be 0, 1, or 2.')
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def multiply(self, multiplier, axis=None):
        """
        In-place multiplication

        :param multiplier: A matrix or vector to be multiplied
        :param axis: The dim along which 'multiplier' is multiplied
        :return: Nothing (as it performs in-place operations)
        """
        if self.finalized:
            if multiplier.ndim == 1:
                if axis == 0:  # multiplier is np.array of length |haplotypes|
                    raise NotImplementedError('The method is not yet implemented for the axis.')
                elif axis == 1:  # multiplier is np.array of length |loci|
                    sz = len(multiplier)
                    multiplier_mat = lil_matrix((sz, sz))
                    multiplier_mat.setdiag(multiplier)
                    for hid in xrange(self.shape[1]):
                        self.data[hid] = self.data[hid] * multiplier_mat
                elif axis == 2:  # multiplier is np.array of length |reads|
                    for hid in xrange(self.shape[1]):
                        self.data[hid].data *= multiplier[self.data[hid].indices]
                else:
                    raise RuntimeError('The axis should be 0, 1, or 2.')
            elif multiplier.ndim == 2:
                if axis == 0:  # multiplier is sp.sparse matrix of shape |reads| x |haplotypes|
                    for hid in xrange(self.shape[1]):
                        self.data[hid].data *= multiplier[self.data[hid].indices, hid]
                elif axis == 1:  # multiplier is sp.sparse matrix of shape |reads| x |loci|
                    for hid in xrange(self.shape[1]):
                        self.data[hid] = self.data[hid].multiply(multiplier)
                elif axis == 2:  # multiplier is np.matrix of shape |haplotypes| x |loci|
                    for hid in xrange(self.shape[1]):
                        multiplier_vec = multiplier[hid, :]
                        multiplier_vec = multiplier_vec.ravel()
                        self.data[hid].data *= multiplier_vec.repeat(np.diff(self.data[hid].indptr))
                else:
                    raise RuntimeError('The axis should be 0, 1, or 2.')
            elif isinstance(multiplier, Sparse3DMatrix):  # multiplier is Sparse3DMatrix object
                    for hid in xrange(self.shape[1]):
                        self.data[hid] = self.data[hid].multiply(multiplier.data[hid])
            else:
                raise RuntimeError('The multiplier should be 1, 2 dimensional numpy array or a Sparse3DMatrix object.')
        else:
            raise RuntimeError('The original matrix must be finalized.')

    #
    # Helper methods
    #

    def combine(self, other):
        if self.finalized and other.finalized:
            dmat = self.__class__()
            dmat.shape = (self.shape[0], self.shape[1], self.shape[2] + other.shape[2])
            for hid in xrange(self.shape[1]):
                dmat.data.append(vstack([self.data[hid], other.data[hid]])) # WARNING: vstack returns coo_matrix!!
            dmat.finalize()
            return dmat
        else:
            raise RuntimeError('Both matrices must be finalized.')

    def save(self, h5file, title=None, index_dtype='uint32', data_dtype=float, incidence_only=True, complib='zlib'):
        if self.finalized:
            h5fh = tables.open_file(h5file, 'w', title=title)
            fil  = tables.Filters(complevel=1, complib=complib)
            h5fh.set_node_attr(h5fh.root, 'incidence_only', incidence_only)
            h5fh.set_node_attr(h5fh.root, 'mtype', 'csc_matrix')
            h5fh.set_node_attr(h5fh.root, 'shape', self.shape)
            for hid in xrange(self.shape[1]):
                hgroup = h5fh.create_group(h5fh.root, 'h%d' % hid, 'Sparse matrix components for Haplotype %d' % hid)
                spmat = self.data[hid]
                i1 = h5fh.create_carray(hgroup, 'indptr', obj=spmat.indptr.astype(index_dtype), filters=fil)
                i2 = h5fh.create_carray(hgroup, 'indices', obj=spmat.indices.astype(index_dtype), filters=fil)
                if not incidence_only:
                    d = h5fh.create_carray(hgroup, 'data', obj=spmat.data.astype(data_dtype), filters=fil)
            h5fh.flush()
            h5fh.close()
        else:
            raise RuntimeError('The matrix is not finalized.')


if __name__ == "__main__":
    pass  # TODO: Give some usage example

