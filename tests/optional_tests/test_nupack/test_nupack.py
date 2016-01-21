import coral as cr
import csv
import os
from nose.tools import assert_equal, assert_true
import numpy as np


class TestNUPACK(object):
    def __init__(self):
        self.dnas = [cr.DNA('GATACTAGCG'),
                     cr.DNA('TACGATT'),
                     cr.DNA('GATACG')]
        self.rnas = [s.transcribe() for s in self.dnas]
        self.nupack = cr.analysis.NUPACK()

    def test_pfunc(self):
        '''Test the simplest (partition function) command pfunc with a single
        sequence input.'''
        # test DNA
        output_dna = self.nupack.pfunc(self.dnas[0])
        assert_equal(output_dna['free_energy'], -4.92069506e-02)
        assert_equal(output_dna['pfunc'], 1.08311357973974e+00)
        # test RNA with 1995 params
        output_rna = self.nupack.pfunc(self.rnas[0])
        assert_equal(output_rna['free_energy'], -9.28516187e-02)
        assert_equal(output_rna['pfunc'], 1.16259513934557e+00)
        # test RNA with 1999 params
        output_rna99 = self.nupack.pfunc(self.rnas[0], material='rna1999')
        assert_equal(output_rna99['free_energy'], -7.97413222e-03)
        assert_equal(output_rna99['pfunc'], 1.01302234408117e+00)

    def test_pfunc_multi(self):
        '''Test the simplest (partition function) command pfunc with the
        -multi option, which returns an ordered complex partition function.'''
        # test DNA
        output_dna = self.nupack.pfunc_multi(self.dnas)
        assert_equal(output_dna['free_energy'], -9.59943928e+00)
        assert_equal(output_dna['pfunc'], 5.81176347268940e+06)
        # test RNA with 1995 params
        output_rna = self.nupack.pfunc_multi(self.rnas)
        assert_equal(output_rna['free_energy'], -5.45632785e+00)
        assert_equal(output_rna['pfunc'], 6.99579788609535e+03)
        # test RNA with 1999 params
        output_rna99 = self.nupack.pfunc_multi(self.rnas, material='rna1999')
        assert_equal(output_rna99['free_energy'], -5.27740504e+00)
        assert_equal(output_rna99['pfunc'], 5.23308895793574e+03)

    def test_pairs(self):
        '''Test the pairs command.'''
        # test DNA
        dna_mat = self._process_ppairs('pairs_dna.tsv', len(self.dnas[0]))
        output_dna = self.nupack.pairs(self.dnas[0])
        assert_true(np.array_equal(dna_mat, output_dna))

        # test RNA
        rna_mat = self._process_ppairs('pairs_rna.tsv', len(self.rnas[0]))
        output_rna = self.nupack.pairs(self.rnas[0])
        assert_true(np.array_equal(rna_mat, output_rna))

        # test RNA 1999
        rna99_mat = self._process_ppairs('pairs_rna99.tsv', len(self.rnas[0]))
        output_rna99 = self.nupack.pairs(self.rnas[0], material='rna1999')
        assert_true(np.array_equal(rna99_mat, output_rna99))

    def test_pairs_multi(self):
        '''Test the pairs command with the -multi flag.'''
        # Test DNA
        dim = sum([len(s) for s in self.dnas])
        dna_ppairs = self._process_ppairs('pairs_multi_dna.ppairs', dim)
        dna_epairs = self._process_ppairs('pairs_multi_dna.epairs', dim)
        dna_output = self.nupack.pairs_multi(self.dnas)
        for expected, output in zip([dna_ppairs, dna_epairs], dna_output):
            assert_true(np.array_equal(expected, output))

        # Test RNA
        dim = sum([len(s) for s in self.rnas])
        rna_ppairs = self._process_ppairs('pairs_multi_rna.ppairs', dim)
        rna_epairs = self._process_ppairs('pairs_multi_rna.epairs', dim)
        rna_output = self.nupack.pairs_multi(self.rnas)
        for expected, output in zip([rna_ppairs, rna_epairs], rna_output):
            assert_true(np.array_equal(expected, output))

        # Test RNA 1999
        dim = sum([len(s) for s in self.rnas])
        rna99_ppairs = self._process_ppairs('pairs_multi_rna99.ppairs', dim)
        rna99_epairs = self._process_ppairs('pairs_multi_rna99.epairs', dim)
        expected_mats = [rna99_ppairs, rna99_epairs]
        rna99_output = self.nupack.pairs_multi(self.rnas, material='rna1999')
        for expected, output in zip(expected_mats, rna99_output):
            assert_true(np.array_equal(expected, output))

    def _process_ppairs(self, filename, dim):
        mat = np.zeros((dim, dim + 1))
        curdir = os.path.dirname(__file__)
        tsvpath = os.path.join(curdir, 'data', filename)
        with open(tsvpath) as f:
            reader = csv.reader(f, delimiter='\t')
            for line in reader:
                i = int(line[0]) - 1
                j = int(line[1]) - 1
                prob = float(line[2])
                mat[i, j] = prob
        return mat
