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
        assert_equal(output_dna, (-4.92069506e-02, 1.08311357973974e+00))
        # test RNA with 1995 params
        output_rna = self.nupack.pfunc(self.rnas[0])
        assert_equal(output_rna, (-9.28516187e-02, 1.16259513934557e+00))
        # test RNA with 1999 params
        output_rna99 = self.nupack.pfunc(self.rnas[0], material='rna1999')
        assert_equal(output_rna99, (-7.97413222e-03, 1.01302234408117e+00))

    def test_pfunc_multi(self):
        '''Test the simplest (partition function) command pfunc with the
        -multi option, which returns an ordered complex partition function.'''
        # test DNA
        output_dna = self.nupack.pfunc_multi(self.dnas)
        assert_equal(output_dna, (-9.59943928e+00, 5.81176347268940e+06))
        # test RNA with 1995 params
        output_rna = self.nupack.pfunc_multi(self.rnas)
        assert_equal(output_rna, (-5.45632785e+00, 6.99579788609535e+03))
        # test RNA with 1999 params
        output_rna99 = self.nupack.pfunc_multi(self.rnas, material='rna1999')
        assert_equal(output_rna99, (-5.27740504e+00, 5.23308895793574e+03))

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

    def test_mfe(self):
        # Test DNA
        dna_output = self.nupack.mfe(sum(self.dnas))
        assert_equal(dna_output['mfe'], -1.210)
        assert_equal(dna_output['dotparens'], '........((((.......))))')
        assert_equal(dna_output['pairlist'],
                     [[8, 22], [9, 21], [10, 20], [11, 19]])

        # Test RNA
        rna_output = self.nupack.mfe(sum(self.rnas))
        assert_equal(rna_output['mfe'], -1.100)
        assert_equal(rna_output['dotparens'], '........((((.......))))')
        assert_equal(rna_output['pairlist'],
                     [[8, 22], [9, 21], [10, 20], [11, 19]])

        # Test RNA 1999
        rna99_output = self.nupack.mfe(sum(self.rnas), material='rna1999')
        assert_equal(rna99_output['mfe'], -0.300)
        assert_equal(rna99_output['dotparens'], '........((((.......))))')
        assert_equal(rna99_output['pairlist'],
                     [[8, 22], [9, 21], [10, 20], [11, 19]])

    def test_mfe_degenerate(self):
        # Test '-degenerate' flag with DNA
        degenerate_input = self.dnas[0] + self.dnas[0] + self.dnas[0]
        degenerate_output = self.nupack.mfe(degenerate_input, degenerate=True)
        # Should generate 2 degenerate equal-MFE structures
        assert_equal(degenerate_output[0]['mfe'], -1.330)
        assert_equal(degenerate_output[0]['dotparens'],
                     '..............((((......))))..')
        assert_equal(degenerate_output[0]['pairlist'],
                     [[14, 27], [15, 26], [16, 25], [17, 24]])
        assert_equal(degenerate_output[1]['mfe'], -1.330)
        assert_equal(degenerate_output[1]['dotparens'],
                     '....((((......))))............')
        assert_equal(degenerate_output[1]['pairlist'],
                     [[4, 17], [5, 16], [6, 15], [7, 14]])

    def test_mfe_multi(self):
        # Test DNA
        dna_output = self.nupack.mfe_multi(self.dnas)
        assert_equal(dna_output['mfe'], -8.773)
        assert_equal(dna_output['dotparens'], '.((.....((+..))...+.))...')
        assert_equal(dna_output['pairlist'],
                     [[1, 19], [2, 18], [8, 13], [9, 12]])

        # Test RNA
        rna_output = self.nupack.mfe_multi(self.rnas)
        assert_equal(rna_output['mfe'], -3.863)
        assert_equal(rna_output['dotparens'], '(.......((+..))...+....).')
        assert_equal(rna_output['pairlist'],
                     [[0, 21], [8, 13], [9, 12]])

        # Test RNA 1999
        rna99_output = self.nupack.mfe_multi(self.rnas, material='rna1999')
        assert_equal(rna99_output['mfe'], -4.263)
        assert_equal(rna99_output['dotparens'], '(.......((+..))...+....).')
        assert_equal(rna99_output['pairlist'],
                     [[0, 21], [8, 13], [9, 12]])

    def test_subopt(self):
        # Test DNA
        dna_output = self.nupack.subopt(self.dnas[0], 2.5)
        # For DNA, 3 are found
        assert_equal(dna_output[0]['mfe'], 0.000)
        assert_equal(dna_output[0]['dotparens'], '..........')
        assert_equal(dna_output[0]['pairlist'], [])
        assert_equal(dna_output[1]['mfe'], 1.940)
        assert_equal(dna_output[1]['dotparens'], '....(....)')
        assert_equal(dna_output[1]['pairlist'], [[4, 9]])
        assert_equal(dna_output[2]['mfe'], 2.500)
        assert_equal(dna_output[2]['dotparens'], '.(...)....')
        assert_equal(dna_output[2]['pairlist'], [[1, 5]])

        # Test RNA
        rna_output = self.nupack.subopt(self.rnas[0], 2.5)
        assert_equal(rna_output[0]['mfe'], 0.000)
        assert_equal(rna_output[0]['dotparens'], '..........')
        assert_equal(rna_output[0]['pairlist'], [])
        assert_equal(rna_output[1]['mfe'], 1.300)
        assert_equal(rna_output[1]['dotparens'], '(.......).')
        assert_equal(rna_output[1]['pairlist'], [[0, 8]])

    def test_subopt_multi(self):
        # Test DNA
        dna_output = self.nupack.subopt_multi(self.dnas, 0.5)
        assert_equal(dna_output[0]['mfe'], -8.773)
        assert_equal(dna_output[0]['dotparens'], '.((.....((+..))...+.))...')
        assert_equal(dna_output[0]['pairlist'],
                     [[1, 19], [2, 18], [8, 13], [9, 12]])
        assert_equal(dna_output[1]['mfe'], -8.323)
        assert_equal(dna_output[1]['dotparens'], '.((...(.((+..)).).+.))...')
        assert_equal(dna_output[1]['pairlist'],
                     [[1, 19], [2, 18], [6, 15], [8, 13], [9, 12]])

        # Test RNA
        rna_output = self.nupack.subopt_multi(self.rnas, 0.5)
        assert_equal(rna_output[0]['mfe'], -3.863)
        assert_equal(rna_output[0]['dotparens'], '(.......((+..))...+....).')
        assert_equal(rna_output[0]['pairlist'], [[0, 21], [8, 13], [9, 12]])
        assert_equal(rna_output[1]['mfe'], -3.663)
        assert_equal(rna_output[1]['dotparens'], '.((.....((+..))...+.))...')
        assert_equal(rna_output[1]['pairlist'],
                     [[1, 19], [2, 18], [8, 13], [9, 12]])

        # Test RNA 1999
        rna99_output = self.nupack.subopt_multi(self.rnas, 0.5,
                                                material='rna1999')
        assert_equal(rna99_output[0]['mfe'], -4.263)
        assert_equal(rna99_output[0]['dotparens'],
                     '(.......((+..))...+....).')
        assert_equal(rna99_output[0]['pairlist'], [[0, 21], [8, 13], [9, 12]])

    def test_count(self):
        # Test DNA
        dna_output = self.nupack.count(self.dnas[0])
        assert_equal(dna_output, 9)

        # Test RNA
        rna_output = self.nupack.count(self.rnas[0])
        assert_equal(rna_output, 9)

        # Test DNA
        rna99_output = self.nupack.count(self.rnas[0])
        assert_equal(rna99_output, 9)

    def test_count_multi(self):
        # Test DNA
        dna_output = self.nupack.count_multi(self.dnas)
        assert_equal(dna_output, 7681)

        # Test RNA
        rna_output = self.nupack.count_multi(self.rnas)
        assert_equal(rna_output, 7681)

        # Test RNA 1999
        rna99_output = self.nupack.count_multi(self.rnas)
        assert_equal(rna99_output, 7681)

    def test_energy(self):
        # Test DNA
        dna_output = self.nupack.energy(self.dnas[0], '..(....)..')
        assert_equal(dna_output, 200003.05000000000000)

        # Test RNA
        rna_output = self.nupack.energy(self.rnas[0], '..(....)..')
        assert_equal(rna_output, 200003.70000000000000)

        # Test RNA 1999
        rna99_output = self.nupack.energy(self.rnas[0], '..(....)..',
                                          material='rna1999')
        assert_equal(rna99_output, 200005.59999999999999)

    def test_energy_multi(self):
        # Test DNA
        dna_output = self.nupack.energy_multi(self.dnas,
                                              '(......(((+..))..)+....).')
        assert_equal(dna_output, 199998.83729716184106)

        # Test RNA
        rna_output = self.nupack.energy_multi(self.rnas,
                                              '(......(((+..))..)+....).')
        assert_equal(rna_output, 200002.13729716184108)

        # Test RNA 1999
        rna99_output = self.nupack.energy_multi(self.rnas,
                                                '(......(((+..))..)+....).',
                                                material='rna1999')
        assert_equal(rna99_output, 200002.43729716184109)

    def test_prob(self):
        # Test DNA
        dna_output = self.nupack.prob(self.dnas[0], '..........')
        assert_equal(dna_output, .9233)

        # Test RNA
        rna_output = self.nupack.prob(self.rnas[0], '..........')
        assert_equal(rna_output, .8601)

        # Test RNA 1999
        rna99_output = self.nupack.prob(self.rnas[0], '..........',
                                        material='rna1999')
        assert_equal(rna99_output, .9871)

    def test_prob_multi(self):
        # Test DNA
        dna_output = self.nupack.prob_multi(self.dnas,
                                            '(.......((+..))...+....).')
        assert_equal(dna_output, .04460)

        # Test RNA
        rna_output = self.nupack.prob_multi(self.rnas,
                                            '(.......((+..))...+....).')
        assert_equal(rna_output, .07534)

        # Test RNA 1999
        rna99_output = self.nupack.prob_multi(self.rnas,
                                              '(.......((+..))...+....).',
                                              material='rna1999')
        assert_equal(rna99_output, .1927)

    def test_defect(self):
        # Test DNA
        dna_output = self.nupack.defect(self.dnas[0], '..(....)..')
        assert_equal(dna_output, [2.152, .2152])

        # Test RNA
        rna_output = self.nupack.defect(self.rnas[0], '..(....)..')
        assert_equal(rna_output, [2.274, .2274])

        # Test RNA 1999
        rna99_output = self.nupack.defect(self.rnas[0], '..(....)..',
                                          material='rna1999')
        assert_equal(rna99_output, [2.025, .2025])

    def test_defect_multi(self):
        # Test DNA
        dna_output = self.nupack.defect_multi(self.dnas,
                                              '(.......((+..))...+....).')
        assert_equal(dna_output, [5.790, .2517])

        # Test RNA
        rna_output = self.nupack.defect_multi(self.rnas,
                                              '(.......((+..))...+....).')
        assert_equal(rna_output, [6.522, .2836])

        # Test RNA 1999
        rna99_output = self.nupack.defect_multi(self.rnas,
                                                '(.......((+..))...+....).',
                                                material='rna1999')
        assert_equal(rna99_output, [4.733, .2058])

    def test_complexes(self):
        dnas = self.dnas[:2]
        rnas = self.rnas[:2]
        # Test DNA
        dna_data = [[1, 0, -4.92069506e-02],
                    [0, 1, -2.01874269e-02],
                    [2, 0, -5.76962013e+00],
                    [1, 1, -5.60626976e+00],
                    [0, 2, -6.02915496e+00],
                    [3, 0, -1.02261804e+01],
                    [2, 1, -1.08277874e+01],
                    [1, 2, -9.83335621e+00],
                    [0, 3, -8.94565103e+00],
                    [4, 0, -1.54019941e+01],
                    [3, 1, -1.57681947e+01],
                    [2, 2, -1.56221971e+01],
                    [1, 3, -1.48468888e+01],
                    [0, 4, -1.41217372e+01]]
        dna_expected = [{'complex': [int(row[0]), int(row[1])],
                         'energy': float(row[2])} for row in dna_data]

        dna_output = self.nupack.complexes(dnas, 4)
        for expect, output in zip(dna_expected, dna_output):
            assert_equal(expect['complex'], output['complex'])
            assert_equal(expect['energy'], output['energy'])

        # Test RNA
        rna_data = [[1, 0, -9.28516187e-02],
                    [0, 1, -4.32138317e-03],
                    [2, 0, -5.33244849e+00],
                    [1, 1, -2.96200638e+00],
                    [0, 2, -3.28345221e+00],
                    [3, 0, -8.66644413e+00],
                    [2, 1, -8.12678954e+00],
                    [1, 2, -5.50808307e+00],
                    [0, 3, -4.51692546e+00],
                    [4, 0, -1.33713694e+01],
                    [3, 1, -1.17159283e+01],
                    [2, 2, -1.08018186e+01],
                    [1, 3, -7.97713466e+00],
                    [0, 4, -6.98583185e+00]]

        rna_expected = [{'complex': [int(row[0]), int(row[1])],
                         'energy': float(row[2])} for row in rna_data]

        rna_output = self.nupack.complexes(rnas, 4)
        for expect, output in zip(rna_expected, rna_output):
            assert_equal(expect['complex'], output['complex'])
            assert_equal(expect['energy'], output['energy'])

        # Test RNA 1999
        rna99_data = [[1, 0, -7.97413222e-03],
                      [0, 1, -1.09376608e-04],
                      [2, 0, -5.66294745e+00],
                      [1, 1, -3.00357694e+00],
                      [0, 2, -3.61831003e+00],
                      [3, 0, -8.58329727e+00],
                      [2, 1, -8.40127714e+00],
                      [1, 2, -5.24448000e+00],
                      [0, 3, -3.89211467e+00],
                      [4, 0, -1.34538247e+01],
                      [3, 1, -1.15873619e+01],
                      [2, 2, -1.07578288e+01],
                      [1, 3, -7.43779859e+00],
                      [0, 4, -6.56448976e+00]]

        rna99_expected = [{'complex': [int(row[0]), int(row[1])],
                           'energy': float(row[2])} for row in rna99_data]

        rna99_output = self.nupack.complexes(rnas, 4, material='rna1999')
        for expect, output in zip(rna99_expected, rna99_output):
            assert_equal(expect['complex'], output['complex'])
            assert_equal(expect['energy'], output['energy'])

        # Test DNA with pairs option
        dim = sum([len(x) for x in dnas])
        dnapairs_data = self._read_cxepairs('complexes_pairs_dna.cx-epairs')
        for i, pairlist in enumerate(dnapairs_data):
            dna_expected[i]['epairs'] = self._pairs_to_np(pairlist, dim)
        dnapairs_output = self.nupack.complexes(dnas, 4, pairs=True)
        # Since there's a numpy matrix in there, have to use numpy comparison
        # on each element
        for expected, output in zip(dna_expected, dnapairs_output):
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['energy'], output['energy'])
            assert_true(np.array_equal(expected['epairs'], output['epairs']))

        # Test RNA with pairs option
        dim = sum([len(x) for x in rnas])
        rnapairs_data = self._read_cxepairs('complexes_pairs_rna.cx-epairs')
        for i, pairlist in enumerate(rnapairs_data):
            rna_expected[i]['epairs'] = self._pairs_to_np(pairlist, dim)
        rnapairs_output = self.nupack.complexes(rnas, 4, pairs=True)
        # Since there's a numpy matrix in there, have to use numpy comparison
        # on each element
        for expected, output in zip(rna_expected, rnapairs_output):
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['energy'], output['energy'])
            assert_true(np.array_equal(expected['epairs'], output['epairs']))

        # Test RNA 1999 with pairs option
        dim = sum([len(x) for x in rnas])
        rna99pairs_dat = self._read_cxepairs('complexes_pairs_rna99.cx-epairs')
        for i, pairlist in enumerate(rna99pairs_dat):
            rna99_expected[i]['epairs'] = self._pairs_to_np(pairlist, dim)
        rna99pairs_output = self.nupack.complexes(rnas, 4, pairs=True,
                                                  material='rna1999')
        # Since there's a numpy matrix in there, have to use numpy comparison
        # on each element
        for expected, output in zip(rna99_expected, rna99pairs_output):
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['energy'], output['energy'])
            assert_true(np.array_equal(expected['epairs'], output['epairs']))

        # Test DNA with the ordered option
        dna_ocx = [[1, 0, -4.92069506e-02],
                   [0, 1, -2.01874269e-02],
                   [2, 0, -5.76962013e+00],
                   [1, 1, -5.60626976e+00],
                   [0, 2, -6.02915496e+00],
                   [3, 0, -1.02261804e+01],
                   [2, 1, -1.08277874e+01],
                   [1, 2, -9.83335621e+00],
                   [0, 3, -8.94565103e+00],
                   [4, 0, -1.54019941e+01],
                   [3, 1, -1.57681947e+01],
                   [2, 2, -1.51522348e+01],
                   [2, 2, -1.52349820e+01],
                   [1, 3, -1.48468888e+01],
                   [0, 4, -1.41217372e+01]]
        dna_ocx_expect = [{'complex': [l[0], l[1]], 'energy': l[2]} for l in
                          dna_ocx]
        dna_ocx_keys = [[1],
                        [2],
                        [1, 1],
                        [1, 2],
                        [2, 2],
                        [1, 1, 1],
                        [1, 1, 2],
                        [1, 2, 2],
                        [2, 2, 2],
                        [1, 1, 1, 1],
                        [1, 1, 1, 2],
                        [1, 1, 2, 2],
                        [1, 2, 1, 2],
                        [1, 2, 2, 2],
                        [2, 2, 2, 2]]
        for i, key in enumerate(dna_ocx_keys):
            dna_ocx_expect[i]['order'] = key
        dna_ocx = self.nupack.complexes(dnas, 4, ordered=True)
        for expect, output in zip(dna_ocx_expect, dna_ocx):
            assert_equal(expect['complex'], output['complex'])
            assert_equal(expect['energy'], output['energy'])
            assert_equal(expect['order'], output['order'])

        # Test RNA with the ordered option
        rna_ocx = [[1, 0, -9.28516187e-02],
                   [0, 1, -4.32138317e-03],
                   [2, 0, -5.33244849e+00],
                   [1, 1, -2.96200638e+00],
                   [0, 2, -3.28345221e+00],
                   [3, 0, -8.66644413e+00],
                   [2, 1, -8.12678954e+00],
                   [1, 2, -5.50808307e+00],
                   [0, 3, -4.51692546e+00],
                   [4, 0, -1.33713694e+01],
                   [3, 1, -1.17159283e+01],
                   [2, 2, -1.05784008e+01],
                   [2, 2, -1.00680851e+01],
                   [1, 3, -7.97713466e+00],
                   [0, 4, -6.98583185e+00]]
        rna_ocx_expect = [{'complex': [l[0], l[1]], 'energy': l[2]} for l in
                          rna_ocx]
        rna_ocx_keys = [[1],
                        [2],
                        [1, 1],
                        [1, 2],
                        [2, 2],
                        [1, 1, 1],
                        [1, 1, 2],
                        [1, 2, 2],
                        [2, 2, 2],
                        [1, 1, 1, 1],
                        [1, 1, 1, 2],
                        [1, 1, 2, 2],
                        [1, 2, 1, 2],
                        [1, 2, 2, 2],
                        [2, 2, 2, 2]]
        for i, key in enumerate(rna_ocx_keys):
            rna_ocx_expect[i]['order'] = key
        rna_ocx = self.nupack.complexes(rnas, 4, ordered=True)
        for expect, output in zip(rna_ocx_expect, rna_ocx):
            assert_equal(expect['complex'], output['complex'])
            assert_equal(expect['energy'], output['energy'])
            assert_equal(expect['order'], output['order'])

        # Test RNA 99 with the ordered option
        rna99_ocx = [[1, 0, -7.97413222e-03],
                     [0, 1, -1.09376608e-04],
                     [2, 0, -5.66294745e+00],
                     [1, 1, -3.00357694e+00],
                     [0, 2, -3.61831003e+00],
                     [3, 0, -8.58329727e+00],
                     [2, 1, -8.40127714e+00],
                     [1, 2, -5.24448000e+00],
                     [0, 3, -3.89211467e+00],
                     [4, 0, -1.34538247e+01],
                     [3, 1, -1.15873619e+01],
                     [2, 2, -1.03498891e+01],
                     [2, 2, -1.03107451e+01],
                     [1, 3, -7.43779859e+00],
                     [0, 4, -6.56448976e+00]]
        rna99_ocx_expect = [{'complex': [l[0], l[1]], 'energy': l[2]} for l
                            in rna99_ocx]
        rna99_ocx_keys = [[1],
                          [2],
                          [1, 1],
                          [1, 2],
                          [2, 2],
                          [1, 1, 1],
                          [1, 1, 2],
                          [1, 2, 2],
                          [2, 2, 2],
                          [1, 1, 1, 1],
                          [1, 1, 1, 2],
                          [1, 1, 2, 2],
                          [1, 2, 1, 2],
                          [1, 2, 2, 2],
                          [2, 2, 2, 2]]
        for i, key in enumerate(rna99_ocx_keys):
            rna99_ocx_expect[i]['order'] = key
        rna99_ocx = self.nupack.complexes(rnas, 4, ordered=True,
                                          material='rna1999')
        for expect, output in zip(rna99_ocx_expect, rna99_ocx):
            assert_equal(expect['complex'], output['complex'])
            assert_equal(expect['energy'], output['energy'])
            assert_equal(expect['order'], output['order'])

        # Test DNA with the ordered and pairs options
        dnapairs_ocx_d = self._read_cxepairs('complexes_pairs_dna.ocx-epairs')
        for i, pairlist in enumerate(dnapairs_ocx_d):
            dna_ocx_expect[i]['epairs'] = self._pairs_to_np(pairlist, dim)
        dna_ocx = self.nupack.complexes(dnas, 4, ordered=True, pairs=True)
        for expected, output in zip(dna_ocx_expect, dna_ocx):
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['energy'], output['energy'])
            assert_equal(expected['order'], output['order'])
            assert_true(np.array_equal(expected['epairs'], output['epairs']))

        # Test RNA with the ordered and pairs options
        rnapairs_ocx_d = self._read_cxepairs('complexes_pairs_rna.ocx-epairs')
        for i, pairlist in enumerate(rnapairs_ocx_d):
            rna_ocx_expect[i]['epairs'] = self._pairs_to_np(pairlist, dim)
        rna_ocx = self.nupack.complexes(rnas, 4, ordered=True, pairs=True)
        for expected, output in zip(rna_ocx_expect, rna_ocx):
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['energy'], output['energy'])
            assert_equal(expected['order'], output['order'])
            assert_true(np.array_equal(expected['epairs'], output['epairs']))

        # Test RNA 1999 with the ordered and pairs options
        r99pairs_ocx_epairs_file = 'complexes_pairs_rna99.ocx-epairs'
        rna99pairs_ocx_d = self._read_cxepairs(r99pairs_ocx_epairs_file)
        for i, pairlist in enumerate(rna99pairs_ocx_d):
            rna99_ocx_expect[i]['epairs'] = self._pairs_to_np(pairlist, dim)
        rna99_ocx = self.nupack.complexes(rnas, 4, ordered=True, pairs=True,
                                          material='rna1999')
        for expected, output in zip(rna99_ocx_expect, rna99_ocx):
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['energy'], output['energy'])
            assert_equal(expected['order'], output['order'])
            assert_true(np.array_equal(expected['epairs'], output['epairs']))

        # Test DNA with the mfe option
        dna_ocx_mfe_expect = self._process_mfe('complexes_mfe_dna.ocx-mfe')
        for expect, mfedat in zip(dna_ocx_expect, dna_ocx_mfe_expect):
            expect['mfe'] = mfedat['mfe']
            expect['dotparens'] = mfedat['dotparens']
            expect['pairlist'] = mfedat['pairlist']
        dna_ocx_mfe = self.nupack.complexes(dnas, 4, mfe=True)
        for expected, output in zip(dna_ocx_expect, dna_ocx_mfe):
            assert_equal(expected['energy'], output['energy'])
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['order'], output['order'])
            assert_equal(expected['mfe'], output['mfe'])
            assert_equal(expected['dotparens'], output['dotparens'])
            assert_equal(expected['pairlist'], output['pairlist'])

        # Test RNA with the mfe option
        rna_ocx_mfe_expect = self._process_mfe('complexes_mfe_rna.ocx-mfe')
        for expect, mfedat in zip(rna_ocx_expect, rna_ocx_mfe_expect):
            expect['mfe'] = mfedat['mfe']
            expect['dotparens'] = mfedat['dotparens']
            expect['pairlist'] = mfedat['pairlist']
        rna_ocx_mfe = self.nupack.complexes(rnas, 4, mfe=True)
        for expected, output in zip(rna_ocx_expect, rna_ocx_mfe):
            assert_equal(expected['energy'], output['energy'])
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['order'], output['order'])
            assert_equal(expected['mfe'], output['mfe'])
            assert_equal(expected['dotparens'], output['dotparens'])
            assert_equal(expected['pairlist'], output['pairlist'])

        # Test RNA 1999 with the mfe option
        rna99_ocx_mfe_expect = self._process_mfe('complexes_mfe_rna99.ocx-mfe')
        for expect, mfedat in zip(rna99_ocx_expect, rna99_ocx_mfe_expect):
            expect['mfe'] = mfedat['mfe']
            expect['dotparens'] = mfedat['dotparens']
            expect['pairlist'] = mfedat['pairlist']
        rna99_ocx_mfe = self.nupack.complexes(rnas, 4, mfe=True,
                                              material='rna1999')
        for expected, output in zip(rna99_ocx_expect, rna99_ocx_mfe):
            assert_equal(expected['energy'], output['energy'])
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['order'], output['order'])
            assert_equal(expected['mfe'], output['mfe'])
            assert_equal(expected['dotparens'], output['dotparens'])
            assert_equal(expected['pairlist'], output['pairlist'])

        # Test DNA with the mfe and pairs options
        dna_ocx_mfe_pairs = self.nupack.complexes(dnas, 4, mfe=True,
                                                  pairs=True)
        for expected, output in zip(dna_ocx_expect, dna_ocx_mfe_pairs):
            assert_equal(expected['energy'], output['energy'])
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['order'], output['order'])
            assert_equal(expected['mfe'], output['mfe'])
            assert_equal(expected['dotparens'], output['dotparens'])
            assert_equal(expected['pairlist'], output['pairlist'])
            assert_true(np.array_equal(expected['epairs'], output['epairs']))

        # TODO: Restore the 'RNA' material version - for whatever reason, it
        # has a disagreeing dotparens structure in it.

        # Test RNA 1999 with the mfe and pairs options
        rna99_ocx_mfe_pairs = self.nupack.complexes(rnas, 4, mfe=True,
                                                    pairs=True,
                                                    material='rna1999')
        for expected, output in zip(rna99_ocx_expect, rna99_ocx_mfe_pairs):
            assert_equal(expected['energy'], output['energy'])
            assert_equal(expected['complex'], output['complex'])
            assert_equal(expected['order'], output['order'])
            assert_equal(expected['mfe'], output['mfe'])
            assert_equal(expected['dotparens'], output['dotparens'])
            assert_equal(expected['pairlist'], output['pairlist'])
            assert_true(np.array_equal(expected['epairs'], output['epairs']))

    def test_complexes_timeonly(self):
        # Test complex size of 4
        dna_4 = self.nupack.complexes_timeonly(self.dnas[:2], 4)
        assert_equal(dna_4, 0.33)
        rna_4 = self.nupack.complexes_timeonly(self.rnas[:2], 4)
        assert_equal(rna_4, 0.33)
        rna99_4 = self.nupack.complexes_timeonly(self.dnas[:2], 4)
        assert_equal(rna99_4, 0.33)

        # Test complex size of 8
        dna_8 = self.nupack.complexes_timeonly(self.dnas[:2], 8)
        assert_equal(dna_8, 18.66)
        rna_8 = self.nupack.complexes_timeonly(self.rnas[:2], 8)
        assert_equal(rna_8, 18.66)
        rna99_8 = self.nupack.complexes_timeonly(self.dnas[:2], 8)
        assert_equal(rna99_8, 18.66)

    def _test_concentrations(self):
        dnas = self.dnas[:2]

        # Test simple case
        expected = [[9.999907e-07],
                    [9.998446e-07],
                    [7.675362e-11],
                    [3.682595e-12],
                    [1.893207e-12],
                    [2.675905e-16],
                    [1.287201e-16],
                    [2.122290e-18],
                    [4.907111e-19],
                    [8.627347e-21],
                    [6.789014e-22],
                    [1.778772e-22],
                    [2.099592e-24],
                    [4.853486e-25]]

        complexes = self.nupack.complexes(dnas, 4)
        arg_output = self.nupack.concentrations(complexes, 1e-6)
        for expect, out in zip(expected, arg_output):
            assert_equal(expect, out['concentration'])

        # Test ordered case
        ordered_expected = [[9.999907e-07],
                            [9.998446e-07],
                            [7.675362e-11],
                            [3.682595e-12],
                            [1.893207e-12],
                            [2.675905e-16],
                            [1.287201e-16],
                            [2.122290e-18],
                            [4.907111e-19],
                            [8.627347e-21],
                            [6.789014e-22],
                            [1.237905e-22],
                            [5.408664e-23],
                            [2.099592e-24],
                            [4.853486e-25]]

        ocomplexes = self.nupack.complexes(dnas, 4)
        ordered_output = self.nupack.concentrations(ocomplexes, 1e-6,
                                                    ordered=True)
        for expect, out in zip(ordered_expected, ordered_output):
            assert_equal(expect, out['concentration'])

        # Test with pairs
        pairs_expected = self._read_pairs('concentrations_pairs.cx-epairs')
        pairs_output = self.nupack.concentrations(dnas, 1e-6, pairs=True)
        for expect, out in zip(pairs_expected, pairs_output):
            assert_equal(expect, out['concentration'])

        # Test with pairs with complexes as argument
        arg_pairs_output = self.nupack.concentrations(dnas, 1e-6, pairs=True,
                                                      complexes=complexes)
        for expect, out in zip(pairs_expected, arg_pairs_output):
            assert_equal(expect, out['concentration'])

        # Test with pairs and high cutoff
        pairs_ct_expected = [[9.933000e-01],
                             [9.868000e-01],
                             [9.999984e-01],
                             [9.843713e-01],
                             [9.932675e-01],
                             [9.885687e-01],
                             [9.980516e-01],
                             [9.922003e-01],
                             [9.965996e-01],
                             [9.963997e-01],
                             [9.999922e-01],
                             [9.999922e-01],
                             [9.965996e-01],
                             [9.978997e-01],
                             [9.984999e-01]]
        pairs_ct_output = self.nupack.concentrations(dnas, 1e-6, pairs=True,
                                                     cutoff=0.9)
        for expect, out in zip(pairs_ct_expected, pairs_ct_output):
            assert_equal(expect, out['fpairs'])

        # TODO: sort options

    def test_distributions(self):
        dnas = self.dnas[:2]
        # Test with no optional arguments
        dist_cx_expected = [[0, 1],
                            [1, 0],
                            [0, 2],
                            [1, 1],
                            [1, 2],
                            [2, 0],
                            [3, 0],
                            [2, 1],
                            [0, 3],
                            [4, 0],
                            [3, 1],
                            [2, 2],
                            [1, 3]]
        dist_evs_expected = [2.000000e+00,
                             1.000000e+00,
                             1.223070e-14,
                             6.288671e-15,
                             1.170654e-29,
                             0.000000e+00,
                             0.000000e+00,
                             0.000000e+00,
                             0.000000e+00,
                             0.000000e+00,
                             0.000000e+00,
                             0.000000e+00,
                             0.000000e+00]

        dist_probcols_expected = [[1.210143e-14, 6.288671e-15, 1.000000e+00],
                                  [6.217249e-15, 1.000000e+00, 0.000000e+00],
                                  [1.000000e+00, 1.223070e-14, 0.000000e+00],
                                  [1.000000e+00, 6.288671e-15, 0.000000e+00],
                                  [1.000000e+00, 1.170654e-29, 0.000000e+00],
                                  [1.000000e+00, 0.000000e+00, 0.000000e+00],
                                  [1.000000e+00, 0.000000e+00, 0.000000e+00],
                                  [1.000000e+00, 0.000000e+00, 0.000000e+00],
                                  [1.000000e+00, 0.000000e+00, 0.000000e+00],
                                  [1.000000e+00, 0.000000e+00, 0.000000e+00],
                                  [1.000000e+00, 0.000000e+00, 0.000000e+00],
                                  [1.000000e+00, 0.000000e+00, 0.000000e+00],
                                  [1.000000e+00, 0.000000e+00, 0.000000e+00]]

        # complexes = self.nupack.complexes(dnas, 4)
        # output = self.nupack.distributions(complexes, [1, 2], 1e-6)
        # for i, out in enumerate(output):
        #     assert_equal(out, {'complex': dist_cx_expected,
        #                        'ev': dist_evs_expected,
        #                        'probcols': dist_probcols_expected})

    def _process_mfe(self, filename):
        curdir = os.path.dirname(__file__)
        mfepath = os.path.join(curdir, 'data', filename)
        with open(mfepath) as f:
            mfefile = f.read()
        commentline = '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %'
        sections = mfefile.split(commentline)
        output = []
        for section in sections[1::2]:
            lines = section.split('\n')
            # Remove first three lines - information not used
            lines.pop(0)
            lines.pop(0)
            lines.pop(0)
            # Remove the last (empty) line
            lines.pop()
            mfe_data = {}
            mfe_data['mfe'] = float(lines[0])
            mfe_data['dotparens'] = lines[1]
            pairlist = []
            if len(lines) > 2:
                for line in lines[2:]:
                    data = line.split('\t')
                    pairlist.append([int(data[0]) - 1, int(data[1]) - 1])
            mfe_data['pairlist'] = pairlist
            output.append(mfe_data)

        return output

    def _read_cxepairs(self, filename):
        curdir = os.path.dirname(__file__)
        epairspath = os.path.join(curdir, 'data', filename)
        with open(epairspath) as f:
            epairsfile = f.read()
        commentline = '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %'
        sections = epairsfile.split(commentline)
        output = []
        for section in sections[1::2]:
            # Skipping every other one eliminates comments and blank lines
            lines = section.split('\n')
            # Remove first three lines - information not used
            lines.pop(0)
            lines.pop(0)
            lines.pop(0)
            # Remove the last (empty) line
            lines.pop()
            # Extract data
            complex_data = []
            for line in lines:
                data = line.split('\t')
                base_i = int(data[0])
                base_j = int(data[1])
                prob = float(data[2])
                complex_data.append((base_i, base_j, prob))
            output.append(complex_data)

        return output

    def _process_ppairs(self, filename, dim):
        curdir = os.path.dirname(__file__)
        tsvpath = os.path.join(curdir, 'data', filename)
        with open(tsvpath) as f:
            pairlist = [row for row in csv.reader(f, delimiter='\t')]

        return self._pairs_to_np(pairlist, dim)

    def _pairs_to_np(self, pairlist, dim):
        '''Given a set of pair probability lines, construct a numpy array.

        :param pairlist: a list of pair probability triples
        :type pairlist: list

        '''
        mat = np.zeros((dim, dim + 1))
        for line in pairlist:
            i = int(line[0]) - 1
            j = int(line[1]) - 1
            prob = float(line[2])
            mat[i, j] = prob
        return mat
