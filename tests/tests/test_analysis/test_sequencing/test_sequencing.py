'''Test sequencing module.'''
import os
from pymbt import analysis, seqio
from nose.tools import assert_equal
# IDEA: unit tests for plotting?


class TestSanger(object):
    '''Test Sanger class.'''
    def __init__(self):
        ref_name = 'reference_sequence.gb'
        reference_path = os.path.join(os.path.dirname(__file__), ref_name)
        results_path = os.path.dirname(__file__)
        self.reference = self.read_reference(reference_path)
        self.results = self.read_results(results_path)

    def read_results(self, path):
        '''Read in sequencing results.'''
        seqs = seqio.read_sequencing(path)
        return seqs

    def read_reference(self, path):
        '''Read in sequencing reference sequence.'''
        seq = seqio.read_dna(path)
        return seq

    def test_sanger(self):
        '''Test basic Sanger methods and init'''
        self.sanger = analysis.Sanger(self.reference, self.results)
        self.sanger.report()
        # There should be one mismatch at position 381
        assert_equal(self.sanger.mismatches, [[(381, 381)], []])
        # There should be a two-base insertion at position 763
        assert_equal(self.sanger.insertions, [[(763, 765)], []])
        # There should be a two-base deletion at position 440
        assert_equal(self.sanger.deletions, [[(440, 441)], []])

    def test_input_single(self):
        '''If users inputs single sequence as a result, handle seamlessly.'''
        analysis.Sanger(self.reference, self.results[0])

    def test_alignment(self):
        '''Ensure that alignment is consistent.'''
        # TODO: this
        pass
