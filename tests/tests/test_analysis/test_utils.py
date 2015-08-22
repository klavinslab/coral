'''
Tests for utils submodule of the analysis module.

'''

from nose.tools import assert_equal, assert_raises
from pymbt import analysis, DNA, RNA, Peptide


def test_utils():
    test_DNA = DNA('ATAGCGATACGAT')
    test_RNA = RNA('AUGCGAUAGCGAU')
    test_peptide = Peptide('msvkkkpvqg')
    test_str = 'msvkkkpvgq'

    assert_equal(analysis.utils.sequence_type(test_DNA), 'dna')
    assert_equal(analysis.utils.sequence_type(test_RNA), 'rna')
    assert_equal(analysis.utils.sequence_type(test_peptide), 'peptide')
    assert_raises(Exception, analysis.utils.sequence_type, test_str)
