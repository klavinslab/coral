'''
Tests for RandomCodons class of analysis module.

'''

from nose.tools import assert_equal, assert_not_equal, assert_raises
import coral as cr


def test_randomcodons():
    '''
    This test is pretty basic right now - not sure how much checking
    can be done for a random DNA base generator.

    '''

    reference_seq = cr.RNA('AUGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUAG')
    reference_peptide = cr.reaction.translate(reference_seq)
    output = cr.random.random_codons(reference_peptide)
    output_peptide = cr.reaction.translate(reference_seq)

    assert_equal(len(output), len(reference_seq) - 3)
    assert_equal(reference_peptide, output_peptide)
    assert_not_equal(reference_seq, output)

    # Setting too high a threshold should raise ValueError
    assert_raises(ValueError, cr.random.random_codons, reference_peptide,
                  frequency_cutoff=1.5)

    # Weighted should work
    w_output = cr.random.random_codons(reference_peptide, weighted=True)
    w_output_peptide = cr.reaction.translate(reference_seq)

    assert_equal(len(w_output), len(reference_seq) - 3)
    assert_equal(reference_peptide, w_output_peptide)
    assert_not_equal(reference_seq, w_output)
