'''Tests for Needleman-Wunsch function.'''

import coral as cr
from nose.tools import assert_equal


def test_needle():
    ref_seq = cr.DNA('ATGCGATACGATA')

    res_seq0 = cr.DNA('ATGCGATA---TA')  # Gapped
    res_seq1 = cr.DNA('ATGCGATAATGCGATA')  # Insertion
    res_seq2 = cr.DNA('ATGCGATATA')  # Deletion
    res_seq3 = cr.DNA('ATGCGATAAGATA')  # Mismatch
    results = [res_seq0, res_seq1, res_seq2, res_seq3]

    exp_seq0 = (cr.DNA('ATGCGATACGATA'), cr.DNA('ATGCGATA---TA'), 9)
    exp_seq1 = (cr.DNA('ATGCGATA---CGATA'), cr.DNA('ATGCGATAATGCGATA'), 12)
    exp_seq2 = (cr.DNA('ATGCGATACGATA'), cr.DNA('ATGCGATA---TA'), 9)
    exp_seq3 = (cr.DNA('ATGCGATACGATA'), cr.DNA('ATGCGATAAGATA'), 11)

    expected = [exp_seq0, exp_seq1, exp_seq2, exp_seq3]

    for seq, exp in zip(results, expected):
        aligned = cr.alignment.needle(ref_seq, seq, gap_open=-1, gap_extend=0)
        assert_equal(aligned, exp)
