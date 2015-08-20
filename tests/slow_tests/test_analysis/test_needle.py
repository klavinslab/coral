'''
Tests for Needleman-Wunsch function 'needle'

'''

from nose.tools import assert_equal
from pymbt import analysis, DNA


def test_needle():
    ref_seq = DNA("ATGCGATACGATA")

    res_seq0 = DNA("ATGCGATA---TA")  # Gapped
    res_seq1 = DNA("ATGCGATAATGCGATA")  # Insertion
    res_seq2 = DNA("ATGCGATATA")  # Deletion
    res_seq3 = DNA("ATGCGATAAGATA")  # Mismatch
    results = [res_seq0, res_seq1, res_seq2, res_seq3]

    exp_seq0 = (DNA("ATGCGATACGATA"), DNA("ATGCGATA---TA"), 9)
    exp_seq1 = (DNA("ATGCGATA---CGATA"), DNA("ATGCGATAATGCGATA"), 12)
    exp_seq2 = (DNA("ATGCGATACGATA"), DNA("ATGCGATA---TA"), 9)
    exp_seq3 = (DNA("ATGCGATACGATA"), DNA("ATGCGATAAGATA"), 11)

    expected = [exp_seq0, exp_seq1, exp_seq2, exp_seq3]

    for seq, exp in zip(results, expected):
        aligned = analysis.needle(ref_seq, seq, gap_open=-1, gap_extend=0)
        assert_equal(aligned, exp)
