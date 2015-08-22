'''
Tests StructureWindows analysis class.

'''

from nose.tools import assert_equal
from pymbt import analysis, DNA


def test_structure_windows():
    '''Tests StructureWindows class in structure_windows.'''

    seq = 'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggc' + \
          'gacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaag' + \
          'ctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgacc' + \
          'accttcggctacggcctgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttc' + \
          'ttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggc' + \
          'aactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctg' + \
          'aagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaac' + \
          'agccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatc' + \
          'cgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatc' + \
          'ggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaa' + \
          'gaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcact' + \
          'ctcggcatggacgagctgtacaagtaa'
    dna_seq = DNA(seq)
    walker = analysis.StructureWindows(dna_seq)
    walker.windows(window_size=60, context_len=90, step=10)
    assert_equal(walker.scores,
                 (0.578570075,
                  0.5928413833333335,
                  0.5535072916666667,
                  0.5425574666666667,
                  0.6028716333333335,
                  0.5907444666666667,
                  0.5532209166666666,
                  0.5882098916666667,
                  0.6471799,
                  0.6957834999999999,
                  0.6209094583333334,
                  0.5929873583333332,
                  0.6117790833333332,
                  0.6116499166666667,
                  0.5987705999999998,
                  0.6439044999999999,
                  0.6817365833333334,
                  0.6488576499999998,
                  0.6900404249999998,
                  0.6657639999999999,
                  0.7083993333333333,
                  0.6360369916666666,
                  0.6452116666666665,
                  0.6395126666666666,
                  0.6288818333333333,
                  0.6351839999999999,
                  0.6463396666666666,
                  0.6717609166666665,
                  0.67853025,
                  0.7012450833333332,
                  0.6620117499999998,
                  0.7250783333333332,
                  0.6995034166666668,
                  0.7386933333333333,
                  0.7494905833333333,
                  0.7247731666666668,
                  0.7510857500000001,
                  0.7458025000000003,
                  0.7434455,
                  0.6702263583333334,
                  0.6390452499999999,
                  0.6503500249999998,
                  0.646285175,
                  0.606586825,
                  0.5707148,
                  0.644573625,
                  0.6644399750000001,
                  0.6716777749999999,
                  0.6807071583333334))
