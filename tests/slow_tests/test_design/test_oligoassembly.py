'''
Tests for the OligoAssembly design class.

'''

from nose.tools import assert_equal, assert_raises, assert_true
from pymbt import design, DNA


def test_oligo_assembly():
    '''
    Tests output of OligoAssembly class.

    '''

    seq = 'atgcgtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgat' + \
          'gttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgcaacatacggaaaactt' + \
          'acccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactact' + \
          'ttcggttatggtgttcaatgctttgcgagatacccagatcatatgaaacagcatgactttttc' + \
          'aagagtgccatgcccgaaggttatgtacaggaaagaactatatttttcaaagatgacgggaac' + \
          'tacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaa' + \
          'ggtattgattttaaagaagatggaaacattcttggacacaaattggaatacaactataactca' + \
          'cacaatgtatacatcatggcagacaaacaaaagaatggaatcaaagttaacttcaaaattaga' + \
          'cacaacattgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggc' + \
          'gatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttcgaaagat' + \
          'cccaacgaaaagagagaccacatggtccttcttgagtttgtaacagctgctgggattacacat' + \
          'ggcatggatgaactatacaaaaggcctgctgcaaacgacgaaaactacgctttagtagcttaa'
    dna_seq = DNA(seq)

    # Test to make sure oligo_number parameter is working
    reference_oligos = ['atgcgtaaaggagaagaacttttcactggagttgtcccaattcttgttg' +
                        'aattagatggtgatgttaatgggcacaaattttctgtcagtgg',
                        'tttccagtagtgcaaataaatttaagggtaagttttccgtatgttgcat' +
                        'caccttcaccctctccactgacagaaaatttgtgcccattaaca',
                        'cggaaaacttacccttaaatttatttgcactactggaaaactacctgtt' +
                        'ccatggccaacacttgtcactactttcggttatggtgttcaatgc',
                        'ttcgggcatggcactcttgaaaaagtcatgctgtttcatatgatctggg' +
                        'tatctcgcaaagcattgaacaccataaccgaaagtagtgacaagt',
                        'ttcaagagtgccatgcccgaaggttatgtacaggaaagaactatatttt' +
                        'tcaaagatgacgggaactacaagacacgtgctgaagtcaagtttgaa',
                        'agaatgtttccatcttctttaaaatcaataccttttaactcgattctat' +
                        'taacaagggtatcaccttcaaacttgacttcagcacgtgtcttgtagt',
                        'tcgagttaaaaggtattgattttaaagaagatggaaacattcttggaca' +
                        'caaattggaatacaactataactcacacaatgtatacatcatggcaga',
                        'tcttcaatgttgtgtctaattttgaagttaactttgattccattctttt' +
                        'gtttgtctgccatgatgtatacattgtgtgagttatagttgtattc',
                        'ggaatcaaagttaacttcaaaattagacacaacattgaagatggaagcg' +
                        'ttcaactagcagaccattatcaacaaaatactccaattggcgatgg',
                        'cgttgggatctttcgaaagggcagattgtgtggacaggtaatggttgtc' +
                        'tggtaaaaggacagggccatcgccaattggagtattttgttga',
                        'tgccctttcgaaagatcccaacgaaaagagagaccacatggtccttctt' +
                        'gagtttgtaacagctgctgggattacacatggcat',
                        'ttaagctactaaagcgtagttttcgtcgtttgcagcaggccttttgtat' +
                        'agttcatccatgccatgtgtaatcccagcagctg']

    oligo_n_assembly = design.OligoAssembly(dna_seq,
                                            tm=72,
                                            length_range=(100, 160),
                                            require_even=True,
                                            start_5=True,
                                            oligo_number=12)
    oligo_n_assembly.design_assembly()

    oligo_n_output = [str(oligo).lower() for oligo in oligo_n_assembly.oligos]
    assert_equal(oligo_n_output, reference_oligos)

    # Test to make sure oligo_number parameter fails with too restrictive of
    # settings
    def design_impossible(test_seq):
        assembly = design.OligoAssembly(test_seq,
                                        tm=72,
                                        length_range=(120, 120),
                                        require_even=True,
                                        start_5=True,
                                        oligo_number=12)
        assembly.design_assembly()

    assert_raises(Exception, design_impossible, dna_seq)


def test_overlapping_overlaps():
    '''
    Sometimes, an assembly produces a result with 'overlapping overlaps' -
    not ideal. This should eventually be replaced by a catchable exception or
    prevented outright.

    '''

    test_seq = DNA('ATCAATACTTATTACGATATATATAT' * 34)
    oligo_n_assembly = design.OligoAssembly(test_seq,
                                            tm=65,
                                            length_range=(80, 150),
                                            require_even=True,
                                            start_5=True,
                                            overlap_min=20,
                                            oligo_number=10)
    oligo_n_assembly.design_assembly()
    assert_true(type(oligo_n_assembly.warning) == str)
