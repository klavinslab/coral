'''
Tests utils submodule of reaction module.

'''

from nose.tools import assert_equal, assert_raises
from pymbt import reaction, DNA


def test_convert_sequence():
    '''Tests DNA translation function.'''

    seq = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGC' + \
          'GACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAG' + \
          'CTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACC' + \
          'ACCTTCGGCTACGGCCTGCAGTGCTTCGCCCGCTACCCCGACCACATGAAGCAGCACGACTTC' + \
          'TTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGC' + \
          'AACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTG' + \
          'AAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAAC' + \
          'AGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATC' + \
          'CGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATC' + \
          'GGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCTACCAGTCCGCCCTGAGCAAA' + \
          'GACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACT' + \
          'CTCGGCATGGACGAGCTGTACAAGTAA'
    dna = DNA(seq)
    prot = 'MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLV' + \
           'TTFGYGLQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRI' + \
           'ELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQN' + \
           'TPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK'
    rna = reaction.utils.convert_sequence(dna, 'rna')
    r_trans = reaction.utils.convert_sequence(rna, 'dna')
    trans = reaction.utils.convert_sequence(rna, 'peptide')
    assert_equal(str(trans), prot)
    assert_equal(str(r_trans), seq)
    assert_raises(ValueError, reaction.utils.convert_sequence, seq, 'rna')

    # Gapped sequence shouldfail
    assert_raises(ValueError, reaction.utils.convert_sequence,
                  DNA('atg-'), 'rna')

    # Sequence without stop codon should still work
    nostop_dna = DNA('atgaaaaaaaaaaaa')
    nostop_rna = reaction.utils.convert_sequence(nostop_dna, 'rna')
    nostop_peptide = reaction.utils.convert_sequence(nostop_rna, 'peptide')
    assert_equal(str(nostop_rna), 'AUGAAAAAAAAAAAA')
    assert_equal(str(nostop_peptide), 'MKKKK')

    assert_raises(ValueError, reaction.utils.convert_sequence, 'duck', 'rna')
