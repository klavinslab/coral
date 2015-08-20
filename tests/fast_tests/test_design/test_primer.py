'''Tests primer design module.'''
import warnings
from nose.tools import assert_equals, assert_not_equal, assert_raises
from pymbt import design, DNA


def test_primer():
    '''Test primer function.'''
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
    dna_seq = DNA(seq)
    primer = design.primer(dna_seq, tm=72, min_len=10, tm_undershoot=1,
                           tm_overshoot=3, end_gc=False,
                           tm_parameters='cloning', overhang=None)
    assert_equals(str(primer), 'ATGGTGAGCAAGGGCGAGGAG')
    # Ensure that overhang is appropriately applied
    overhang_primer = design.primer(dna_seq, tm=72, min_len=10,
                                    tm_undershoot=1, tm_overshoot=3,
                                    end_gc=False,
                                    tm_parameters='cloning',
                                    overhang=DNA('GATCGATAT'))
    assert_equals(str(overhang_primer), 'GATCGATATATGGTGAGCAAGGGCGAGGAG')
    # If sequence is too short (too low of Tm), raise ValueError
    too_short = DNA('at')
    assert_raises(ValueError, design.primer, too_short, tm=72)
    # Should design different primers (sometimes) if ending on GC is preferred
    diff_template = DNA('GATCGATCGATACGATCGATATGCGATATGATCGATAT')
    nogc = design.primer(diff_template, tm=72, min_len=10,
                         tm_undershoot=1, tm_overshoot=3, end_gc=False,
                         tm_parameters='cloning', overhang=None)
    withgc = design.primer(diff_template, tm=72, min_len=10,
                           tm_undershoot=1, tm_overshoot=3,
                           end_gc=True, tm_parameters='cloning', overhang=None)
    assert_not_equal(nogc, withgc)
    # Should raise ValueError if it's impossible to create an end_gc primer
    end_at_template = DNA('ATGCGATACGATACGCGATATGATATATatatatat' +
                          'ATAAaaaaaaaaaattttttttTTTTTTTTTTTTTT' +
                          'TTTTTTTTTT')
    assert_raises(ValueError, design.primer, end_at_template,
                  end_gc=True, tm=72)
    # If there's structure, should issue a warning
    structure_template = DNA('ATGCGATCGATAGGCGA')
    structure_template += structure_template.reverse_complement()
    with warnings.catch_warnings(True) as w:
        design.primer(structure_template, structure=True, tm=72)
        assert len(w) > 0


def test_primers():
    '''Test primers function.'''
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
    dna_seq = DNA(seq)
    primers_list = design.primers(dna_seq, tm=72, min_len=10,
                                  tm_undershoot=1, tm_overshoot=3,
                                  end_gc=False, tm_parameters='cloning',
                                  overhangs=None)
    primers = [str(x.primer()) for x in primers_list]
    assert_equals(primers, ['ATGGTGAGCAAGGGCGAGGAG',
                            'TTACTTGTACAGCTCGTCCATGCCG'])
