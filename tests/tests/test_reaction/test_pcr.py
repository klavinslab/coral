'''Test functionality of PCR class of reaction module.'''

import os
from coral import design, reaction, seqio, DNA, Primer
from nose.tools import assert_equal, assert_true, assert_raises


def test_basic():
    to_amplify = 'atgtctaaaggtgaagaattattcactggtgttgtcccaatgctgctggtattacc' + \
                 'catggtattgatgaattgtacaaatag'
    template = DNA(to_amplify)
    forward, reverse = design.primers(template)

    amplicon = reaction.pcr(template, forward, reverse)
    assert_equal(amplicon, template)

def test_over_origin():
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))
    assert_true(template.topology == "circular")
    primer1 = design.primer(template[-200:])
    primer2 = design.primer(template.reverse_complement()[-200:])
    over_origin = reaction.pcr(template, primer1, primer2)
    expected = template[-200:] + template[0:200]
    assert_equal(str(over_origin), str(expected))


def test_primer_bind_error():
    to_amplify = 'atgtctaaaggtgaagaattattcactggtgttgtcccaatgctgctggtattacc' + \
                 'catggtattgatgaattgtacaaatag'
    template = DNA(to_amplify)
    primer1, primer2 = design.primers(template)
    # Mess up the second primer so it doesn't bind
    primer2.anneal[10:] = "AAAAAAAAAA"
    assert_raises(reaction._pcr.PrimerBindError, reaction.pcr,
                  template, primer1, primer2)

def test_primer_overlap():
    ''' Tests case in which primers overlap (i.e. primer
    dimers) '''

    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))
    p1 = design.primer(template[100:])
    p2 = design.primer(template[:113].reverse_complement())
    reaction.pcr(template, p1, p2, min_bases=10)

def test_primers_are_in_same_direction_error():
    ''' Tests case in which both forward and reverse
        primer bind in the same direction '''

    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))
    p1 = design.primer(template[100:])
    p2 = design.primer(template[150:])
    assert_raises(reaction._pcr.PrimerBindError, reaction.pcr,
                  template, p1, p2)
    assert_raises(reaction._pcr.PrimerBindError, reaction.pcr,
                  template, p2, p1)
    p1 = design.primer(template.reverse_complement()[100:])
    p2 = design.primer(template.reverse_complement()[150:])
    assert_raises(reaction._pcr.PrimerBindError, reaction.pcr,
                  template, p1, p2)
    assert_raises(reaction._pcr.PrimerBindError, reaction.pcr,
                  template, p2, p1)

def test_AmbiguousPrimingError():
    ''' Tests case in which there are multiple primer
        binding sites '''
    s = 'ACGTGCTGTGATGTCGTGTGA'
    s2 = 'AGGCTGGCTGGAGGTTCG'
    template = DNA(s) + DNA(s) + DNA(s2).reverse_complement()
    p1 = Primer(DNA(s), 65)
    p2 = Primer(DNA(s2).reverse_complement(), 65)
    assert_raises(reaction._pcr.AmbiguousPrimingError, reaction.pcr,
                  template, p1, p2)

def test_overhang():
    ''' Tests if overhangs on primers are correctly
        appended to pcr fragments. '''

    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))

    overhang = DNA("AGCGGGGGGGGGCTGGGGCTGAT")
    p1 = design.primer(template[100:])
    p1 = Primer(overhang + p1.primer(), 65) #add overhang

    rev_overhang = DNA("GGGGGGGGGGGGGGGGGGG")
    p2 = design.primer(template[:300].reverse_complement())
    p2 = Primer(rev_overhang + p2.primer(), 65)

    expected = overhang + template[100:300] + rev_overhang.reverse_complement()
    amplicon = reaction.pcr(template, p1, p2)
    assert_equal(str(expected), str(amplicon))

    amplicon = reaction.pcr(template, p2, p1)
    assert_equal(str(expected), str(amplicon))
