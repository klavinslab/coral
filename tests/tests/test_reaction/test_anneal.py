'''Test functionality of PCR class of reaction module.'''

import os
from coral import design, reaction, seqio, DNA, Primer
from nose.tools import assert_equal, assert_true, assert_raises

def test_basic():
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))

    ''' test forward priming '''
    seq = DNA('cgccagggttttcccagtcacgac')
    primer = Primer(seq, 50.6)
    matches = reaction.anneal(template, primer)
    fwd_matches, rev_matches = matches
    loc = template.locate(seq)
    assert_true(len(fwd_matches.keys()) == len(loc[0]))
    assert_true(len(rev_matches.keys()) == len(loc[1]))
    for m in loc[0]:
        assert_true(m+len(seq) in fwd_matches)
    for m in loc[1]:
        assert_true(m+len(seq) in rev_matches)

    ''' test reverse priming '''
    seq = DNA('ACAAGAGAGATTGGGAAGGAAAGGATCA')
    primer = Primer(seq, 50.6)
    matches = reaction.anneal(template, primer)
    fwd_matches, rev_matches = matches
    loc = template.locate(seq)
    assert_true(len(fwd_matches.keys()) == len(loc[0]))
    assert_true(len(rev_matches.keys()) == len(loc[1]))
    for m in loc[0]:
        assert_true(m+len(seq) in fwd_matches)
    for m in loc[1]:
        assert_true(m+len(seq) in rev_matches)

def test_near_index():
    ''' test binding near index for circular templates '''

    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))
    seq = DNA('aggccctttcgtctcgcgcgttt')
    primer = Primer(seq, 50.6)
    matches = reaction.anneal(template, primer)
    fwd_matches, rev_matches = matches
    loc = template.locate(seq)
    assert_true(len(fwd_matches.keys()) == len(loc[0]))
    assert_true(len(rev_matches.keys()) == len(loc[1]))
    for m in loc[0]:
        assert_true(m+len(seq) in fwd_matches)
    for m in loc[1]:
        assert_true(m+len(seq) in rev_matches)

def test_overhang():
    ''' test forward priming '''
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))
    seq = DNA('cgccagggttttcccagtcacgac')
    overhang = DNA('ggggggg')
    seq2 = overhang + seq
    primer = Primer(seq2, 50.6)
    matches = reaction.anneal(template, primer)
    fwd_matches, rev_matches = matches
    loc = template.locate(seq)
    assert_true(len(fwd_matches.keys()) == len(loc[0]))
    assert_true(len(rev_matches.keys()) == len(loc[1]))
    for m in loc[0]:
        assert_true(m+len(seq) in fwd_matches)
        assert_true(fwd_matches[m+len(seq)][0] == primer)
        assert_true(fwd_matches[m+len(seq)][1] == seq.to_ss())
        assert_true(fwd_matches[m+len(seq)][2] == overhang.to_ss())
    for m in loc[1]:
        assert_true(m+len(seq) in rev_matches)
        assert_true(rev_matches[m+len(seq)][0] == primer)
        assert_true(rev_matches[m+len(seq)][1] == seq.to_ss())
        assert_true(rev_matches[m+len(seq)][2] == overhang.to_ss())

    ''' test forward priming '''
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))
    seq = DNA('ACAAGAGAGATTGGGAAGGAAAGGATCA')
    overhang = DNA('ggggggg')
    seq2 = overhang + seq
    primer = Primer(seq2, 50.6)
    matches = reaction.anneal(template, primer)
    fwd_matches, rev_matches = matches
    loc = template.locate(seq)
    assert_true(len(fwd_matches.keys()) == len(loc[0]))
    assert_true(len(rev_matches.keys()) == len(loc[1]))
    for m in loc[0]:
        assert_true(m+len(seq) in fwd_matches)
        assert_true(fwd_matches[m+len(seq)][0] == primer)
        assert_true(fwd_matches[m+len(seq)][1] == seq.to_ss())
        assert_true(fwd_matches[m+len(seq)][2] == overhang.to_ss())
    for m in loc[1]:
        assert_true(m+len(seq) in rev_matches)
        assert_true(rev_matches[m+len(seq)][0] == primer)
        assert_true(rev_matches[m+len(seq)][1] == seq.to_ss())
        assert_true(rev_matches[m+len(seq)][2] == overhang.to_ss())


def test_multiple_priming():
    ''' test multiple binding sites '''

    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))
    seq = DNA('cgccagggttttcccagtcacgac')
    template = template.linearize()
    template = template + seq + DNA("AGGCGTATGC") + seq[5:]
    template = template + DNA("GGGGGGG") + seq.reverse_complement() + DNA("GGAAAG")
    template = template.circularize()
    primer = Primer(seq, 50.6)
    matches = reaction.anneal(template, primer)
    fwd_matches, rev_matches = matches
    loc = template.locate(seq)
    assert_true(len(fwd_matches.keys()) == len(loc[0]))
    assert_true(len(rev_matches.keys()) == len(loc[1]))
    for m in loc[0]:
        assert_true(m+len(seq) in fwd_matches)
    for m in loc[1]:
        assert_true(m+len(seq) in rev_matches)
    print matches

def test_no_priming():
    ''' test no priming '''

    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))
    seq = DNA('ggaggagggcggcgaggcgagcgacggaggggga')
    primer = Primer(seq, 50.6)
    matches = reaction.anneal(template, primer)
    fwd_matches, rev_matches = matches
    loc = template.locate(seq)
    assert_true(len(fwd_matches.keys()) == len(loc[0]))
    assert_true(len(rev_matches.keys()) == len(loc[1]))
    for m in loc[0]:
        assert_true(m+len(seq) in fwd_matches)
    for m in loc[1]:
        assert_true(m+len(seq) in rev_matches)


def test_min_primer_length():
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))

    ''' test forward priming '''
    seq = DNA('cgccagggttttcccagtcacgac')
    seq = seq[:15]
    primer = Primer(seq, 50.6)
    assert_raises(reaction._anneal.PrimerLengthError, reaction.anneal, template, primer, min_bases = 16)

def test_min_tm():
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))

    ''' test forward priming '''
    seq = DNA('CTTCTATCGAACAA')
    primer = Primer(seq, 40)
    matches = reaction.anneal(template, primer, min_tm=60.0)
    assert_true(len(matches[0]) == 0)
    matches = reaction.anneal(template, primer, min_tm=30.0)
    assert_true(len(matches[0]) > 0)

def test_primertypeerror():
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                              "pMODKan-HO-pACT1GEV.ape"))

    dna_seq = DNA('cgccagggttttcccagtcacgac')
    primer = Primer(dna_seq, 65.1)

    assert_raises(reaction._anneal.GenericAnnealError, reaction.anneal,
                  template, dna_seq)

    assert_raises(reaction._anneal.GenericAnnealError, reaction.anneal,
                  Primer(template, 50.6), primer)
