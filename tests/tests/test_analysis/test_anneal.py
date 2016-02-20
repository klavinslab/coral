'''Test functionality of PCR class of reaction module.'''
import os
import coral as cr
from coral import analysis, seqio, DNA, Primer
from nose.tools import assert_true, assert_raises


def test_basic():
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                                           "pMODKan-HO-pACT1GEV.ape"))
    # Test forward priming.
    seq = DNA('cgccagggttttcccagtcacgac')
    primer = Primer(seq, 50.6)
    matches = analysis.anneal(template, primer)
    fwd_matches, rev_matches = matches
    fwd_indices = [match[0] for match in fwd_matches]
    # fwd_lens = [match[1] for match in fwd_matches]
    rev_indices = [match[0] for match in rev_matches]
    # rev_lens = [match[1] for match in rev_matches]

    loc = template.locate(seq)

    assert_true(len(fwd_matches) == len(loc[0]))
    assert_true(len(rev_matches) == len(loc[1]))

    # Top strand matches
    for match in loc[0]:
        assert_true(match + len(seq) in fwd_indices)
    # Top strand matches
    for match in loc[1]:
        assert_true(match + len(seq) in rev_indices)

    # Test reverse priming
    seq = DNA('ACAAGAGAGATTGGGAAGGAAAGGATCA')
    primer = Primer(seq, 50.6)
    matches = analysis.anneal(template, primer)
    fwd_matches, rev_matches = matches
    fwd_indices = [match[0] for match in fwd_matches]
    rev_indices = [match[0] for match in rev_matches]

    loc = template.locate(seq)

    assert_true(len(fwd_indices) == len(loc[0]))
    assert_true(len(rev_indices) == len(loc[1]))
    for match in loc[0]:
        assert_true(match + len(seq) in fwd_indices)
    for match in loc[1]:
        assert_true(match + len(seq) in rev_indices)


def test_near_index():
    '''Test binding near index for circular templates.'''
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                                           "pMODKan-HO-pACT1GEV.ape"))
    template = template.circularize()
    seq = DNA('aggccctttcgtctcgcgcgttt')
    primer = Primer(seq, 50.6)
    matches = analysis.anneal(template, primer)
    fwd_matches, rev_matches = matches
    fwd_indices = [match[0] for match in fwd_matches]
    rev_indices = [match[0] for match in rev_matches]

    loc = template.locate(seq)

    print fwd_matches
    print loc[0]
    print loc[1]
    assert_true(len(fwd_matches) == len(loc[0]))
    assert_true(len(rev_matches) == len(loc[1]))
    for match in loc[0]:
        expected = match + len(seq)
        if expected > len(template):
            expected -= len(template)
        assert_true(expected in fwd_indices)
    for match in loc[1]:
        expected = match + len(seq)
        if expected > len(template):
            expected -= len(template)
        assert_true(expected in rev_indices)


def test_overhang():
    ''' test forward priming '''
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                                           "pMODKan-HO-pACT1GEV.ape"))
    seq = DNA('cgccagggttttcccagtcacgac')
    overhang = DNA('ggggggg')
    seq2 = overhang + seq
    primer = Primer(seq2, 50.6)
    matches = analysis.anneal(template, primer)
    fwd_matches, rev_matches = matches

    fwd_indices = [match[0] for match in fwd_matches]
    rev_indices = [match[0] for match in rev_matches]

    loc = template.locate(seq)

    assert_true(len(fwd_indices) == len(loc[0]))
    assert_true(len(rev_indices) == len(loc[1]))
    # FIXME: Add match length check for all these cases.
    for match in loc[0]:
        assert_true(match + len(seq) in fwd_indices)
    for match in loc[1]:
        assert_true(match + len(seq) in rev_indices)

    # Test forward priming.
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                                           "pMODKan-HO-pACT1GEV.ape"))
    seq = DNA('ACAAGAGAGATTGGGAAGGAAAGGATCA')
    overhang = DNA('ggggggg')
    seq2 = overhang + seq
    primer = Primer(seq2, 50.6)
    matches = analysis.anneal(template, primer)
    fwd_matches, rev_matches = matches
    fwd_indices = [match[0] for match in fwd_matches]
    rev_indices = [match[0] for match in rev_matches]

    loc = template.locate(seq)

    assert_true(len(fwd_indices) == len(loc[0]))
    assert_true(len(rev_indices) == len(loc[1]))

    for match in loc[0]:
        assert_true(match + len(seq) in fwd_indices)
    for match in loc[1]:
        assert_true(match + len(seq) in rev_indices)


def test_multiple_priming():
    ''' test multiple binding sites '''

    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                                           "pMODKan-HO-pACT1GEV.ape"))
    template = template.circularize()
    seq = DNA('cgccagggttttcccagtcacgac')
    template = template.linearize()
    template = template + seq + DNA("AGGCGTATGC") + seq
    template = (template + DNA("GGGGGGG") + seq.reverse_complement() +
                DNA("GGAAAG"))
    template = template.circularize()
    primer = Primer(seq, 50.6)
    matches = analysis.anneal(template, primer)
    fwd_matches, rev_matches = matches
    fwd_indices = [match[0] for match in fwd_matches]
    rev_indices = [match[0] for match in rev_matches]

    loc = template.locate(seq)

    assert_true(len(fwd_matches) == len(loc[0]))
    assert_true(len(rev_matches) == len(loc[1]))
    for match in loc[0]:
        assert_true(match + len(seq) in fwd_indices)
    for match in loc[1]:
        assert_true(match + len(seq) in rev_indices)


def test_no_priming():
    ''' test no priming '''

    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                                           "pMODKan-HO-pACT1GEV.ape"))
    seq = DNA('ggaggagggcggcgaggcgagcgacggaggggga')
    primer = Primer(seq, 50.6)
    matches = analysis.anneal(template, primer)
    fwd_matches, rev_matches = matches
    loc = template.locate(seq)
    assert_true(len(fwd_matches) == len(loc[0]))
    assert_true(len(rev_matches) == len(loc[1]))
    for match in loc[0]:
        assert_true(match + len(seq) in fwd_matches)
    for match in loc[1]:
        assert_true(match + len(seq) in rev_matches)


def test_min_primer_length():
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                                           "pMODKan-HO-pACT1GEV.ape"))

    # Test forward priming
    seq = DNA('cgccagggttttcccagtcacgac')
    seq = seq[:15]
    primer = Primer(seq, 50.6)
    assert_raises(analysis._sequence.anneal.PrimerLengthError, analysis.anneal,
                  template, primer, min_len=16)


def test_min_tm():
    current_path = os.path.dirname(__file__)
    template = seqio.read_dna(os.path.join(current_path,
                                           "pMODKan-HO-pACT1GEV.ape"))

    # Test forward priming
    # Tm should be ~40 C
    seq = DNA('CTTCTATCGAACAA')
    primer = Primer(seq, seq.tm())
    matches = analysis.anneal(template, primer, min_tm=60.0)
    assert_true(len(matches[0]) == 0)
    matches = analysis.anneal(template, primer, min_tm=30.0)
    assert_true(len(matches[0]) > 0)


def test_primer_larger_than_template():
    template = cr.design.random_dna(50)
    overhangs = [cr.design.random_dna(200), cr.DNA('')]
    expected = overhangs[0] + template
    primer1, primer2 = cr.design.primers(template, overhangs=overhangs,
                                         min_len=14)
    amplicon = cr.reaction.pcr(template, primer1, primer2)

    assert_true(expected == amplicon)
