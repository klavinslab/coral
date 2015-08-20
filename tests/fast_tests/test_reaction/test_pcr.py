'''Test functionality of PCR class of reaction module.'''

import os
from pymbt import design, reaction, seqio, DNA
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
    primer2.anneal[0:10] = "AAAAAAAAAA"

    assert_raises(reaction._pcr.PrimerBindError, reaction.pcr,
                  template, primer1, primer2)
