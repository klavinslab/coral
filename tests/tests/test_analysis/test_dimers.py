'''Test dimers submodule of analysis module.'''
from pymbt import analysis, DNA, Primer
from nose.tools import assert_equal


def test_dimers():
    '''Test dimers function.'''
    anneal_f = DNA('gatcgatcgatacgatcgatatgcgat', stranded='ss')
    tm_f = 71.86183729637946
    primer_f = Primer(anneal_f, tm_f)

    anneal_r = DNA('atatcgatcatatcgcatatcgatcgtatcgat', stranded='ss')
    tm_r = 72.14300162714233
    primer_r = Primer(anneal_r, tm_r)

    dimer_output = analysis.dimers(primer_f, primer_r)

    assert_equal(dimer_output, 0.8529446)
