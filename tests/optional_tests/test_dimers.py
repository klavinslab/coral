'''Test dimers submodule of analysis module.'''
import coral as cr
from nose.tools import assert_equal


def test_dimers():
    '''Test dimers function.'''
    anneal_f = cr.ssDNA('gatcgatcgatacgatcgatatgcgat')
    tm_f = 71.86183729637946
    primer_f = cr.Primer(anneal_f, tm_f)

    anneal_r = cr.ssDNA('atatcgatcatatcgcatatcgatcgtatcgat')
    tm_r = 72.14300162714233
    primer_r = cr.Primer(anneal_r, tm_r)

    dimer_output = cr.structure.dimers(primer_f, primer_r)

    assert_equal(dimer_output, 0.8529446)
