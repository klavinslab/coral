'''Test gibson cloning module.'''
from nose.tools import assert_equal, assert_raises
import coral as cr


def test_gibson_primers():
    '''Test gibson_primers function.'''
    # Fuse tdh3 promoter sequence to yfp (trimmed for readability)
    tdh3_3prime = cr.DNA('aaccagttccctgaaattattcccctacttgactaataagtat'
                         'ataaagacggtaggtattgattgtaattctgtaaatctatttc'
                         'ttaaacttc')
    yfp_nterm = cr.DNA('atggtgagcaagggcgaggagctgttcaccggggtggtgcccatc'
                       'ctggtcgagctggacggcgacgtaaacggccacaagttcagcgtg'
                       'tccggcgagggcgagggcgatgccacctacggcaagctgaccctg'
                       'aag')
    # Expected annealing sequences and their Tms
    fwd_anneal = cr.DNA('atggtgagcaagggcg')
    fwd_tm = 64.64172107821065
    rev_anneal = cr.DNA('gaagtttaagaaatagatttacagaattacaatcaatac')
    rev_tm = 64.24536287254085
    # Expected overlaps
    all_right = cr.DNA('TCGCCCTTGCTCACCAT')
    all_left = cr.DNA('GGTATTGATTGTAATTCTGTAAATCTATTTCTTAAACTTC')
    mixed_fwd = cr.DNA('TTCTTAAACTTC')
    mixed_rev = cr.DNA('CCTTGCTCACCAT')
    # Design primers - with homology all on left side, right side, or mixed
    # All on the 'right' - i.e. fwd primer
    right = cr.cloning.gibson_primers(tdh3_3prime, yfp_nterm, 'right')
    right_rev = cr.Primer(rev_anneal, tm=rev_tm, overhang=all_right)
    right_fwd = cr.Primer(fwd_anneal, tm=fwd_tm)
    assert_equal(right, (right_rev, right_fwd))
    # All on the 'left' - i.e. rev primer
    left = cr.cloning.gibson_primers(tdh3_3prime, yfp_nterm, 'left')
    left_rev = cr.Primer(rev_anneal, tm=rev_tm)
    left_fwd = cr.Primer(fwd_anneal, tm=fwd_tm, overhang=all_left)
    assert_equal(left, (left_rev, left_fwd))
    # On both primers
    mixed = cr.cloning.gibson_primers(tdh3_3prime, yfp_nterm, 'mixed')
    mixed_primer1 = cr.Primer(rev_anneal, tm=rev_tm, overhang=mixed_rev)
    mixed_primer2 = cr.Primer(fwd_anneal, tm=fwd_tm, overhang=mixed_fwd)
    assert_equal(mixed, (mixed_primer1, mixed_primer2))

    assert_raises(ValueError, cr.cloning.gibson_primers, tdh3_3prime,
                  yfp_nterm, 'duck')
