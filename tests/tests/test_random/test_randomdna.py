'''
Tests for RandomBases class of analysis module.

'''

from nose.tools import assert_equal
import coral as cr


def test_randomdna():
    '''
    This test is pretty basic right now - not sure how much checking
    can be done for a random DNA base generator.

    '''

    output = cr.random.random_dna(200)
    assert_equal(len(output), 200)
