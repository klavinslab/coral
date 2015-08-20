'''
Tests for the Tm analysis class.

'''

from nose.tools import assert_equal
from pymbt import analysis, DNA


def test_finnzymes():
    '''
    Tests finnzymes method output.

    '''

    melt = analysis.tm(DNA('ATGCGATAGCGATAGC'), parameters='cloning')
    assert_equal(melt, 55.2370030020752)
