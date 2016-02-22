'''Tests for the Tm analysis class.'''

from nose.tools import assert_equal
import coral as cr


def test_finnzymes():
    '''Tests finnzymes method output.'''

    melt = cr.thermo.tm(cr.DNA('ATGCGATAGCGATAGC'), parameters='cloning')
    assert_equal(melt, 55.2370030020752)
