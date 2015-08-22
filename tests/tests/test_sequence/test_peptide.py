'''
Tests for the Peptide sequence class.

'''

from pymbt import Peptide
from nose.tools import assert_equal, assert_false, assert_true, assert_raises


class TestPeptide(object):
    '''
    Testing class for sequence.Peptide

    '''

    def __init__(self):
        self.test_peptide = Peptide('mkgp')

    def test_locate(self):
        assert_equal(self.test_peptide.locate('mk'), [0])
        assert_equal(self.test_peptide.locate('gp'), [2])
        assert_equal(len(self.test_peptide.locate('augg')), 0)

    def test_copy(self):
        assert_equal(self.test_peptide, self.test_peptide.copy())

    def test_getitem(self):
        assert_equal(str(self.test_peptide[0]), 'M')
        assert_equal(str(self.test_peptide[1]), 'K')
        assert_equal(str(self.test_peptide[2]), 'G')
        assert_equal(str(self.test_peptide[3]), 'P')
        assert_equal(str(self.test_peptide[-1]), 'P')

    def test_delitem(self):
        copy0 = self.test_peptide.copy()
        del copy0[0]
        assert_equal(str(copy0), 'KGP')
        copy1 = self.test_peptide.copy()
        del copy1[1]
        assert_equal(str(copy1), 'MGP')
        copy2 = self.test_peptide.copy()
        del copy2[2]
        assert_equal(str(copy2), 'MKP')
        copy3 = self.test_peptide.copy()
        del copy3[3]
        assert_equal(str(copy3), 'MKG')
        copy_1 = self.test_peptide.copy()
        del copy_1[-1]
        assert_equal(str(copy_1), 'MKG')

    def test_setitem(self):
        copy0 = self.test_peptide.copy()
        copy0[0] = 'q'
        assert_equal(str(copy0), 'QKGP')
        copy1 = self.test_peptide.copy()
        copy1[1] = 'q'
        assert_equal(str(copy1), 'MQGP')
        copy2 = self.test_peptide.copy()
        copy2[2] = 'q'
        assert_equal(str(copy2), 'MKQP')
        copy3 = self.test_peptide.copy()
        copy3[3] = 'q'
        assert_equal(str(copy3), 'MKGQ')
        copy_1 = self.test_peptide.copy()
        copy_1[-1] = 'q'
        assert_equal(str(copy_1), 'MKGQ')

    def test_str(self):
        assert_equal(str(self.test_peptide), 'MKGP')

    def test_len(self):
        assert_equal(len(self.test_peptide), 4)

    def test_add(self):
        assert_equal(str((self.test_peptide + self.test_peptide)),
                     'MKGPMKGP')

    def test_radd(self):
        assert_equal(str(sum([self.test_peptide, self.test_peptide])),
                     'MKGPMKGP')

        def radd_800(seq):
            return 800 + seq

        assert_raises(TypeError, radd_800, self.test_peptide)

    def test_mul(self):
        assert_equal(str(self.test_peptide * 4), 'MKGPMKGPMKGPMKGP')

        def mul_float(seq):
            return seq * 7.56

        assert_raises(TypeError, mul_float, self.test_peptide)

    def test_eq(self):
        assert_true(self.test_peptide == Peptide('mkgp'))

    def test_ne(self):
        assert_true(self.test_peptide != Peptide('mkqp'))

    def test_contains(self):
        assert_true('m' in self.test_peptide)
        assert_true('k' in self.test_peptide)
        assert_true('g' in self.test_peptide)
        assert_true('p' in self.test_peptide)
        assert_false('q' in self.test_peptide)
