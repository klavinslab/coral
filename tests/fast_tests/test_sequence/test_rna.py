'''
Tests for the RNA sequence class.

'''

from pymbt import RNA
from nose.tools import assert_equal, assert_false, assert_true, assert_raises


class TestRNA(object):
    '''
    Testing class for sequence.RNA

    '''

    def __init__(self):
        self.test_rna = RNA('augc')

    def test_reverse_complement(self):
        assert_equal(str(self.test_rna.reverse_complement()), 'GCAU')

    def test_locate(self):
        assert_equal(self.test_rna.locate('au'), [0])
        assert_equal(self.test_rna.locate('gc'), [2])
        assert_equal(len(self.test_rna.locate('augg')), 0)

    def test_copy(self):
        assert_equal(self.test_rna, self.test_rna.copy())

    def test_getitem(self):
        assert_equal(str(self.test_rna[0]), 'A')
        assert_equal(str(self.test_rna[1]), 'U')
        assert_equal(str(self.test_rna[2]), 'G')
        assert_equal(str(self.test_rna[3]), 'C')
        assert_equal(str(self.test_rna[-1]), 'C')

    def test_delitem(self):
        copy0 = self.test_rna.copy()
        del copy0[0]
        assert_equal(str(copy0), 'UGC')
        copy1 = self.test_rna.copy()
        del copy1[1]
        assert_equal(str(copy1), 'AGC')
        copy2 = self.test_rna.copy()
        del copy2[2]
        assert_equal(str(copy2), 'AUC')
        copy3 = self.test_rna.copy()
        del copy3[3]
        assert_equal(str(copy3), 'AUG')
        copy_1 = self.test_rna.copy()
        del copy_1[-1]
        assert_equal(str(copy_1), 'AUG')

    def test_setitem(self):
        copy0 = self.test_rna.copy()
        copy0[0] = 'u'
        assert_equal(str(copy0), 'UUGC')
        copy1 = self.test_rna.copy()
        copy1[1] = 'a'
        assert_equal(str(copy1), 'AAGC')
        copy2 = self.test_rna.copy()
        copy2[2] = 'a'
        assert_equal(str(copy2), 'AUAC')
        copy3 = self.test_rna.copy()
        copy3[3] = 'a'
        assert_equal(str(copy3), 'AUGA')
        copy_1 = self.test_rna.copy()
        copy_1[-1] = 'a'
        assert_equal(str(copy_1), 'AUGA')

    def test_str(self):
        assert_equal(str(self.test_rna), 'AUGC')

    def test_len(self):
        assert_equal(len(self.test_rna), 4)

    def test_add(self):
        assert_equal(str((self.test_rna + self.test_rna)), 'AUGCAUGC')

    def test_radd(self):
        assert_equal(str(sum([self.test_rna, self.test_rna])), 'AUGCAUGC')

        def radd_800(seq):
            return 800 + seq

        assert_raises(TypeError, radd_800, self.test_rna)

    def test_mul(self):
        assert_equal(str((self.test_rna * 4)), 'AUGCAUGCAUGCAUGC')

        def mul_float(seq):
            return seq * 7.56

        assert_raises(TypeError, mul_float, self.test_rna)

    def test_eq(self):
        assert_true(self.test_rna == RNA('augc'))

    def test_ne(self):
        assert_true(self.test_rna != RNA('aagc'))

    def test_contains(self):
        assert_true('a' in self.test_rna)
        assert_true('u' in self.test_rna)
        assert_true('g' in self.test_rna)
        assert_true('c' in self.test_rna)
        assert_false('t' in self.test_rna)
