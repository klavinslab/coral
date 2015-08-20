'''Tests for the DNA sequence class.'''
from pymbt import reaction, DNA, Feature, RestrictionSite
from nose.tools import assert_equal, assert_false, assert_true, assert_raises
from nose.tools import assert_not_equal


class TestDNA(object):
    '''Testing class for DNA'''
    def __init__(self):
        self.test_dna = DNA('atgc')

    def test_reverse_complement(self):
        assert_equal(str(self.test_dna.reverse_complement()), 'GCAT')

    def test_linearize(self):
        circ_dna = self.test_dna.circularize()
        assert_equal(str(circ_dna.linearize()), 'ATGC')
        assert_equal(str(circ_dna.linearize(0)), 'ATGC')
        assert_equal(str(circ_dna.linearize(1)), 'TGCA')
        assert_equal(str(circ_dna.linearize(2)), 'GCAT')
        assert_equal(str(circ_dna.linearize(3)), 'CATG')
        assert_equal(str(circ_dna.linearize(-1)), 'CATG')
        assert_raises(ValueError, circ_dna.linearize().linearize)

    def test_to_ss_ds(self):
        assert_equal(self.test_dna.to_ds(), self.test_dna)
        ss_dna = self.test_dna.to_ss()
        assert_equal(ss_dna.stranded, 'ss')
        ds_dna = self.test_dna.to_ds()
        assert_equal(ds_dna.stranded, 'ds')
        assert_equal(ds_dna.top(), str(ss_dna))

        ds_to_ss_to_ds = self.test_dna.to_ss().to_ds()
        assert_equal(self.test_dna, ds_to_ss_to_ds)

        empty_top = reaction.three_resect(self.test_dna, 400)
        assert_equal(empty_top.to_ds(), self.test_dna)

    def test_locate(self):
        assert_equal(self.test_dna.locate('a'), [[0], [2]])
        assert_equal(self.test_dna.locate('at'), [[0], [2]])
        assert_equal(self.test_dna.locate('gc'), [[2], [0]])
        assert_equal(self.test_dna.locate('atgg'), [[], []])
        # Circular DNA tests
        assert_equal(self.test_dna.circularize().locate('a'), [[0], [2]])
        assert_equal(self.test_dna.circularize().locate('at'), [[0], [2]])
        assert_equal(self.test_dna.circularize().locate('gc'), [[2], [0]])
        assert_equal(self.test_dna.circularize().locate('atgg'), [[], []])

    def test_copy(self):
        assert_equal(self.test_dna, self.test_dna.copy())

    def test_palindrome(self):
        palindromic_seq_even = DNA('ATGCGCAT')
        nonpalindromic_seq_even = DNA('ATGCGCAA')
        almost_palindrome_odd = DNA('ATGCCAT')
        assert_true(palindromic_seq_even.is_palindrome())
        assert_false(nonpalindromic_seq_even.is_palindrome())
        assert_false(almost_palindrome_odd.is_palindrome())

    def test_getitem(self):
        assert_equal(str(self.test_dna[0]), 'A')
        assert_equal(str(self.test_dna[1]), 'T')
        assert_equal(str(self.test_dna[2]), 'G')
        assert_equal(str(self.test_dna[3]), 'C')
        assert_equal(str(self.test_dna[-1]), 'C')

    def test_delitem(self):
        copy0 = self.test_dna.copy()
        del copy0[0]
        assert_equal(str(copy0), 'TGC')
        copy1 = self.test_dna.copy()
        del copy1[1]
        assert_equal(str(copy1), 'AGC')
        copy2 = self.test_dna.copy()
        del copy2[2]
        assert_equal(str(copy2), 'ATC')
        copy3 = self.test_dna.copy()
        del copy3[3]
        assert_equal(str(copy3), 'ATG')
        copy_1 = self.test_dna.copy()
        del copy_1[-1]
        assert_equal(str(copy_1), 'ATG')

    def test_setitem(self):
        copy0 = self.test_dna.copy()
        copy0[0] = 't'
        assert_equal(str(copy0), 'TTGC')
        copy1 = self.test_dna.copy()
        copy1[1] = 'a'
        assert_equal(str(copy1), 'AAGC')
        copy2 = self.test_dna.copy()
        copy2[2] = 'a'
        assert_equal(str(copy2), 'ATAC')
        copy3 = self.test_dna.copy()
        copy3[3] = 'a'
        assert_equal(str(copy3), 'ATGA')
        copy_1 = self.test_dna.copy()
        copy_1[-1] = 'a'
        assert_equal(str(copy_1), 'ATGA')

        def set_gap(seq):
            seq[2] = '-'
        assert_raises(ValueError, set_gap, self.test_dna)

    def test_repr(self):
        expected_repr = 'linear dsDNA:\nATGC\nTACG'
        assert_equal(repr(self.test_dna), expected_repr)

        expected_circ_repr = 'circular dsDNA:\nATGC\nTACG'
        assert_equal(repr(self.test_dna.circularize()), expected_circ_repr)

        repr_1 = 'linear dsDNA:\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC ... '
        repr_2 = 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n'
        repr_3 = 'TACGTACGTACGTACGTACGTACGTACGTACGTACGTACG ... '
        repr_4 = 'TACGTACGTACGTACGTACGTACGTACGTACGTACGTACG'
        expected_long_repr = repr_1 + repr_2 + repr_3 + repr_4
        assert_equal(repr(self.test_dna * 50), expected_long_repr)

    def test_str(self):
        assert_equal(str(self.test_dna), 'ATGC')

    def test_len(self):
        assert_equal(len(self.test_dna), 4)

    def test_add(self):
        assert_equal(str(self.test_dna + self.test_dna), 'ATGCATGC')
        assert_equal(str(self.test_dna.to_ss() +
                         self.test_dna.to_ss()), 'ATGCATGC')

    def test_radd(self):
        assert_equal(str(sum([self.test_dna, self.test_dna])), 'ATGCATGC')

        def radd_800(seq):
            return 800 + seq

        assert_raises(TypeError, radd_800, self.test_dna)

    def test_mul(self):
        assert_equal(str(self.test_dna * 4), 'ATGCATGCATGCATGC')

        def mul_float(seq):
            return seq * 7.56

        assert_raises(TypeError, mul_float, self.test_dna)

        # TODO: reimplement this test using manual sequence input
        # def mul_incompatible(seq):
        #     return seq * 3

        # incompatible_seq = self.test_dna.copy()
        # incompatible_seq = incompatible_seq.five_resect(1)
        # incompatible_seq = incompatible_seq.reverse_complement()
        # incompatible_seq = incompatible_seq.five_resect(1)
        # incompatible_seq = incompatible_seq.reverse_complement()
        # assert_raises(Exception, mul_incompatible, incompatible_seq)

    def test_eq(self):
        assert_true(self.test_dna == DNA('atgc'))

    def test_ne(self):
        assert_true(self.test_dna != DNA('aagc'))

    def test_contains(self):
        assert_true('a' in self.test_dna)
        assert_true('t' in self.test_dna)
        assert_true('g' in self.test_dna)
        assert_true('c' in self.test_dna)
        assert_false('u' in self.test_dna)

    def test_flip(self):
        flipped = self.test_dna.flip()
        assert_equal(str(self.test_dna), flipped._bottom)
        assert_equal(self.test_dna._bottom, str(flipped))


def test_stranded_init():
    ss_dna = DNA('atgc', stranded='ss')
    assert_true(all([base == '-' for base in ss_dna.bottom()]))

    ds_dna = DNA('atgc')
    assert_equal(str(ds_dna), ds_dna.reverse_complement().bottom())


def test_stranded_complemented():
    ss_dna = DNA('atgc', stranded='ss')
    r_ss_dna = ss_dna.reverse_complement()
    assert_equal(r_ss_dna.top(), 'GCAT')
    assert_equal(r_ss_dna.bottom(), '----')


class TestFeatures(object):
    '''Test features model using DNA object.'''
    def __init__(self):
        self.dna = DNA('ATGC') * 50
        self.apply_features()

    def apply_features(self):
        misc_feature = Feature('Misc Feature', 1, 20, 'misc_feature')
        misc_1_feature = Feature('Misc Feature', 1, 20, 'misc_feature',
                                 strand=1)
        coding_feature = Feature('Coding Feature', 21, 40, 'CDS')
        primer_feature = Feature('Primer Feature', 41, 60, 'primer_bind')
        promoter_feature = Feature('Promoter Feature', 61, 80, 'promoter')
        terminator_feature = Feature('Terminator Feature', 81, 100,
                                     'terminator')
        rbs_feature = Feature('RBS Feature', 101, 120, 'RBS')
        origin_feature = Feature('Origin Feature', 121, 140, 'rep_origin')
        utr3_feature = Feature("3'UTR Feature", 141, 160, "3'UTR")
        origin_feature2 = Feature('Origin Feature', 161, 180, 'rep_origin')

        input_features = [misc_feature, misc_1_feature, coding_feature,
                          primer_feature, promoter_feature, terminator_feature,
                          rbs_feature, origin_feature, utr3_feature,
                          origin_feature2]
        self.dna.features = input_features

    def test_good_features(self):
        for feature in self.dna.features:
            assert_true(feature.copy() in self.dna.features)

    def test_rev_comp(self):
        rev = self.dna.reverse_complement()
        for feature, rev_feature in zip(self.dna.features, rev.features):
            assert_not_equal(feature.strand, rev_feature.strand)
            assert_equal(len(self.dna) - feature.start, rev_feature.stop)
            assert_equal(len(self.dna) - feature.stop, rev_feature.start)

    def test_extract(self):
        test_utr3_feature = [feature for feature in self.dna.features if
                            feature.name == "3'UTR Feature"][0]
        extracted = self.dna.extract(test_utr3_feature)
        assert_equal(str(extracted), 'TGCATGCATGCATGCATGC')

    def test_getitem(self):
        subsequence = self.dna[30:100]
        remaining_features = [Feature('Primer Feature', 11, 30, 'primer_bind'),
                              Feature('Promoter Feature', 31, 50, 'promoter'),
                              Feature('Terminator Feature', 51, 70,
                                      'terminator')]

        assert_equal(subsequence.features, remaining_features)
        assert_false(self.dna[10].features)
        new_seq = DNA('ATGC', features=[Feature('A', 0, 0, 'misc_feature')])
        assert_equal(new_seq[0].features[0],
                     Feature('A', 0, 0, 'misc_feature'))

    def test_ne(self):
        '''Test != operator'''
        assert_true(self.dna.features[0] != self.dna.features[4])
        assert_false(self.dna.features[0] != self.dna.features[0])


class TestRestrictionSite(object):
    '''Test RestrictionSite class.'''
    def __init__(self):
        self.ecorv = RestrictionSite(DNA('GATATC'), (3, 3), name='EcoRV')
        self.foki = RestrictionSite(DNA('GGATG'), (14, 18), name='FokI')

    def test_cuts_outside(self):
        '''Test cuts_outside method.'''
        assert_false(self.ecorv.cuts_outside())
        assert_true(self.foki.cuts_outside())

    def test_len(self):
        '''Test len function.'''
        assert_equal(len(self.ecorv), 6)
        assert_equal(len(self.foki), 5)
