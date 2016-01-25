'''Test functionality of PCR class of reaction module.'''
from nose.tools import assert_equal, assert_raises
import coral as cr


# TODO: refactor into class that shares sequence reading, etc.
class TestPCR(object):
    '''Testing class for reaction.pcr.'''
    def __init__(self):
        # Part BBa_R0010 (pLac promoter)
        bba_r0010 = ('caatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcag'
                     'ctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgca'
                     'attaatgtgagttagctcactcattaggcaccccaggctttacacttta'
                     'tgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttca'
                     'caca')
        self.template = cr.DNA(bba_r0010, circular=False)

    def pcr_equal(self, expected, template, primer1, primer2):
        '''Boilerplate to assert that a pcr reaction matches the expected
        sequence.'''
        amplicon = cr.reaction.pcr(template, primer1, primer2)
        assert_equal(expected, amplicon)

    '''
    Test all the different scenarios:
    Linear template:
      1. Primers point towards each other (whole amplicon and subsets)
      2. Extract region contained between them + their sizes
    Circular template (primers can point towards or away from each other):
      3. Origin (index 0) is between the primers
      4. Origin (index 0) isn't involved at all (amplicon outside it)
      5. Origin isn't between the primers, but forward primer overlaps it.
      6. Origin isn't between the primers, but reverse primer overlaps it.
      7. Origin isn't between the primers, but both primers overlap it.
    '''
    def test_linear_towards_one_another(self):
        '''Amplify entire template, no overhangs.'''
        fwd, rev = cr.design.primers(self.template)
        self.pcr_equal(self.template, self.template, fwd, rev)

    def test_subset(self):
        '''Amplify part of the template, no overhangs.'''
        fwd, rev = cr.design.primers(self.template[30:-30])
        self.pcr_equal(self.template[30:-30], self.template, fwd, rev)

    def test_linear_away_from_one_another(self):
        '''If primers point away from each other on a linear template, raise
        eception.'''
        fwd = cr.design.primer(self.template[-60:])
        rev = cr.design.primer(self.template[:60].reverse_complement())
        assert_raises(cr.reaction._pcr.PrimingError, cr.reaction.pcr,
                      self.template, fwd, rev)

    def test_circular_over_origin(self):
        '''Amplify a circular template over the origin, no overhangs.'''
        template = self.template.circularize()
        fwd = cr.design.primer(template[-60:])
        rev = cr.design.primer(template[:60].reverse_complement())
        expected = template[-60:] + template[:60]
        self.pcr_equal(expected, template, fwd, rev)

    def test_circular_outside_origin(self):
        '''Amplify entire template, no overhangs.'''
        template = self.template.circularize()
        fwd, rev = cr.design.primers(template[30:-30])
        self.pcr_equal(self.template[30:-30], template, fwd, rev)

    def test_circular_fwd_overlap(self):
        '''Test for when the forward primer overlaps the origin (index 0).'''
        template = self.template.circularize()
        fwd = cr.design.primer(template.rotate(4))
        rev = cr.design.primer(template[:60].reverse_complement())
        expected = template[-4:] + template[:60]
        self.pcr_equal(expected, template, fwd, rev)

    def test_circular_rev_overlap(self):
        '''Test for when the reverse primer overlaps the origin (index 0).'''
        template = self.template.circularize()
        fwd = cr.design.primer(template[60:])
        rev = cr.design.primer(template.rotate(-4).reverse_complement())
        expected = template[60:] + template[:4]
        self.pcr_equal(expected, template, fwd, rev)

    def test_circular_fwd_rev_overlap(self):
        '''Test for when both primers overlap the origin (index 0).'''
        template = self.template.circularize()
        fwd = cr.design.primer(template.rotate(4))
        rev = cr.design.primer(template.rotate(-4).reverse_complement())
        expected = template[-4:] + template.linearize() + template[:4]
        self.pcr_equal(expected, template, fwd, rev)

    def test_primer_order(self):
        '''Amplification should occur regardless of the order in which primers
        are specified.'''
        fwd, rev = cr.design.primers(self.template)
        fwd, rev = rev, fwd
        self.pcr_equal(self.template, self.template, fwd, rev)

    def test_primers_overlap(self):
        '''Tests case where primers overlap.'''
        # ~10 bp of overlap
        fwd_seq = self.template[100:120]
        rev_seq = self.template[110:130].reverse_complement()
        fwd = cr.Primer(fwd_seq, fwd_seq.tm())
        rev = cr.Primer(rev_seq, rev_seq.tm())
        # Should amplify 30 bp region
        self.pcr_equal(self.template[100:130], self.template, fwd, rev)

    def test_overhang(self):
        '''Tests that primer overhangs are added correctly to the amplicon.'''
        template = self.template[30:-30]
        fwd_overhang = cr.DNA('AGCGGGGGGGGGCTGGGGCTGAT')
        rev_overhang = cr.DNA('GGGTGGGGGGGGGGGGGGG')
        fwd = cr.design.primer(template, overhang=fwd_overhang)
        rev = cr.design.primer(template.reverse_complement(),
                               overhang=rev_overhang)
        expected = (fwd_overhang +
                    template +
                    rev_overhang.reverse_complement())
        self.pcr_equal(expected, self.template, fwd, rev)

    def test_bad_primer2(self):
        '''Test that a bad primer raises a PrimingError.'''
        fwd, rev = cr.design.primers(self.template)
        rev.anneal[-10:] = 'AAAAAAAAAA'
        assert_raises(cr.reaction._pcr.PrimingError, cr.reaction.pcr,
                      self.template, fwd, rev)

    def test_primers_in_same_direction(self):
        '''Test that primers that bind the same strand raise an error.'''
        fwd1 = cr.design.primer(self.template)
        fwd2 = cr.design.primer(self.template[50:])
        assert_raises(cr.reaction._pcr.PrimingError, cr.reaction.pcr,
                      self.template, fwd1, fwd2)

    def test_ambiguous_priming(self):
        '''Test that ambiguous primer binding sites raises an error.'''
        template = self.template[:50] + self.template
        fwd, rev = cr.design.primers(template)
        assert_raises(cr.reaction._pcr.PrimingError, cr.reaction.pcr,
                      template, fwd, rev)
