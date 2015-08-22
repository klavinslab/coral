'''Test restriction reaction module.'''

from nose.tools import assert_equal
from pymbt import reaction, DNA, RestrictionSite


class TestDigest(object):
    '''Test digest function.'''
    def __init__(self):
        # Contains NcoI site
        self.dna = DNA('TGACCATGGAAA')

    def test_not_found(self):
        '''If site not found, should return input sequence in list.'''
        ecorv = RestrictionSite(DNA('GATATC'), (3, 3), name='EcoRV')
        assert_equal(self.dna, reaction.digest(self.dna, ecorv)[0])

    def test_ncoi_cut(self):
        '''Test standard TypeII cutter.'''
        ncoi = RestrictionSite(DNA('CCATGG'), (1, 5), name='NcoI')
        assert_equal(reaction.digest(self.dna, ncoi),
                     [DNA('TGAC----', bottom='CATGGTCA'),
                      DNA('CATGGAAA', bottom='TTTC----')])
        assert_equal(reaction.digest(self.dna.circularize(), ncoi),
                     [DNA('CATGGAAATGAC----', bottom='CATGGTCATTTC----')])

    def test_ecorv_cut(self):
        '''Test blunt-end cutter.'''
        ecorv = RestrictionSite(DNA('GATATC'), (3, 3), name='EcoRV')
        assert_equal(reaction.digest(DNA('GATATC'), ecorv),
                     [DNA('GAT'), DNA('ATC')])

    def test_psti_cut(self):
        '''Test 3\' cutter.'''
        psti = RestrictionSite(DNA('CTGCAG'), (5, 1), name='PstI')
        assert_equal(reaction.digest(DNA('ACTGCAGA'), psti),
                     [DNA('ACTGCA', bottom='----GT'),
                      DNA('----GA', bottom='TCTGCA')])
