'''Provides the Structure class for convenient structural analysis
questions.'''


from coral.analysis import ViennaRNA


class Structure(object):

    def __init__(self, mode='viennarna'):
        allowed = ['viennarna']
        if mode not in allowed:
            raise ValueError('Accepted values for mode are {}'.format(allowed))
        self.mode = mode
        self.calculator = ViennaRNA()

    def mfe(self, strand):
        if self.mode == 'viennarna':
            return self.calculator.fold(strand)['mfe']
