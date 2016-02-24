'''Peptide module.'''
from ._sequence import Sequence
from . import alphabets


class Peptide(object):
    '''Peptide sequence.'''

    def __init__(self, peptide, alphabet=alphabets.peptide, features=None,
                 run_checks=True):
        '''
        :param peptide: Input sequence (peptide).
        :type peptide: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :returns: coral.Peptide instance.

        '''
        self.sequence = Sequence(peptide, alphabet=alphabet,
                                 run_checks=run_checks, any_char='X')
        self.alphabet = alphabet
        if features is None:
            self.features = []
        else:
            self.features = features

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely editable copy of the current sequence.
        :rtype: coral.Peptide

        '''
        return type(self)(str(self.sequence), alphabet=self.alphabet,
                          features=self.features, run_checks=False)

    @classmethod
    def extract(self, feature, remove_subfeatures=False):
        return self.sequence.extract(self, feature,
                                     remove_subfeatures=remove_subfeatures)

    def locate(self, pattern):
        return self.sequence.locate(pattern)

    def __add__(self, other):
        # Merge the sequences
        copy = self.copy()
        ocopy = other.copy()
        copy.sequence += ocopy.sequence

        # Merge the features
        for feature in ocopy.features:
            feature.move(len(copy))
        copy.features += ocopy.features

        return copy

    def __contains__(self, key):
        return key in self.sequence

    def __delitem__(self, key):
        self.sequence.__delitem__(key)

    def __eq__(self, other):
        return self.sequence == other.sequence

    def __getitem__(self, key):
        new_instance = type(self)(str(self.sequence[key]),
                                  alphabet=self.alphabet,
                                  features=self.features, run_checks=False)
        return new_instance

    def __len__(self):
        return len(self.sequence)

    def __mul__(self, n):
        # TODO: Keep features as well?
        copy = self.copy()
        copy.sequence = self.sequence.__mul__(n)
        return copy

    def __ne__(self, other):
        return self.sequence != other.sequence

    def __radd__(self, other):
        if other == 0 or other is None:
            # For compatibility with sum()
            return self
        else:
            copy = self.copy()
            try:
                return copy.sequence.__radd__(other.sequence)
            except AttributeError:
                raise TypeError('Cannot add {} to {}'.format(self, other))

    def __repr__(self):
        '''String to print when object is called directly.'''
        return self.sequence.__repr__()

    def __setitem__(self, key, value):
        self.sequence.__setitem__(key, value)

    def __str__(self):
        return str(self.sequence)
