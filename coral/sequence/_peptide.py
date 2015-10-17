'''Peptide module.'''
from ._sequence import Sequence


class Peptide(object):
    '''Peptide sequence.'''
    def __init__(self, peptide, features=None, run_checks=True):
        '''
        :param peptide: Input sequence (peptide).
        :type peptide: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :returns: coral.Peptide instance.

        '''
        self._sequence = Sequence(peptide, 'peptide', run_checks=run_checks,
                                  any_char='X')
        if features is None:
            self.features = []
        else:
            self.features = features

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely editable copy of the current sequence.
        :rtype: coral.Peptide

        '''
        return type(self)(str(self._sequence), features=self.features,
                          run_checks=False)

    @classmethod
    def extract(self, feature, remove_subfeatures=False):
        return self._sequence.extract(self, feature,
                                      remove_subfeatures=remove_subfeatures)

    def locate(self, pattern):
        return self._sequence.locate(pattern)

    def __add__(self, other):
        # Merge the sequences
        copy = self.copy()
        ocopy = other.copy()
        copy._sequence += ocopy._sequence

        # Merge the features
        for feature in ocopy.features:
            feature.move(len(copy))
        copy.features += ocopy.features

        return copy

    def __contains__(self, key):
        return key in self._sequence

    def __delitem__(self, key):
        self._sequence.__delitem__(key)

    def __eq__(self, other):
        return self._sequence == other._sequence

    def __getitem__(self, key):
        new_instance = type(self)(str(self._sequence[key]),
                                  features=self.features, run_checks=False)
        return new_instance

    def __len__(self):
        return len(self._sequence)

    def __mul__(self, n):
        # TODO: Keep features as well?
        copy = self.copy()
        copy._sequence = self._sequence.__mul__(n)
        return copy

    def __ne__(self, other):
        return self._sequence != other._sequence

    def __radd__(self, other):
        if other == 0 or other is None:
            # For compatibility with sum()
            return self
        else:
            copy = self.copy()
            try:
                return copy._sequence.__radd__(other._sequence)
            except AttributeError:
                raise TypeError('Cannot add {} to {}'.format(self, other))

    def __repr__(self):
        '''String to print when object is called directly.'''
        header = 'Peptide:'
        sequence = self._sequence.__repr__()
        return ' '.join([header, sequence])

    def __setitem__(self, key, value):
        self._sequence.__setitem__(key, value)

    def __str__(self):
        return str(self._sequence)
