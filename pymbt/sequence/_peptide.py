'''Peptide module.'''
from ._sequence import BaseSequence


class Peptide(BaseSequence):
    '''Peptide sequence.'''
    def __init__(self, peptide, run_checks=True):
        '''
        :param peptide: Input sequence (peptide).
        :type peptide: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :returns: pymbt.Peptide instance.

        '''
        super(Peptide, self).__init__(peptide, 'peptide',
                                      run_checks=run_checks)

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely editable copy of the current sequence.
        :rtype: pymbt.Peptide

        '''
        # Significant performance improvements by skipping alphabet check
        return type(self)(self._sequence, run_checks=False)

    def extract(self, name, pure=False):
        return super(Peptide, self).extract(self, name, 'X', pure=True)

    def __contains__(self, query):
        '''Defines `query in sequence` operator.

        :param query: query string or DNA sequence
        :type query: str or pymbt.DNA

        '''
        return super(Peptide, self).__contains__(query, 'X')
