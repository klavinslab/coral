'''RNA sequences classes.'''
import collections

import coral as cr
from . import alphabets
from ._nucleicacid import NucleicAcid


class RNA(NucleicAcid):
    '''ssRNA sequence.'''

    def __init__(self, sequence, alphabet=alphabets.rna, circular=False,
                 skip_checks=False):
        '''
        :param sequence: Input sequence (RNA).
        :type sequence: str
        :param alphabet: Alphabet to use for this RNA sequence (defaults to
                         \'AUGCN-\').
        :type alphabet: cr.Alphabet
        :param skip_checks: Skips input checking (alphabet check), useful for
                            computationally intense tasks.
        :type skip_checks: bool
        :returns: coral.RNA instance.

        '''
        super(RNA, self).__init__(sequence, alphabet=alphabet,
                                  circular=circular, skip_checks=skip_checks,
                                  any_char='N')

    def copy(self):
        return type(self)(self.seq, alphabet=self.alphabet,
                          circular=self.circular, skip_checks=True)

    def mw(self):
        '''Calculate the molecular weight.

        :returns: The molecular weight of the current sequence in amu.
        :rtype: float

        '''
        counter = collections.Counter(self.seq.lower())
        # TODO: use any_char, not n
        # TODO: add to constants module
        mw_a = counter['a'] * 326.2
        mw_g = counter['g'] * 305.2
        mw_c = counter['c'] * 345.2
        mw_u = counter['u'] * 306.2

        mw_n = counter['n'] * (313.2 + 289.2 + 329.2 + 306.2) / 5
        return mw_a + mw_u + mw_g + mw_c + mw_n + 159.0

    def reverse_transcribe(self):
        '''Reverse transcribe to DNA.

        :returns: The reverse transcribed (DNA) version of the current RNA.
        :rtype: coral.DNA

        '''
        return cr.reaction.reverse_transcribe(self)

    def translate(self):
        '''Translate sequence into a peptide.

        :returns: A translated peptide from the current sequence.
        :rtype: coral.Peptide

        '''
        return cr.reaction.translate(self)
