'''RNA sequences classes.'''
import coral as cr
from . import alphabets
from ._nucleicacid import NucleicAcid


class RNA(NucleicAcid):
    '''ssRNA sequence.'''
    def __init__(self, rna, alphabet=alphabets.rna, circular=False,
                 run_checks=True):
        '''
        :param rna: Input sequence (RNA).
        :type rna: str
        :param alphabet: Alphabet to use for this RNA sequence (defaults to
                         \'AUGCN-\').
        :type alphabet: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :returns: coral.RNA instance.

        '''
        super(RNA, self).__init__(rna, 'rna', alphabet=alphabet,
                                  circular=circular, run_checks=run_checks,
                                  any_char='N')
        self.material = 'rna'

    def copy(self):
        return type(self)(self.seq, alphabet=self.alphabet,
                          circular=self.circular, run_checks=False)

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
