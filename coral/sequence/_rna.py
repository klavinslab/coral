'''RNA sequences classes.'''
import coral.reaction
from coral.sequence._nucleicacid import NucleicAcid


class RNA(NucleicAcid):
    '''ssRNA sequence.'''

    def __init__(self, rna, circular=False, run_checks=True):
        '''
        :param rna: Input sequence (RNA).
        :type rna: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :returns: coral.RNA instance.

        '''
        super(RNA, self).__init__(rna, 'rna', circular=circular,
                                  run_checks=run_checks, any_char='N')
        self.ds = False

    def copy(self):
        return type(self)(self.seq, circular=self.circular, run_checks=False)

    def reverse_transcribe(self):
        '''Reverse transcribe to DNA.

        :returns: The reverse transcribed (DNA) version of the current RNA.
        :rtype: coral.DNA

        '''
        return coral.reaction.reverse_transcribe(self)

    def translate(self):
        '''Translate sequence into a peptide.

        :returns: A translated peptide from the current sequence.
        :rtype: coral.Peptide

        '''
        return coral.reaction.translate(self)
