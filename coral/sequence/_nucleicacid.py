'''Base sequence classes.'''
import collections
import coral
from coral.sequence._sequence import Sequence
from coral.constants.molecular_bio import COMPLEMENTS


class NucleicAcid(Sequence):
    '''Abstract sequence container for a single nucleic acid sequence
    molecule.'''

    def __init__(self, sequence, material, circular=False, run_checks=True,
                 any_char='N'):
        '''
        :param sequence: Input sequence.
        :type sequence: str
        :param material: Material type (dna, rna)
        :type material: str
        :param circular: The topology of the sequence - if the ends connect,
                         (a circular sequence), set to True. Otherwise, set to
                         False. Enables operations like .rotate().
        :type circular: bool
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :param any_char: Character representing \'any\', e.g. N for DNA.
        :type any_char: str
        :type run_checks: bool
        :returns: coral.sequence.Sequence instance.

        '''
        super(NucleicAcid, self).__init__(sequence, material,
                                          run_checks=run_checks,
                                          any_char=any_char)
        self.circular = circular

    def copy(self):
        return type(self)(self.seq, self.material, circular=self.circular,
                          run_checks=False)

    def circularize(self):
        '''Circularize the sequence, if linear.

        :returns: A circularized version of the current sequence.
        :rtype: coral.sequence._sequence.Sequence

        '''
        copy = self.copy()
        copy.circular = True
        return copy

    def complement(self):
        copy = self.copy()
        code = dict(COMPLEMENTS[self.material])
        copy.seq = ''.join([code[str(base)] for base in copy])
        return copy

    def gc(self):
        '''Find the frequency of G and C in the current sequence.'''
        gc = len([base for base in self.seq if base == 'C' or base == 'G'])
        return float(gc) / len(self)

    def is_palindrome(self):
        seq_len = len(self.seq)
        if seq_len % 2 == 0:
            # Sequence has even number of bases, can test non-overlapping seqs
            wing = seq_len / 2
            l_wing = self[0: wing]
            r_wing = self[wing:]
            if l_wing == r_wing.reverse_complement():
                return True
            else:
                return False
        else:
            # Sequence has odd number of bases and cannot be a palindrome
            return False

    def is_rotation(self, other):
        '''Determine whether two sequences are the same, just at different
        rotations.

        :param other: The sequence to check for rotational equality.
        :type other: coral.sequence._sequence.Sequence

        '''
        if len(self) != len(other):
            return False

        for i in range(len(self)):
            if self.rotate(i) == other:
                return True

        return False

    def linearize(self, index=0):
        '''Linearize the Sequence at an index.

        :param index: index at which to linearize.
        :type index: int
        :returns: A linearized version of the current sequence.
        :rtype: coral.sequence._sequence.Sequence
        :raises: ValueError if the input is a linear sequence.

        '''
        if not self.circular and index != 0:
            raise ValueError('Cannot relinearize a linear sequence.')
        copy = self.copy()
        # Snip at the index
        if index:
            return copy[index:] + copy[:index]
        copy.circular = False

        return copy

    def locate(self, pattern):
        '''Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str
        :returns: Indices of pattern matches.
        :rtype: list of ints

        '''
        if self.circular:
            if len(pattern) >= 2 * len(self):
                raise ValueError('Search pattern longer than searchable ' +
                                 'sequence.')
            seq = self + self[:len(pattern) - 1]
            return super(NucleicAcid, seq).locate(pattern)
        else:
            return super(NucleicAcid, self).locate(pattern)

    def mw(self):
        '''Calculate the molecular weight.

        :returns: The molecular weight of the current sequence in amu.
        :rtype: float

        '''
        counter = collections.Counter(self.seq.lower())
        mw_a = counter['a'] * 313.2
        mw_t = counter['t'] * 304.2
        mw_g = counter['g'] * 289.2
        mw_c = counter['c'] * 329.2
        mw_u = counter['u'] * 306.2

        if self.material == 'dna':
            return mw_a + mw_t + mw_g + mw_c + 79.0
        else:
            return mw_a + mw_u + mw_g + mw_c + 159.0

    def rotate(self, n):
        '''Rotate Sequence by n bases.

        :param n: Number of bases to rotate.
        :type n: int
        :returns: The current sequence reoriented at `index`.
        :rtype: coral.sequence._sequence.Sequence
        :raises: ValueError if applied to linear sequence or `index` is
                 negative.

        '''
        if not self.circular and n != 0:
            raise ValueError('Cannot rotate a linear sequence')
        else:
            rotated = self[-n:] + self[:-n]
            return rotated.circularize()

    def rotate_to(self, index):
        '''Orient Sequence to index (only applies to circular sequences).

        :param index: Position at which to re-zero the Sequence.
        :type index: int
        :returns: The current sequence reoriented at `index`.
        :rtype: coral.sequence._sequence.Sequence
        :raises: ValueError if applied to linear sequence or `index` is
                 negative.

        '''
        return self.rotate(-index)

    def reverse(self):
        return self[::-1]

    def reverse_complement(self):
        copy = self.copy()
        copy.seq = str(self.reverse().complement())
        return copy

    def tm(self, parameters='cloning'):
        '''Find the melting temperature.

        :param parameters: The tm method to use (cloning, santalucia98,
                       breslauer86)
        :type parameters: str

        '''
        return coral.analysis.tm(self, parameters=parameters)
