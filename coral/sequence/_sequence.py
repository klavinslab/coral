'''Base sequence classes.'''
import re
from coral.constants.genbank import TO_CORAL
from coral.constants.molecular_bio import ALPHABETS, COMPLEMENTS


class Sequence(object):
    '''Abstract representation of single chain of molecular sequences, e.g.
       a single DNA or RNA strand or Peptide.'''
    def __init__(self, sequence, material, run_checks=True, any_char='N'):
        '''
        :param sequence: Input sequence.
        :type sequence: str
        :param material: Material type (dna, rna, peptide)
        :type material: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :param any_char: Character representing "any", e.g. N for DNA.
        :type any_char: str
        :type run_checks: bool
        :returns: coral.sequence.Sequence instance.

        '''
        if run_checks:
            self.seq = process_seq(sequence, material)
        else:
            self.seq = sequence

        self.material = material
        self.any_char = any_char

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely editable copy of the current sequence.

        '''
        # Significant performance improvements by skipping alphabet check
        return type(self)(self.seq, self.material, run_checks=False)

    def locate(self, pattern):
        '''Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str
        :returns: Indices of pattern matches.
        :rtype: list of ints

        '''
        if len(pattern) > len(self):
            raise ValueError('Pattern too long.')
        pattern = str(pattern).upper()
        re_pattern = '(?=' + pattern + ')'
        return [index.start() for index in
                re.finditer(re_pattern, self.seq)]

    def __add__(self, other):
        '''Defines addition.

        :param other: Instance with which to sum.
        :type other: coral.sequence.Sequence
        :returns: Concatenated sequence.
        :rtype: coral.sequence.Sequence

        '''
        if type(self) != type(other):
            try:
                other = type(self)(other)
            except AttributeError:
                raise TypeError('Cannot add {} to {}'.format(self, other))

        return type(self)(self.seq + other.seq, self.material,
                          run_checks=False)

    def __contains__(self, query):
        '''`x in y`.

        :param query: Query (i.e. exact pattern) sequence to find.
        :type query: str
        :returns: Whether the query is found in the current sequence.
        :rtype: bool

        '''
        query_str = str(query).upper()
        query_str = re.sub(self.any_char, '.', query_str)
        if re.search(query_str, str(self)):
            return True
        else:
            return False

    def __delitem__(self, key):
        '''Deletes sequence at location key.

        :param key: Index to delete
        :type key: int
        :returns: The current sequence with the moiety at `key` removed.
        :rtype: coral.sequence.Sequence

        '''
        sequence_list = list(self.seq)
        del sequence_list[key]
        self.seq = ''.join(sequence_list)

    def __eq__(self, other):
        '''Define == operator. True if sequences are the same.

        :param other: Other sequence.
        :type other: coral.sequence._sequence.Sequence
        :returns: Whether two sequences have the same base string (sequence).
        :rtype: bool

        '''
        if str(self) == str(other):
            return True
        else:
            return False

    def __getitem__(self, key):
        '''Indexing and slicing of sequences.

        :param key: int or slice for subsetting.
        :type key: int or slice
        :returns: Slice of the current sequence.
        :rtype: coral.sequence.Sequence

        '''
        # If empty sequence, return empty sequence
        if not self.seq or self.seq is None:
            if not isinstance(key, slice):
                raise IndexError()
            return self.copy()
        new_seq = self.seq[key]
        return type(self)(new_seq, self.material, run_checks=False)

    def __len__(self):
        '''Calculate sequence length.

        :returns: The length of the sequence.
        :rtype: int

        '''
        return len(self.seq)

    def __mul__(self, n):
        '''Concatenate copies of the sequence.

        :param n: Factor by which to multiply the sequence.
        :type n: int
        :returns: The current sequence repeated n times.
        :rtype: coral.sequence.Sequence
        :raises: TypeError if n is not an integer.

        '''
        # Input checking
        if n != int(n):
            raise TypeError('Multiplication by non-integer.')
        return sum([x for x in _decompose(self, n)])

    def __ne__(self, other):
        '''Define != operator.

        :param other: Other sequence.
        :type other: coral.sequence._sequence.Sequence
        :returns: The opposite of ==.
        :rtype: bool

        '''
        try:
            return not (self == other)
        except TypeError:
            return False

    def __radd__(self, other):
        '''Add unlike types (enables sum function).

        :param other: Object of any other type.
        :param other: coral.sequence._sequence.Sequence
        :returns: Concatenated sequence.
        :rtype: coral.sequence.Sequence

        '''
        if other == 0 or other is None:
            # For compatibility with sum()
            return self
        elif type(self) != type(other):
            try:
                other = type(self)(other)
            except AttributeError:
                raise TypeError('Cannot add {} to {}'.format(self, other))
        return self + other

    def __repr__(self):
        '''String to print when object is called directly.'''
        display_bases = 40
        if len(self.seq) < 90:
            sequence = self.seq
        else:
            sequence = ''.join([self.seq[:display_bases], ' ... ',
                                self.seq[-display_bases:]])
        return str(sequence)

    def __setitem__(self, index, new_value):
        '''Sets index value to new value.

        :param index: Index to modify.
        :type index: int
        :param new_value: Value to input.
        :type new_value: str or coral.sequence._sequence.Sequence
        :returns: The current sequence with the moiety at `index` replaced
                  by `new_value`.
        :rtype: coral.sequence.Sequence

        '''
        sequence_list = list(self.seq)
        sequence_list[index] = str(type(self)(new_value, self.material))
        self.seq = ''.join(sequence_list)

    def __str__(self):
        '''Cast to string.

        :returns: A string of the current sequence
        :rtype: str

        '''
        return self.seq


class NucleicAcidSequence(Sequence):
    '''Abstract sequence container for a single nucleic acid sequence
       molecule.'''
    def reverse_complement(self):
        copy = self.copy()
        copy.seq = reverse_complement(self.seq, self.material)
        return copy

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


def _decompose(string, n):
    '''Given string and multiplier n, find m**2 decomposition.

    :param string: input string
    :type string: str
    :param n: multiplier
    :type n: int
    :returns: generator that produces m**2 * string if m**2 is a factor of n
    :rtype: generator of 0 or 1

    '''
    binary = [int(x) for x in bin(n)[2:]]
    new_string = string
    counter = 1
    while counter <= len(binary):
        if binary[-counter]:
            yield new_string
        new_string += new_string
        counter += 1


class Feature(object):
    '''Represent A DNA feature - annotate and extract sequence by metadata.'''
    def __init__(self, name, start, stop, feature_type, gene='', locus_tag='',
                 qualifiers={}, strand=0, gaps=[]):
        '''
        :param name: Name of the feature. Used during feature extraction.
        :type name: str
        :param start: Where the feature starts
        :type start: int
        :param stop: Where the feature stops
        :type stop: int
        :param feature_type: The type of the feature. Allowed types:
                                'coding', 'primer', 'promoter', 'terminator',
                                'rbs'
        :type name: str
        :param strand: Watson (0) or Crick (1) strand of the feature.
        :type strand: int
        :param gaps: Gap locations if the feature has gaps.
        :type gaps: list of coordinates (2-tuple/list)
        :param gene: gene attribute (Genbank standard) - gene name (e.g. galK
                     on MG1655 genome).
        :type gene: str
        :param locus_tag: locus_tag attribute (Genbank standard) - systematic
                          locus name (e.g. b0757 for galK on MG1655 genome).
        :type locus_tag: str
        :param qualifiers: Complete Genbank qualifiers key:value pairs
        :type qualifiers: dict
        :returns: coral.Feature instance.
        :raises: ValueError if `feature_type` is not in
                 coral.constants.genbank.TO_CORAL.

        '''
        self.name = name
        self.start = int(start)
        self.stop = int(stop)
        self.modified = False
        self.gene = gene
        self.locus_tag = locus_tag
        self.qualifiers = qualifiers
        self.strand = strand
        self.gaps = gaps

        allowed_types = list(sorted(TO_CORAL.keys()))

        if feature_type in allowed_types:
            self.feature_type = feature_type
        else:
            msg1 = 'feature_type '
            msg2 = 'must be one of the following: {}'.format(allowed_types)
            raise ValueError(msg1 + msg2)

    def move(self, bases):
        '''Move the start and stop positions.

        :param bases: bases to move - can be negative
        :type bases: int

        '''
        self.start += bases
        self.stop += bases

    def copy(self):
        '''Return a copy of the Feature.

        :returns: A safely editable copy of the current feature.
        :rtype: coral.Feature

        '''
        return type(self)(self.name, self.start, self.stop, self.feature_type,
                          gene=self.gene, locus_tag=self.locus_tag,
                          qualifiers=self.qualifiers, strand=self.strand)

    def __repr__(self):
        '''Represent a feature.'''
        if self.modified:
            part1 = '(Modified) {} "{}" feature '.format(self.name,
                                                         self.feature_type)
        else:
            part1 = '{} "{}" feature '.format(self.name, self.feature_type)
        part2 = '({0} to {1}) on strand {2}'.format(self.start, self.stop,
                                                    self.strand)
        return part1 + part2

    def __eq__(self, other):
        '''Define equality.

        :returns: Whether the name and feature type are the same.
        :rtype: bool

        '''
        if self.name != other.name:
            return False
        if self.feature_type != other.feature_type:
            return False
        if self.start != other.start:
            return False
        if self.stop != other.stop:
            return False
        if self.strand != other.strand:
            return False
        if self.gaps != other.gaps:
            return False

        return True

    def __ne__(self, other):
        '''Define inequality.'''
        if not self == other:
            return True
        else:
            return False


def reverse_complement(sequence, material):
    '''Reverse complement a sequence.

    :param sequence: Sequence to reverse complement
    :type sequence: str
    :param material: dna, rna, or peptide.
    :type material: str
    '''
    code = dict(COMPLEMENTS[material])
    reverse_sequence = sequence[::-1]
    return ''.join([code[str(base)] for base in reverse_sequence])


def check_alphabet(seq, material):
    '''Verify that a given string is valid DNA, RNA, or peptide characters.

    :param seq: DNA, RNA, or peptide sequence.
    :type seq: str
    :param material: Input material - 'dna', 'rna', or 'pepide'.
    :type sequence: str
    :returns: Whether the `seq` is a valid string of `material`.
    :rtype: bool
    :raises: ValueError if `material` isn't "dna", "rna", or "peptide".
             ValueError if `seq` contains invalid characters for its
             material type.

    '''
    errs = {'dna': 'DNA', 'rna': 'RNA', 'peptide': 'peptide'}
    if material == 'dna' or material == 'rna' or material == 'peptide':
        alphabet = ALPHABETS[material]
        err_msg = errs[material]
    else:
        msg = 'Input material must be "dna", "rna", or "peptide".'
        raise ValueError(msg)
    # This is a bottleneck when modifying sequence - hence the run_checks
    # optional parameter in sequence objects..
    # First attempt with cython was slower. Could also try pypy.
    if re.search('[^' + alphabet + ']', seq):
        raise ValueError('Encountered a non-%s character' % err_msg)


def process_seq(seq, material):
    '''Validate and process sequence inputs.

    :param seq: input sequence
    :type seq: str
    :param material: DNA, RNA, or peptide
    :type: str
    :returns: Uppercase version of `seq` with the alphabet checked by
              check_alphabet().
    :rtype: str

    '''
    check_alphabet(seq, material)
    seq = seq.upper()
    return seq


def palindrome(seq):
    '''Test whether a sequence is palindrome.

    :param seq: Sequence to analyze (DNA or RNA).
    :type seq: coral.DNA or coral.RNA
    :returns: Whether a sequence is a palindrome.
    :rtype: bool

    '''
    seq_len = len(seq)
    if seq_len % 2 == 0:
        # Sequence has even number of bases, can test non-overlapping seqs
        wing = seq_len / 2
        l_wing = seq[0: wing]
        r_wing = seq[wing:]
        if l_wing == r_wing.reverse_complement():
            return True
        else:
            return False
    else:
        # Sequence has odd number of bases and cannot be a palindrome
        return False
