'''Generate a random DNA sequence.'''
import random
import coral as cr
import coral.sequence.alphabets as alphabets


def random_sequence(n, alphabet):
    '''Generate a random sequence given a custom alphabet of symbols (e.g. A,
    T, G, C, X for 5-letter DNA code).

    :param n: Number of letters - the length of the output.
    :type n: int
    :returns: Random sequence of length n. Note that it will return a Sequence
              object, not a convenient container like DNA or RNA. It can be
              used directly by those objects, however - constrained, of course,
              by its alphabet: cr.DNA(random_sequence(n, myalphabet)).
    :rtype: coral.Sequence

    '''
    sequence = ''.join([random.choice(alphabet.symbols) for i in range(n)])
    return cr.Sequence(sequence, alphabet, skip_checks=True)


def random_dna(n):
    '''Generate a random DNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random DNA sequence of length n.
    :rtype: coral.DNA

    '''
    return cr.DNA(random_sequence(n, alphabets.dna_unambiguous))


def random_ssdna(n):
    '''Generate a random ssDNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random ssDNA sequence of length n.
    :rtype: coral.ssDNA

    '''
    return cr.ssDNA(random_sequence(n, alphabets.dna_unambiguous))


def random_rna(n):
    '''Generate a random RNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random RNA sequence of length n.
    :rtype: coral.RNA

    '''
    return cr.RNA(random_sequence(n, alphabets.rna_unambiguous))


def random_peptide(n):
    '''Generate a random DNA sequence.

    :param n: Output sequence length.
    :type n: int
    :returns: Random Peptide sequence of length n.
    :rtype: coral.Peptide

    '''
    return cr.Peptide(random_sequence(n, alphabets.peptide_unambiguous))
