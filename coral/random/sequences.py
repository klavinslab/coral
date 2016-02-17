'''Generate a random DNA sequence.'''
import random
import coral as cr


def random_dna(n):
    '''Generate a random DNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random DNA sequence of length n.
    :rtype: coral.DNA

    '''
    return cr.DNA(''.join([random.choice('ATGC') for i in range(n)]))


def random_ssdna(n):
    '''Generate a random ssDNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random ssDNA sequence of length n.
    :rtype: coral.ssDNA

    '''
    return cr.ssDNA(''.join([random.choice('ATGC') for i in range(n)]))


def random_rna(n):
    '''Generate a random RNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random RNA sequence of length n.
    :rtype: coral.RNA

    '''
    return cr.RNA(''.join([random.choice('ATGC') for i in range(n)]))


def random_peptide(n):
    '''Generate a random DNA sequence.

    :param n: Output sequence length.
    :type n: int
    :returns: Random Peptide sequence of length n.
    :rtype: coral.Peptide

    '''
    return cr.Peptide(''.join([random.choice('ATGC') for i in range(n)]))
