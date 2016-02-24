'''Generate a random DNA sequence.'''
import random
import coral as cr
import coral.sequence.alphabets as alphabets


def random_dna(n):
    '''Generate a random DNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random DNA sequence of length n.
    :rtype: coral.DNA

    '''
    residues = [random.choice(alphabets.dna_unambiguous) for i in range(n)]
    return cr.DNA(''.join(residues))


def random_ssdna(n):
    '''Generate a random ssDNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random ssDNA sequence of length n.
    :rtype: coral.ssDNA

    '''
    residues = [random.choice(alphabets.dna_unambiguous) for i in range(n)]
    return cr.ssDNA(''.join(residues))


def random_rna(n):
    '''Generate a random RNA sequence.

    :param n: Number of bases.
    :type n: int
    :returns: Random RNA sequence of length n.
    :rtype: coral.RNA

    '''
    residues = [random.choice(alphabets.rna_unambiguous) for i in range(n)]
    return cr.RNA(''.join(residues))


def random_peptide(n):
    '''Generate a random DNA sequence.

    :param n: Output sequence length.
    :type n: int
    :returns: Random Peptide sequence of length n.
    :rtype: coral.Peptide

    '''
    residues = [random.choice(alphabets.peptide) for i in range(n)]
    return cr.Peptide(''.join(residues))
