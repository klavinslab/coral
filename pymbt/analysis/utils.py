'''Utils for analysis module.'''
import pymbt


def sequence_type(seq):
    '''Validates a pymbt.sequence data type.

    :param sequence_in: input DNA sequence.
    :type sequence_in: any
    :returns: The material - 'dna', 'rna', or 'peptide'.
    :rtype: str
    :raises: ValueError

    '''
    if isinstance(seq, pymbt.DNA):
        material = 'dna'
    elif isinstance(seq, pymbt.RNA):
        material = 'rna'
    elif isinstance(seq, pymbt.Peptide):
        material = 'peptide'
    else:
        raise ValueError('Input was not a recognized pymbt.sequence object.')
    return material
