'''Utils for analysis module.'''
import coral


def sequence_type(seq):
    '''Validates a coral.sequence data type.

    :param sequence_in: input DNA sequence.
    :type sequence_in: any
    :returns: The material - 'dna', 'rna', or 'peptide'.
    :rtype: str
    :raises: ValueError

    '''
    if isinstance(seq, coral.DNA):
        material = 'dna'
    elif isinstance(seq, coral.RNA):
        material = 'rna'
    elif isinstance(seq, coral.Peptide):
        material = 'peptide'
    else:
        raise ValueError('Input was not a recognized coral.sequence object.')
    return material
