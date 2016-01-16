'''CRISPR-Cas design methods. The current goal is to capture the design methods
listed in the AddGene guide: http://www.addgene.org/CRISPR/guide/'''
import coral


def locate_pams(sequence, pam='NGG'):
    '''Locate PAM sequences within a target sequence.

    :param sequence: Sequence to search for PAMs.
    :type sequence: coral.DNA
    :param pam: PAM sequence (
    :type pam: str or coral.DNA
    '''
    return coral.DNA(str(sequence)).locate(pam)
