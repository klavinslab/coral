'''Write genbank sequences.'''
import textwrap


def write_fasta(sequence, handle):
    '''Write coral.DNA, coral.ssDNA, coral.RNA, or coral.Peptide instances to a
    FASTA-format file.

    :param sequence: The sequence to write.
    :type sequence: coral.DNA, coral.ssDNA, or coral.RNA
    :param handle: File handle (i.e. the output of open()).
    :type handle: file

    '''
    lines = []
    lines.append('>' + sequence.name)
    lines += textwrap.wrap(str(sequence), 79)
    handle.write('\n'.join(lines))
