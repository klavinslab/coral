'''Utilities for reactions.'''
import pymbt
from pymbt.constants.molecular_bio import ALPHABETS, CODONS


def convert_sequence(seq, to_material):
    '''Translate a DNA sequence into peptide sequence.

    The following conversions are supported:
        Transcription (seq is DNA, to_material is 'rna')
        Reverse transcription (seq is RNA, to_material is 'dna')
        Translation (seq is RNA, to_material is 'peptide')

    :param seq: DNA or RNA sequence.
    :type seq: pymbt.DNA or pymbt.RNA
    :param to_material: material to which to convert ('rna', 'dna', or
                        'peptide').
    :type to_material: str
    :returns: sequence of type pymbt.sequence.[material type]

    '''
    if isinstance(seq, pymbt.DNA) and to_material == 'rna':
        # Transcribe

        # Can't transcribe a gap
        if '-' in seq:
            raise ValueError('Cannot transcribe gapped DNA')
        # Convert DNA chars to RNA chars
        origin = ALPHABETS['dna'][:-1]
        destination = ALPHABETS['rna']
        code = dict(zip(origin, destination))
        converted = ''.join([code.get(str(k), str(k)) for k in seq])
        # Instantiate RNA object
        converted = pymbt.RNA(converted)
    elif isinstance(seq, pymbt.RNA):
        if to_material == 'dna':
            # Reverse transcribe
            origin = ALPHABETS['rna']
            destination = ALPHABETS['dna'][:-1]
            code = dict(zip(origin, destination))
            converted = ''.join([code.get(str(k), str(k)) for k in seq])
            # Instantiate DNA object
            converted = pymbt.DNA(converted)
        elif to_material == 'peptide':
            # Translate
            seq_list = list(str(seq))
            # Convert to peptide until stop codon is found.
            converted = []
            while True:
                if len(seq_list) >= 3:
                    base_1 = seq_list.pop(0)
                    base_2 = seq_list.pop(0)
                    base_3 = seq_list.pop(0)
                    codon = ''.join(base_1 + base_2 + base_3).upper()
                    amino_acid = CODONS[codon]
                    # Stop when stop codon is found
                    if amino_acid == '*':
                        break
                    converted.append(amino_acid)
                else:
                    break
            converted = ''.join(converted)
            converted = pymbt.Peptide(converted)
    else:
        msg1 = 'Conversion from '
        msg2 = '{0} to {1} is not supported.'.format(seq.__class__.__name__,
                                                     to_material)
        raise ValueError(msg1 + msg2)

    return converted
