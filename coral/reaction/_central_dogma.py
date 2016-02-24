'''The central dogma of biology - transcription and translation.'''
import coral as cr


def transcribe(dna):
    '''Transcribe DNA to RNA (no post-transcriptional processing).

    :param seq: Sequence to transcribe (DNA).
    :type seq: coral.DNA
    :returns: Transcribed sequence - an RNA sequence.
    :rtype: coral.RNA

    '''
    return cr.RNA(str(dna).replace('T', 'U'))


def translate(rna):
    '''Translate RNA to peptide.

    :param rna: Sequence to translate (RNA).
    :type rna: coral.RNA
    :returns: Translated sequence - a peptide.
    :rtype: coral.Peptide

    '''
    # Translate
    seq_list = list(str(rna))
    # Convert to peptide until stop codon is found.
    converted = []
    while True:
        if len(seq_list) >= 3:
            base_1 = seq_list.pop(0)
            base_2 = seq_list.pop(0)
            base_3 = seq_list.pop(0)
            codon = ''.join(base_1 + base_2 + base_3).upper()
            amino_acid = cr.constants.codons.CODONS[codon]
            # Stop when stop codon is found
            if amino_acid == '*':
                break
            converted.append(amino_acid)
        else:
            break
    converted = ''.join(converted)
    converted = cr.Peptide(converted)

    return converted


def reverse_transcribe(rna):
    '''Reverse transcribe RNA to DNA.

    :param rna: Sequence to reverse transcribe (RNA).
    :type rna: coral.RNA
    :returns: Reverse-transcribed sequence - a DNA sequence.
    :rtype: coral.DNA

    '''
    return cr.DNA(str(rna).replace('U', 'T'))


def coding_sequence(sequence):
    '''Extract coding sequence from an RNA or DNA template.

    :param sequence: Sequence from which to extract a coding sequence.
    :type sequence: coral.RNA
    :param material: Type of sequence ('dna' or 'rna')
    :type material: str
    :returns: The first coding sequence (start codon -> stop codon) matched
              from 5' to 3'.
    :rtype: coral.RNA
    :raises: ValueError if rna argument has no start codon.
             ValueError if rna argument has no stop codon in-frame with the
             first start codon.

    '''
    if sequence.material == 'rna':
        rna = sequence
    elif sequence.material == 'dna':
        rna = transcribe(sequence)
    else:
        raise ValueError('Sequence must be coral RNA or DNA sequence.')

    codons_left = len(rna) // 3
    start_codon = cr.RNA('aug')
    stop_codons = [cr.RNA('uag'), cr.RNA('uga'), cr.RNA('uaa')]
    start = None
    stop = None
    valid = [None, None]
    index = 0
    while codons_left:
        codon = rna[index:index + 3]
        if valid[0] is None:
            if codon in start_codon:
                start = index
                valid[0] = True
        else:
            if codon in stop_codons:
                stop = index + 3
                valid[1] = True
                break
        index += 3
        codons_left -= 1

    if valid[0] is None:
        raise ValueError('Sequence has no start codon.')
    elif stop is None:
        raise ValueError('Sequence has no stop codon.')
    coding_rna = rna[start:stop]

    return coding_rna
