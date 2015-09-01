'''The central dogma of biology - transcription and translation.'''
import coral
from . import utils


def transcribe(dna):
    '''Transcribe DNA to RNA (no post-transcriptional processing).

    :param seq: Sequence to transcribe (DNA).
    :type seq: coral.DNA
    :returns: Transcribed sequence - an RNA sequence.
    :rtype: coral.RNA

    '''
    return utils.convert_sequence(dna, 'rna')


def translate(rna):
    '''Translate RNA to peptide.

    :param rna: Sequence to translate (RNA).
    :type rna: coral.RNA
    :returns: Translated sequence - a peptide.
    :rtype: coral.Peptide

    '''
    return utils.convert_sequence(rna, 'peptide')


def reverse_transcribe(rna):
    '''Reverse transcribe RNA to DNA.

    :param rna: Sequence to reverse transcribe (RNA).
    :type rna: coral.RNA
    :returns: Reverse-transcribed sequence - a DNA sequence.
    :rtype: coral.DNA

    '''
    return utils.convert_sequence(rna, 'dna')


def coding_sequence(rna):
    """Extract coding sequence from an RNA template.

    :param seq: Sequence from which to extract a coding sequence.
    :type seq: coral.RNA
    :param material: Type of sequence ('dna' or 'rna')
    :type material: str
    :returns: The first coding sequence (start codon -> stop codon) matched
              from 5' to 3'.
    :rtype: coral.RNA
    :raises: ValueError if rna argument has no start codon.
             ValueError if rna argument has no stop codon in-frame with the
             first start codon.

    """
    if isinstance(rna, coral.DNA):
        rna = transcribe(rna)
    codons_left = len(rna) // 3
    start_codon = coral.RNA('aug')
    stop_codons = [coral.RNA('uag'), coral.RNA('uga'), coral.RNA('uaa')]
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
