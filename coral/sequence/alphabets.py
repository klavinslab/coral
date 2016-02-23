'''Defines alphabest available for sequences.'''


class AlphabetError(ValueError):
    '''Raise if sequence doesn\'t match the alphabet.'''


dna_unambiguous = 'ATGC'
dna = dna_unambiguous + 'N-'
dna_ry = dna + 'RY'
dna_double = dna_ry + 'RY'
dna_triple = dna_double + 'BDHV'

rna_unambiguous = 'AUGC'
rna = rna_unambiguous + 'N-'

peptide = 'ACDEFGHIKLMNPQRSTVWYX'
peptide_ambiguous = peptide + 'X*.'
