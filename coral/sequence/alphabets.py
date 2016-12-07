'''Defines alphabest available for sequences.'''


class AlphabetError(ValueError):
    '''Raise if sequence doesn\'t match the alphabet.'''


# If you know the complements, you know the alphabet
# Sometimes, the complements don't apply, so you need just an alphabet. This
# prevents the 'complement' method from ever running. i.e. we should always
# know the state of the alphabet and the complements.

class Alphabet(object):
    '''Container for biological sequence alphabets - supports automatic
    creation of valid characters from a sequence-complement dictionary (e.g.
    A:T pairs).

    '''
    def __init__(self, complements=None, symbols=None):
        '''
        :param complements: Dict of valid complement sequences, e.g.
                            {'A': 'T'}. If not supplied, user must indicate
                            an alphabet.
        :type complements: dict
        :param symbols: A string of valid characters the biological sequence
                        can be made up of. Case and order doesn't matter.
        :type symbols: str

        '''
        if complements is None and symbols is None:
            raise ValueError('Initialized without complements or symbols.')
        elif complements is None:
            # User put in symbols, but not complements
            pass
        elif symbols is None:
            # User put in complements, but not symbols. Generate by flattening
            # list
            symbols = ''.join(set([y for x in complements.items() for y in x]))

        # The value of 'complements' can end up being None, but the symbols
        # attribute will always be not-None (at minimum, an empty string)

        self.complements = complements
        self.symbols = symbols

    def __repr__(self):
        return 'complements: {}\nsymbols: {}'.format(self.complements,
                                                     self.symbols)


def _merge_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z


# Unambiguous DNA: good ol' A, T, G, and C
dna_unambiguous = Alphabet(complements={
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
})

# Default DNA: allows ambiguous N and gaps (-), but nothing else
dna = Alphabet(complements=_merge_dicts(dna_unambiguous.complements, {
    'N': 'N',
    '-': '-'
}))

dna_ry = Alphabet(complements=_merge_dicts(dna.complements, {
    'R': 'Y',
    'Y': 'R'
}))


dna_double = Alphabet(symbols=dna_ry.symbols + 'RrYy')
dna_double = Alphabet(symbols=dna_double.symbols + 'BbDdHhVv')

rna_unambiguous = Alphabet(complements={
    'A': 'U',
    'U': 'A',
    'G': 'C',
    'C': 'G'
})

rna = Alphabet(complements=_merge_dicts(rna_unambiguous.complements, {
    'N': 'N',
    'n': 'n',
    '-': '-'
}))

peptide_unambiguous = Alphabet(symbols='ACDEFGHIKLMNPQRSTVWYX')
peptide = Alphabet(symbols=peptide_unambiguous.symbols + '*')
peptide_ambiguous = Alphabet(symbols=peptide.symbols + 'X.')
