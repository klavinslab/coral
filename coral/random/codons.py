'''Generate a random DNA sequence.'''
import random
import coral
from coral.constants.codons import CODON_FREQ_BY_AA


def random_codons(peptide, frequency_cutoff=0.0, weighted=False, table=None):
    '''Generate randomized codons given a peptide sequence.

    :param peptide: Peptide sequence for which to generate randomized
                    codons.
    :type peptide: coral.Peptide
    :param frequency_cutoff: Relative codon usage cutoff - codons that
                             are rarer will not be used. Frequency is
                             relative to average over all codons for a
                             given amino acid.
    :param frequency_cutoff: Codon frequency table to use.
    :param weighted: Use codon table
    :type weighted: bool
    :param table: Codon frequency table to use. Table should be organized
                  by amino acid, then be a dict of codon: frequency.
                  Only relevant if weighted=True or frequency_cutoff > 0.
                  Tables available:

                  constants.molecular_bio.CODON_FREQ_BY_AA['sc'] (default)
    :type table: dict
    :returns: Randomized sequence of codons (DNA) that code for the input
              peptide.
    :rtype: coral.DNA
    :raises: ValueError if frequency_cutoff is set so high that there are no
             codons available for an amino acid in the input peptide.

    '''
    if table is None:
        table = CODON_FREQ_BY_AA['sc']

    # Process codon table using frequency_cutoff
    def _cutoff(table, cutoff):
        '''Generate new codon frequency table given a mean cutoff.

        :param table: codon frequency table of form
                      {amino acid: codon: frequency}.
        :type table: dict
        :param cutoff: value between 0 and 1.0 for mean frequency cutoff.
        :type cutoff: float
        :returns: A codon frequency table with some codons removed.
        :rtype: dict

        '''
        new_table = {}
        # IDEA: cutoff should be relative to most-frequent codon, not average?
        for amino_acid, codons in table.iteritems():
            average_cutoff = cutoff * sum(codons.values()) / len(codons)
            new_table[amino_acid] = {}
            for codon, frequency in codons.iteritems():
                if frequency > average_cutoff:
                    new_table[amino_acid][codon] = frequency
        return new_table

    new_table = _cutoff(table, frequency_cutoff)
    # Select codons randomly or using weighted distribution
    rna = ''
    for amino_acid in str(peptide):
        codons = new_table[amino_acid.upper()]
        if not codons:
            raise ValueError('No {} codons at freq cutoff'.format(amino_acid))
        if weighted:
            cumsum = []
            running_sum = 0
            for codon, frequency in codons.iteritems():
                running_sum += frequency
                cumsum.append(running_sum)
            random_num = random.uniform(0, max(cumsum))
            for codon, value in zip(codons, cumsum):
                if value > random_num:
                    selection = codon
                    break
        else:
            selection = random.choice(codons.keys())
        rna += selection
    return coral.RNA(rna)
