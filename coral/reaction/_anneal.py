'''Primer Annealing Event Simulation.'''
import coral


class PrimerLengthError(Exception):
    '''Primer does not meet minimum length requirements.'''


class GenericAnnealError(Exception):
    '''Generic primer annealing error.'''


def anneal(template, primer, min_tm=50.0, min_bases=14):
    '''Simulates a primer binding event. Will find the maximum subset
    of bases in the primer that binds to the template. Implicitly finds
    overhang sequences. **Note**: Primer binding locations are now indicative
    of the 3' END of the primer instead of the begining of the annealing
    sequence.

    :param temlate: DNA template for which to bind a primer
    :type template: coral.DNA
    :param primer: primer to bind to template
    :type primer: coral.Primer
    :param min_tm: the minimum required temperature for primer binding
    :type min_tm: float
    :param min_bases: minimum number of bases allowed for binding
    :type min_bases: int
    :returns: dictionary of primer binding locations for each strands, e.g.
              {location: (3' end), annealing sequence, overhang sequence}
    :rtype: int, coral.DNA
    :raises: Exception if primer length is too small
             Exception if primer does not bind
             Exception if primer bind
    '''
    if not type(primer) == coral.Primer:
        raise GenericAnnealError('Primer must be of type coral.Primer')

    if not type(template) == coral.DNA:
        raise GenericAnnealError('Template must be of type coral.DNA')

    # TODO: add possibility for primer basepair mismatch
    # TODO: add
    if len(primer) < min_bases:
        raise PrimerLengthError('Primer length for primer seqeunce \
            {} does not exceed minimum number of bases {}'.format(primer,
                                                                  min_bases))

    fwd_matches = {}
    rev_matches = {}
    for i in range(len(primer) - min_bases + 1)[::-1]:
        primer_dna = primer.overhang + primer.anneal
        annealing = primer_dna[i:]
        overhang = primer_dna[:i]
        anneal_temp = annealing.tm()
        if anneal_temp > min_tm:
            p_matches = template.locate(annealing)
            for match in p_matches[0]:
                fwd_matches[match + len(annealing)] = (primer, annealing,
                                                       overhang)
            for match in p_matches[1]:
                rev_matches[match + len(annealing)] = (primer, annealing,
                                                       overhang)
    return fwd_matches, rev_matches
