'''Primer Annealing Event Simulation.'''
import coral


class PrimerLengthError(Exception):
    '''Primer does not meet minimum length requirements.'''


class AnnealError(Exception):
    '''Generic primer annealing error.'''


def anneal(template, primer, min_tm=50.0, min_bases=14):
    '''Simulates a primer binding event. Will find the maximum subset
    of bases in the primer that binds to the template, including overhang
    sequences. **Note**: Primer binding locations indicate the 3' end of the
    primer, not the begining of the annealing sequence.

    :param template: DNA template for which to bind a primer.
    :type template: coral.DNA
    :param primer: Primer to bind to template.
    :type primer: coral.Primer
    :param min_tm: The cutoff melting temperature for primer binding - a binder
                   with a lower Tm will be rejected.
    :type min_tm: float
    :param min_bases: The cutoff for bases required for binding - a binder with
                      fewer bases will be rejected.
    :type min_bases: int
    :returns: A 2-tuple of dictionaries representing matches on the top and
              bottom strands. Each dictionary has keys of this format:
              {location: (annealing sequence, overhang sequence)},
              where location the integer indicating the location of the 3'
              end of the primer bound to the template, the annealing sequence
              is a primer of the annealing sequence, and the overhang sequence
              is a coral.DNA instance.
    :rtype: tuple
    :raises: PrimerLengthError if primer length is too small.
             AnnealError if inputs are of the wrong type.

    '''
    # FIXME: this violates duck typing
    if not type(primer) == coral.Primer:
        raise AnnealError('Primer must be of type coral.Primer')

    if not type(template) == coral.DNA:
        raise AnnealError('Template must be of type coral.DNA')

    # TODO: add possibility for primer basepair mismatch
    if len(primer) < min_bases:
        msg = 'Primer match length does not exceed minimum number of bases.'
        raise PrimerLengthError(msg)

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
    return (fwd_matches, rev_matches)
