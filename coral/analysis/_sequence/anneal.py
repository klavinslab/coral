'''Primer Annealing Event Simulation.'''


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
    :returns: A length 2 list (top and bottom strands) of matches. Each
              match is itself a 2-tuple indicating (1) the location on the
              template of the 3' end of the primer binding site and (2) the
              length of the match (number of bases), e.g. [[25, 15],[]] would
              indicate a single top-strand match at template position 25 with
              15 bases of 3' primer homology.
    :rtype: list
    :raises: PrimerLengthError if primer length is too small.
             AnnealError if inputs are of the wrong type.

    '''
    # TODO: add possibility for primer basepair mismatch
    if len(primer) < min_bases:
        msg = 'Primer match length does not exceed minimum number of bases.'
        raise PrimerLengthError(msg)

    # Overwriting dictionary keys ensures uniqueness
    fwd_matches = {}
    rev_matches = {}
    for i in range(len(primer) - min_bases + 1)[::-1]:
        primer_dna = primer.overhang + primer.anneal
        annealing = primer_dna[i:]
        anneal_temp = annealing.tm()
        anneal_len = len(annealing)
        if anneal_temp > min_tm:
            p_matches = template.locate(annealing)
            for match in p_matches[0]:
                fwd_matches[match + anneal_len] = anneal_len
            for match in p_matches[1]:
                rev_matches[match + anneal_len] = anneal_len

    # Convert dictionaries to lists
    fwds = [[key, val] for key, val in fwd_matches.iteritems()]
    revs = [[key, val] for key, val in rev_matches.iteritems()]

    return [fwds, revs]
