'''Primer Annealing Event Simulation'''

class PrimerLengthError(Exception):
    """Primer does not meet minimum length requirements."""

class PrimerTmError(Exception):
    """Primer annealing temperature is too low."""

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
    :returns: primer binding locations (3' end) and the overhang sequence
    :rtype: int, coral.DNA
    :raises: Exception if primer length is too small
             Exception if primer does not bind
             Exception if primer bind
    '''

    ### TODO: add possibility for primer basepair mismatch
    if len(primer) < min_bases:
        raise PrimerLengthError("Primer length for primer seqeunce \
            {} does not exceed minimum number of bases {}".format(primer, min_bases))
    for i in range(len(primer)-min_bases+1):
        primer_dna = primer.overhang + primer.anneal
        anneal = primer_dna[i:]
        anneal_temp = anneal.tm()
        p_matches = template.locate(anneal)
        p_bind = sum([len(match) for match in p_matches])
        if p_bind > 0:
            if anneal_temp < min_tm:
                raise PrimerTmError("Primer binds but \
                        melting temperature is too low. \
                        Calculated anneal tempearture \
                        {} does not exceed {}".format(anneal_temp, min_tm))
            overhang = primer.anneal[:i]
            fwd_matches = []
            rev_matches = []
            for m in p_matches[0]:
                fwd_matches.append((m + len(anneal), primer, anneal, overhang))
            for m in p_matches[1]:
                rev_matches.append((m + len(anneal), primer, anneal, overhang))
            return fwd_matches, rev_matches
    return [], []
