'''PCR reaction(s).'''


class AmbiguousPrimingError(Exception):
    """Primer binds to more than one place on a template."""


class PrimerBindError(Exception):
    """Primer did not bind correctly."""

def primer_binding(template, primer, min_tm=50.0, min_bases=14):
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

    if len(primer) < min_bases:
        raise Exception("Primer length does not exceed minimum number of bases {}".format(min_bases))
    for i in range(len(primer)-min_bases+1):
        anneal = primer.anneal[i:]
        anneal_temp = anneal.tm()
        p_matches = template.locate(anneal)
        p_bind = sum([len(match) for match in p_matches])
        if p_bind > 0:
            if anneal_temp < min_tm:
                raise Exception("Primer binds but does not melting temperature is too low.")
            overhang = primer.anneal[:i]
            fwd_matches = []
            rev_matches = []
            for m in p_matches[0]:
                fwd_matches.append((m + len(anneal), primer, anneal, overhang))
            for m in p_matches[1]:
                rev_matches.append((m + len(anneal), primer, anneal, overhang))
            return fwd_matches, rev_matches
    return [], []

def pcr(template, primer1, primer2, min_tm=50.0, min_bases=14):
    '''Simulate a PCR (no support for ambiguous PCRs).

    :param template: DNA template from which to PCR.
    :type template: coral.DNA
    :param primer1: First PCR primer.
    :type primer1: coral.Primer
    :param primer2: First PCR primer.
    :type primer2: coral.Primer
    :returns: A dsDNA Amplicon.
    :rtype: coral.DNA
    :raises: Exception if a primer binds more than once on the template.
             Exception if primers bind in overlapping sequence of the template.
             Exception if the PCR would work on a circular version of the
             Exception if there are no forward primer binding sites
             Exception if there are no reverse primer binding sites
             template (implies that input was linear).

    '''
    # FIXME: using the wrong primers/template produces a useless error.
    # make the error useful!
    # Find match in top or bottom strands for each primer
    p1_matches = primer_binding(template, primer1, min_tm=min_tm, min_bases=min_bases)
    p2_matches = primer_binding(template, primer2, min_tm=min_tm, min_bases=min_bases)

    # HEY FIX THIS - should find an ambiguity
    # Make sure there's no ambiguities
    p1_bind = sum([len(match) for match in p1_matches])
    p2_bind = sum([len(match) for match in p2_matches])

    def msg(location):
        return "Top strand: {}, Bottom strand: {}".format(location[0],
                                                          location[1])

    if p1_bind > 1:
        raise AmbiguousPrimingError("Primer 1, {}".format(msg(p1_matches)))
    if p2_bind > 1:
        raise AmbiguousPrimingError("Primer 2, {}".format(msg(p2_matches)))
    if not p1_bind and not p2_bind:
        raise PrimerBindError("Neither primer binds the template")
    if not p1_bind:
        raise PrimerBindError("Primer 1 does not bind the template")
    if not p2_bind:
        raise PrimerBindError("Primer 2 does not bind the template")

    # Make 'reverse' index useful for slicing
    fwds = p1_matches[0] + p2_matches[0]
    revs = p2_matches[1] + p1_matches[1]
    if fwds == []:
        raise PrimerBindError('No forward primers found to bind to template')
    if revs == []:
        raise PrimerBindError('No reverse primers found to bind to template')
    fwd = fwds[0]
    rev = revs[0]
    fwd_loc = fwd[0]
    rev_loc = len(template) - rev[0]
    # TODO: circular search will muck things up. If primers are at the very
    # beginning or end of the plasmid coordinates, things will get weird
    # TODO: Should actually just evaluate the length of the product prior
    # to adding primer overhangs, compare to length of anneal regions.
    #   But would need to keep track of sign - fwd - rev != rev - fwd
    '''Notes on PCR amplification decisions:
        If rev == fwd, primers should ampify entire plasmid
        If rev - fwd >= max(len(rev), len(fwd)), amplify sequence normally
        If rev - fwd < max(len(rev), len(fwd)), raise exception - who knows
        how this construct will amplify

    '''
# primer 1 doesn't necessarily always correspond to being the forward primer
    # Subset
    if rev_loc > fwd_loc:
        amplicon = template[fwd_loc:rev_loc]
    elif rev_loc < fwd_loc and rev_loc > fwd_loc - len(fwd[1]):
        raise PrimerBindError('Primers overlap, no solution is implemented')
    else:
        if template.topology == 'circular':
            amplicon = template[fwd_loc:] + template[:rev_loc]
        else:
            raise Exception('Primers would amplify if template were circular.')
#    if primer1.overhang:
#        amplicon = primer1.overhang.to_ds() + amplicon
#    if primer2.overhang:
#        amplicon += primer2.overhang.to_ds().reverse_complement()
    amplicon = fwd[1].primer().to_ds() + amplicon + rev[1].primer().reverse_complement().to_ds()
    return amplicon
