'''PCR reaction(s).'''
import coral


class PrimingError(Exception):
    '''Raise if there is an error during priming (e.g. during PCR).'''


def pcr(template, primer1, primer2, min_tm=50.0, min_bases=14):
    '''Simulate a PCR.

    :param template: DNA template from which to PCR.
    :type template: coral.DNA
    :param primer1: First PCR primer.
    :type primer1: coral.Primer
    :param primer2: First PCR primer.
    :type primer2: coral.Primer
    :param min_tm: Minimum melting temperature (Tm) at which primers must bind
                   to the template.
    :type min_tm: float
    :param min_bases: Minimum amount of template homology required at the 3'
                      end of each primer.
    :type min_bases: int
    :returns: A dsDNA Amplicon.
    :rtype: coral.DNA
    :raises: PrimingError if a primer binds more than once on the template,
             primers bind in overlapping sequence of the template, there are no
             forward primer binding sites or reverse priming sites, or if the
             PCR would work only on a circular version of the template (if
             template is linear).

    '''
    # Find match in top or bottom strands for each primer
    p1_matches = coral.analysis.anneal(template, primer1, min_tm=min_tm,
                                       min_bases=min_bases)
    p2_matches = coral.analysis.anneal(template, primer2, min_tm=min_tm,
                                       min_bases=min_bases)
    p1_binding_locations = p1_matches[0].keys() + p1_matches[1].keys()
    p2_binding_locations = p2_matches[0].keys() + p2_matches[1].keys()

    # FIXME: Make sure there's no ambiguities

    if len(p1_binding_locations) > 1:
        primer_msg = 'Multiple primer 1 binding locations: {}'
        raise PrimingError(primer_msg.format(p1_binding_locations))
    if len(p2_binding_locations) > 1:
        primer_msg = 'Multiple primer 2 binding locations: {}'
        raise PrimingError(primer_msg.format(p2_binding_locations))
    if not p1_binding_locations and not p2_binding_locations:
        raise PrimingError('Neither primer binds the template')
    if not p1_binding_locations:
        raise PrimingError('Primer 1 does not bind the template')
    if not p2_binding_locations:
        raise PrimingError('Primer 2 does not bind the template')

    # Now see if the primers bind on opposite strands of the template
    fwds = p1_matches[0].copy()
    fwds.update(p2_matches[0])
    revs = p1_matches[1].copy()
    revs.update(p2_matches[1])
    if not fwds:
        raise PrimingError('No primers bind the template\'s top strand.')
    if not revs:
        raise PrimingError('No primers bind the template\'s bottom strand.')

    fwd = fwds[fwds.keys()[0]]
    rev = revs[revs.keys()[0]]
    fwd_loc = fwds.keys()[0]
    rev_loc = len(template) - revs.keys()[0]
    # TODO: circular search will muck things up. If primers are at the very
    # beginning or end of the plasmid coordinates, things will get weird
    # TODO: Should actually just evaluate the length of the product prior
    # to adding primer overhangs, compare to length of anneal regions.
    #   But would need to keep track of sign - fwd - rev != rev - fwd
    '''
    Notes on PCR amplification decisions:
        If rev == fwd, primers should ampify entire plasmid
        If rev - fwd >= max(len(rev), len(fwd)), amplify sequence normally
        If rev - fwd < max(len(rev), len(fwd)), raise exception - who knows
        how this construct will amplify
    '''
    # Note: primer 1 doesn't necessarily always correspond to being the forward
    # primer.
    fwd_primer = fwd[0]
    rev_primer = rev[0]
    if rev_loc < fwd_loc and rev_loc > fwd_loc - len(fwd_primer):
        raise PrimingError('Primers overlap, no solution is implemented')

    if rev_loc > fwd_loc:
        amplicon = template[fwd_loc:rev_loc]
    else:
        if template.topology == 'circular':
            amplicon = template[fwd_loc:] + template[:rev_loc]
        else:
            raise Exception('Primers would amplify only if template were \
                            circular.')
    amplicon = (fwd_primer.primer().to_ds() + amplicon +
                rev_primer.primer().reverse_complement().to_ds())

    return amplicon
