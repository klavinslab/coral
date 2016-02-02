'''PCR reaction(s).'''
import coral


class PrimingError(Exception):
    '''Raise if there is an error during priming (e.g. during PCR).'''


def pcr(template, primer1, primer2, min_tm=50.0, min_primer_len=14):
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
    :param min_primer_len: Minimum amount of template homology required at the
                           3' end of each primer.
    :type min_primer_len: int
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
                                       min_len=min_primer_len)
    p2_matches = coral.analysis.anneal(template, primer2, min_tm=min_tm,
                                       min_len=min_primer_len)
    p1_binding_locations = [m[0] for strand in p1_matches for m in strand]
    p2_binding_locations = [m[0] for strand in p2_matches for m in strand]

    # Ensure unique top and bottom matches
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

    # Check that primers bind on opposite strands of the template
    tops = p1_matches[0] + p2_matches[0]
    bottoms = p1_matches[1] + p2_matches[1]
    if not tops:
        raise PrimingError('No primers bind the template\'s top strand.')
    if not bottoms:
        raise PrimingError('No primers bind the template\'s bottom strand.')

    # Figure out which primer matches the top strand
    if p1_matches[0]:
        # primer1 is top
        fwd = primer1
        rev = primer2
    else:
        # primer2 matches top strand
        fwd = primer2
        rev = primer1

    # Now we can simulate the PCR. If primer locations are overlapping, we
    # throw an error. If the primers won't amplify a product (e.g. a linear
    # template with primers facing away from one another), throw a different
    # error. Otherwise, amplify the product, including any overhangs.

    # 3' locations, annealing region length
    fwd_3, fwd_len = tops[0]
    rev_3, rev_len = bottoms[0]

    # 5' locations
    fwd_5 = fwd_3 - fwd_len
    rev_5 = rev_3 - rev_len

    # location on the top strand where the 'reverse' primer ends (its 3' end)
    rev_3_top = len(template) - rev_3
    rev_5_top = len(template) - rev_5
    # TODO: Use % operator?
    if rev_5_top > len(template):
        rev_5_top = rev_5_top - len(template)

    # overhangs
    fwd_overhang = fwd.primer()[:-fwd_len]
    rev_overhang = rev.primer()[:-rev_len]

    # TODO: what about searching substrings over circulate templates?
    # Cases:
    # 1)  Primers point towards one another - overlapping is fine
    #       -> isolate region between 5' annealing regions and tack on the
    #          rest of the overhang.
    # 2)  Primers point away from one another, template is linear
    #       -> error
    # 3)  Primers point away from one another, template is circular
    #       a) Primers don't overlap
    #           -> rotate to 'top' primer start, extract sequence
    #       b) Primers overlap
    #           -> Extract whole sequence as linear fragment, tack on rest of
    #              'bottom' primer. May disrupt features.
    if template.circular:
        # Circular template - primers always point towards one another
        if rev_3_top > fwd_3:
            # Inter-primer region doesn't go over the origin (index 0)
            # However, the 'anneal' region may extend over it.
            # FIXME: handle case where 'anneal' region extends over origin
            # FIXME: simplify - just generate 'before' and 'after', then
            # combine with preamplicon later
            if rev_3_top + rev_len > len(template):
                # Reverse primer extends over the origin
                if fwd_5 - fwd_len < 0:
                    # Forward primer extends over the origin
                    preamplicon = template.linearize()
                    # Add extra anneal regions
                    before = template[fwd_5:]
                    after = template[:rev_5_top]
                    preamplicon = before + preamplicon + after
                else:
                    # Only the reverse primer extends over the origin
                    preamplicon = template[fwd_5:]
                    after = template[:rev_5_top]
                    preamplicon = preamplicon + after
            elif fwd_5 - fwd_len < 0:
                # Only the forward primer extends over the origin
                before = template[fwd_5:]
                preamplicon = before + template[:rev_5_top]
            else:
                # Extract like normal
                preamplicon = template[fwd_5:len(template) - rev_5]
        else:
            # Inter-primer region goes over the origin (index 0)
            preamplicon_len = len(template) - fwd_5 + rev_5_top
            preamplicon = template.rotate(-fwd_5)[:preamplicon_len]
    else:
        # Linear template
        if rev_3_top < fwd_5 or fwd_3 > rev_5_top:
            # Primers point away from one another.
            raise PrimingError('Primers point away from one another.')
        else:
            # Primers point towards each other.
            preamplicon = template[fwd_5:len(template) - rev_5]

    # Add overhangs
    amplicon = (fwd_overhang.to_ds() +
                preamplicon +
                rev_overhang.to_ds().reverse_complement())

    return amplicon
