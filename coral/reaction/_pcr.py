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

    fwd_loc, fwd_len = tops[0]
    rev_loc, rev_len = bottoms[0]

    # 5' locations
    fwd_5 = fwd_loc - fwd_len
    rev_5 = rev_loc - rev_len

    # FIXME: what about searching substrings over circulate templates?
    if (len(template) - rev_loc) - fwd_loc <= 0:
        # The primers are either overlapping or are non-overlapping and
        # pointing away from one another
        if fwd_loc < rev_5 and rev_loc > fwd_5:
            # They can overlap and amplify one another
            msg = 'Primer dimer amplification unimplemented.'
            raise NotImplementedError(msg)
        else:
            # They point away from each other, which can only be allowed
            # if the template is circular.
            if template.topology != 'circular':
                raise PrimingError('Primer extension in opposite directions.')

    # All of the checks have passed, so the amplicon can be generated!
    if fwd_loc > (len(template) - rev_5):
        # The primers point away from each other and the template is circular
        # - amplify over the origin
        # TODO: clarify code - len(template) - fwd_loc isn't obvious
        middle_len = (len(template) - fwd_loc) + (len(template) - rev_loc)
        middle = template.rotate(fwd_loc)[:middle_len]
    else:
        middle = template[fwd_loc:len(template) - rev_loc]

    amplicon = (fwd.primer().to_ds() + middle +
                rev.primer().to_ds().reverse_complement())

    return amplicon
