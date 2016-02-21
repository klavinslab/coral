'''Primer Annealing Event Simulation.'''


class PrimerLengthError(Exception):
    '''Primer does not meet minimum length requirements.'''


class AnnealError(Exception):
    '''Generic primer annealing error.'''


def anneal(template, primer, min_tm=50.0, min_len=10):
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
    :param min_len: The cutoff for bases required for binding - a binder with
                    fewer bases will be rejected.
    :type min_len: int
    :returns: A length 2 list (top and bottom strands) of matches. Each
              match is itself a 2-tuple indicating (1) the location on the
              template of the 3' end of the primer binding site and (2) the
              length of the match (number of bases), e.g. [[(25, 15)],[]] would
              indicate a single top-strand match at template position 25 with
              15 bases of 3' primer homology.
    :rtype: list
    :raises: PrimerLengthError if primer length is too small.
             AnnealError if inputs are of the wrong type.

    '''
    # TODO: add possibility for primer basepair mismatch
    if len(primer) < min_len:
        msg = 'Primer length is shorter than min_len argument.'
        raise PrimerLengthError(msg)
    if len(template) < min_len:
        msg = 'Template is shorter than the min_len argument.'
        raise AnnealError(msg)

    # Strategy: locate all min-length matches, then extend them until they
    # no longer match. This provides an advantage over the previous strategy of
    # updating a dictionary with indices from coral.DNA.locate() as keys, as
    # the latter's indices may actually move for a given primer as it passes
    # over the origin
    def update_match_linear(base, location_length, anneal_seq):
        '''Increase the location and length of binding site, if applicable.'''
        # TODO: this is very inefficient - should stop updating once the first
        # mismatch occurs.
        location, length = location_length

        if location == 0:
            return location_length

        location_next = location - 1
        length_next = length + 1
        seq = base[location_next:location_next + length_next]

        if seq == anneal_seq:
            return (location_next, length_next)
        else:
            return location_length

    def update_match_circular(base, location_length, anneal_seq):
        '''Increase the location and length of binding site, if applicable.'''
        # TODO: this is very inefficient - should stop updating once the first
        # mismatch occurs.
        base_len = len(base)
        location, length = location_length
        if location == 0:
            location_next = base_len - 1
        else:
            location_next = location - 1
        length_next = length + 1

        if (location_next + length_next) > base_len:
            upstream = base[location_next:]
            downstream = base[:length_next - (base_len - location_next)]
            seq = upstream + downstream
        else:
            # No need to 'rotate' sequence
            seq = base[location_next:location_next + length_next]

        if seq == anneal_seq:
            return (location_next, length_next)
        else:
            return location_length

    if template.circular:
        update_fun = update_match_circular
    else:
        update_fun = update_match_linear

    # Maximum annealing length to test (can't exceed template length)
    max_len = min(len(template), len(primer))

    primer_dna = primer.to_ds()
    anneal_len = min_len
    anneal_seq = primer_dna[-anneal_len:]
    binding_data = []
    for k, strand_locs in enumerate(template.locate(anneal_seq)):
        matches = zip(strand_locs, [min_len] * len(strand_locs))
        for i in range(anneal_len + 1, max_len + 1):
            anneal_seq = primer_dna[-i:]
            for j, match in enumerate(matches):
                if k == 0:
                    matches[j] = update_fun(template.top, match, anneal_seq)
                else:
                    matches[j] = update_fun(template.bottom, match,
                                            anneal_seq)
        binding_data.append(matches)

    # Now, filter out all the matches that are too short
    for i in reversed(range(len(primer_dna) + 1)):
        min_len = i + 1
        tm = primer_dna[-min_len:].tm()
        if tm < min_tm:
            break

    for strand in binding_data:
        for i in reversed(range(len(strand))):
            if strand[i][1] < min_len:
                strand.pop(i)

    # Finally, adjust the position to be the 3' end
    for strand in binding_data:
        for i, match in enumerate(strand):
            length = match[1]
            loc_new = match[0] + length
            if loc_new > len(template):
                # Circularly permute
                loc_new = loc_new - len(template)
            strand[i] = [loc_new, length]

    # Overwriting dictionary keys ensures uniqueness
    # fwd_matches = {}
    # rev_matches = {}
    # for i in range(len(primer) - min_len + 1)[::-1]:
    #     primer_dna = primer.overhang + primer.anneal
    #     annealing = primer_dna[i:]
    #     anneal_temp = annealing.tm()
    #     anneal_len = len(annealing)
    #     if anneal_temp > min_tm:
    #         p_matches = template.locate(annealing)
    #         for match in p_matches[0]:
    #             fwd_matches[match + anneal_len] = anneal_len
    #         for match in p_matches[1]:
    #             rev_matches[match + anneal_len] = anneal_len

    # # Convert dictionaries to lists
    # fwds = [[key, val] for key, val in fwd_matches.iteritems()]
    # revs = [[key, val] for key, val in rev_matches.iteritems()]

    return binding_data
