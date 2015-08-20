'''Gibson design module.'''
import pymbt


class LengthError(Exception):
    '''If primer would be longer than max length, throw this exception'''
    pass


class TmError(Exception):
    '''If the assembly overlap would fall below the minimum Tm, throw this.'''
    pass


def gibson_primers(dna1, dna2, overlap='mixed', maxlen=80, overlap_tm=65.0,
                   insert=None, primer_kwargs={}):
    '''Design Gibson primers given two DNA sequences (connect left to right)

    :param dna1: First piece of DNA for which to design primers. Once Gibsoned,
                 would be connected at its right side to dna2.
    :type dna1: pymbt.DNA
    :param dna2: First piece of DNA for which to design primers. Once Gibsoned,
                 would be connected at its right side to dna2.
    :type dna2: pymbt.DNA
    :param overlap: Specifies location of overlap. 'left' puts it on the 'dna1'
                    side (i.e. the primer to amplify dna2). 'right' puts it on
                    the dna2 side, and 'mixed' does a ~50:50 split
    :type overlap: str
    :param maxlen: Maximum length of each primer.
    :type maxlen: int
    :param overlap_tm: Minimum Tm of overlap
    :type overlap_tm: float
    :param insert: A DNA insert to add with primers and use as assembly
                   homology. This overrides the 'split' argument.
    :type insert: pymbt.DNA
    :param primer_kwargs: keyword arguments to pass to design_primer()
    :type primer_kwargs: dict
    :returns: Reverse, then forward primer for bridging the two sequences.
              Note that the forward primer binds dna2, reverse dna1.
    :rtype: A sequence.Primer tuple
    :raises: ValueError if split parameter is an invalid string.

    '''
    # Annealing sequences
    # DNA 2 primer is a forward primer
    fwd_anneal = pymbt.design.primer(dna2, **primer_kwargs)
    # DNA 1 primer is a reverse primer
    rev_anneal = pymbt.design.primer(dna1.flip(), **primer_kwargs)
    # Overhangs
    if insert is None:
        # No insert, so follow split argument
        if overlap == 'left':
            # If splitting left, put overhang on forward primer
            overlap_revcomp = pymbt.design.primer(dna1.flip(), tm=overlap_tm,
                                                  tm_undershoot=0)
            fwd_overhang = overlap_revcomp.primer().reverse_complement()
            rev_overhang = None
        elif overlap == 'right':
            # If splitting right, put overhang on reverse primer
            overlap = pymbt.design.primer(dna2, tm=overlap_tm, tm_undershoot=0)
            fwd_overhang = None
            rev_overhang = overlap.primer().reverse_complement()
        elif overlap == 'mixed':
            # If mixed, grow size of both until overlap Tm is reached
            overlap_l = dna1[0:0]  # Empty sequence.DNA
            overlap_r = dna2[0]  # First base
            overlap_melt = pymbt.analysis.tm(overlap_r)  # Initial overlap Tm
            while overlap_melt < overlap_tm:
                rlen = len(overlap_r)
                llen = len(overlap_l)
                if rlen > llen:
                    # Increase left side of overlap
                    overlap_l = dna1[-(rlen + 1):]
                else:
                    # Increase right side of overlap
                    overlap_r = dna2[:(llen + 1)]
                overlap = overlap_l + overlap_r
                overlap_melt = pymbt.analysis.tm(overlap)
            fwd_overhang = overlap_l
            rev_overhang = overlap_r.reverse_complement()
        else:
            raise ValueError('split argument must be left, right, or mixed')
        # Generate primers using anneal, overhang, and tm data
        fwd = pymbt.Primer(fwd_anneal.primer(), tm=fwd_anneal.tm,
                           overhang=fwd_overhang)
        rev = pymbt.Primer(rev_anneal.primer(), tm=rev_anneal.tm,
                           overhang=rev_overhang)
    else:
        # There's an insert to use as the overhang
        overlap = insert
        fwd_overhang = insert.to_ss()
        rev_overhang = insert.reverse_complement().to_ss()
        # Generate primers using anneal, overhang, and tm data
        fwd = pymbt.Primer(fwd_anneal.primer(), tm=fwd_anneal.tm,
                           overhang=fwd_overhang)
        rev = pymbt.Primer(rev_anneal.primer(), tm=rev_anneal.tm,
                           overhang=rev_overhang)
        left_trim = 0
        # If either primer is too long, try trimming the overhang
        while len(fwd) > maxlen:
            # Generate new overlap
            overlap = insert[left_trim:]
            # Tm must be above overlap_tm
            if pymbt.analysis.tm(overlap) < overlap_tm:
                raise TmError('Right primer is too long with this Tm setting.')
            # Regenerate forward overhang
            fwd_overhang = overlap.to_ss()
            # Regenerate primer with new overhang
            fwd = pymbt.Primer(fwd_anneal.primer(), tm=fwd_anneal.tm,
                               overhang=fwd_overhang)
            # Increase 'trimming' index
            left_trim += 1
        right_trim = 0
        while len(rev) > maxlen:
            # Generate new overlap
            overlap = insert[:len(insert) - right_trim]
            # Tm must be above overlap_tm
            if pymbt.analysis.tm(overlap) < overlap_tm:
                raise TmError('Left primer is too long with this Tm setting.')
            # Regenerate reverse overhang
            rev_overhang = overlap.reverse_complement().to_ss()
            rev = pymbt.Primer(rev_anneal.primer(), tm=rev_anneal.tm,
                               overhang=rev_overhang)
            # Increase 'trimming' index
            right_trim += 1
    # Check primer lengths
    if any([len(primer) > maxlen for primer in (fwd, rev)]):
        raise LengthError('At least one of the primers is longer than maxlen.')

    return rev, fwd


def gibson(seq_list, circular=True, overlaps='mixed', overlap_tm=65,
           maxlen=80, primer_kwargs={}):
    '''Design Gibson primers given a set of sequences

    :param seq_list: List of DNA sequences to stitch together
    :type seq_list: list containing pymbt.DNA
    :param circular: If true, designs primers for making a circular construct.
                     If false, designs primers for a linear construct.
    :type circular: bool
    :param splits: Specifies locations of overlap. Must be either a single
                   entry of the same type as the 'split' parameter in
                   gibson_primers or a list of those types of the appropriate
                   length (for circular construct, len(seq_list), for
                   linear construct, len(seq_list) - 1)
    :type splits: str or list of str
    :param overlap_tm: Minimum Tm of overlap
    :type overlap_tm: float
    :param maxlen: Maximum length of each primer.
    :type maxlen: int
    :param primer_kwargs: keyword arguments to pass to design.primer
    :type primer_kwargs: dict
    :returns: Forward and reverse primers for amplifying every fragment.
    :rtype: a list of sequence.Primer tuples
    :raises: ValueError if split parameter is an invalid string or wrong size.

    '''

    # Input checking
    if circular:
        n_overlaps = len(seq_list)
    else:
        n_overlaps = len(seq_list) - 1

    if type(overlaps) is str:
        overlaps = [overlaps] * n_overlaps
    else:
        if len(overlaps) != n_overlaps:
            raise ValueError('Incorrect number of "overlaps" entries.')
        else:
            for overlap in overlaps:
                if overlap not in ['left', 'right', 'mixed']:
                    raise ValueError('Invalid "overlaps" setting.')

    # If here, inputs were good
    # Design primers for linear constructs:
    primers_list = []
    for i, (left, right) in enumerate(zip(seq_list[:-1], seq_list[1:])):
        primers_list.append(gibson_primers(left, right, overlaps[i],
                                           overlap_tm=overlap_tm,
                                           primer_kwargs=primer_kwargs))
    if circular:
        primers_list.append(gibson_primers(seq_list[-1], seq_list[0],
                                           overlaps[-1],
                                           overlap_tm=overlap_tm,
                                           primer_kwargs=primer_kwargs))
    else:
        primer_f = pymbt.design.primer(seq_list[0], **primer_kwargs)
        primer_r = pymbt.design.primer(seq_list[-1].reverse_complement(),
                                       **primer_kwargs)
        primers_list.append((primer_r, primer_f))

    # Primers are now in order of 'reverse for seq1, forward for seq2' config
    # Should be in 'forward and reverse primers for seq1, then seq2', etc
    # Just need to rotate one to the right
    flat = [y for x in primers_list for y in x]
    flat = [flat[-1]] + flat[:-1]
    grouped_primers = [(flat[2 * i], flat[2 * i + 1]) for i in
                       range(len(flat) / 2)]

    return grouped_primers


def gibson_equimolar(lengths, concs):
    pass
