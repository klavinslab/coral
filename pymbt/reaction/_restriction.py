'''Restriction endonuclease reactions.'''
import pymbt.reaction


def digest(dna, restriction_enzyme):
    '''Restriction endonuclease reaction.

    :param dna: DNA template to digest.
    :type dna: pymbt.DNA
    :param restriction_site: Restriction site to use.
    :type restriction_site: RestrictionSite
    :returns: list of digested DNA fragments.
    :rtype: pymbt.DNA list

    '''
    pattern = restriction_enzyme.recognition_site
    located = dna.locate(pattern)
    if not located[0] and not located[1]:
        return [dna]
    # Bottom strand indices are relative to the bottom strand 5' end.
    # Convert to same type as top strand
    pattern_len = len(pattern)
    r_indices = [len(dna) - index - pattern_len for index in
                 located[1]]
    # If sequence is palindrome, remove redundant results
    if pattern.is_palindrome():
        r_indices = [index for index in r_indices if index not in
                     located[0]]
    # Flatten cut site indices
    cut_sites = sorted(located[0] + r_indices)
    # Go through each cut site starting at highest one
    # Cut remaining template once, generating remaining + new
    current = [dna]
    for cut_site in cut_sites[::-1]:
        new = _cut(current, cut_site, restriction_enzyme)
        current.append(new[1])
        current.append(new[0])
    current.reverse()
    # Combine first and last back together if digest was circular
    if dna.topology == 'circular':
        current[0] = current.pop() + current[0]
    return current


def _cut(dna, index, restriction_enzyme):
    '''Cuts template once at the specified index.

    :param dna: DNA to cut
    :type dna: pymbt.DNA
    :param index: index at which to cut
    :type index: int
    :param restriction_enzyme: Enzyme with which to cut
    :type restriction_enzyme: pymbt.RestrictionSite
    :returns: 2-element list of digested sequence, including any overhangs.
    :rtype: list

    '''
    # TODO: handle case where cut site is outside of recognition sequence,
    # for both circular and linear cases where site is at index 0
    # Find absolute indices at which to cut
    cut_site = restriction_enzyme.cut_site
    top_cut = index + cut_site[0]
    bottom_cut = index + cut_site[1]

    # Isolate left and ride sequences
    to_cut = dna.pop()
    max_cut = max(top_cut, bottom_cut)
    min_cut = min(top_cut, bottom_cut)
    left = to_cut[:max_cut]
    right = to_cut[min_cut:]

    # If applicable, leave overhangs
    diff = top_cut - bottom_cut
    if not diff:
        # Blunt-end cutter, no adjustment necessary
        pass
    elif diff > 0:
        # 3' overhangs
        left = pymbt.reaction.five_resect(left.flip(), diff).flip()
        right = pymbt.reaction.five_resect(right, diff)
    else:
        # 5' overhangs
        left = pymbt.reaction.three_resect(left, abs(diff))
        right = pymbt.reaction.three_resect(right.flip(), abs(diff)).flip()

    return [left, right]
