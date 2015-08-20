"""Resection (need a new name!) - exonuclease activity."""


def five_resect(dna, n_bases):
    """Remove bases from 5' end of top strand.

    :param dna: Sequence to resect.
    :type dna: pymbt.DNA
    :param n_bases: Number of bases cut back.
    :type n_bases: int
    :returns: DNA sequence resected at the 5' end by n_bases.
    :rtype: pymbt.DNA

    """
    new_instance = dna.copy()
    new_top = '-' * min(len(dna), n_bases) + str(dna)[n_bases:]
    new_instance._sequence = new_top
    new_instance = _remove_end_gaps(new_instance)
    if n_bases >= len(dna):
        new_instance._sequence = ''.join(['-' for i in range(len(dna))])
        new_instance.stranded = 'ss'
    return new_instance


def three_resect(dna, n_bases):
    """Remove bases from 3' end of top strand.

    :param dna: Sequence to resect.
    :type dna: pymbt.DNA
    :param n_bases: Number of bases cut back.
    :type n_bases: int
    :returns: DNA sequence resected at the 3' end by n_bases.
    :rtype: pymbt.DNA

    """
    new_instance = dna.copy()

    new_top = str(dna)[:-n_bases] + '-' * min(len(dna), n_bases)
    new_instance._sequence = new_top
    new_instance = _remove_end_gaps(new_instance)
    if n_bases >= len(dna):
        new_instance._sequence = ''.join(['-' for i in range(len(dna))])
        new_instance.stranded = 'ss'
    return new_instance


def _remove_end_gaps(sequence):
    """Removes double-stranded gaps from ends of the sequence.

    :returns: The current sequence wiht terminal double-strand gaps ('-')
              removed.
    :rtype: pymbt.DNA

    """
    top = sequence.top()
    bottom_rev = sequence.bottom()[::-1]

    top_lstrip = top.lstrip('-')
    bottom_lstrip = bottom_rev.lstrip('-')
    lstrip_len = max(len(top_lstrip), len(bottom_lstrip))
    top = top[-lstrip_len:]
    bottom_rev = bottom_rev[-lstrip_len:]

    top_rstrip = top.rstrip('-')
    bottom_rstrip = bottom_rev.rstrip('-')
    rstrip_len = max(len(top_rstrip), len(bottom_rstrip))
    top = top[0:rstrip_len]
    bottom_rev = bottom_rev[0:rstrip_len]

    bottom = bottom_rev[::-1]
    sequence._sequence = top
    sequence._bottom = bottom
    return sequence
