"""Resection (need a new name!) - exonuclease activity."""
import coral


def five_resect(dna, n_bases):
    """Remove bases from 5' end of top strand.

    :param dna: Sequence to resect.
    :type dna: coral.DNA
    :param n_bases: Number of bases cut back.
    :type n_bases: int
    :returns: DNA sequence resected at the 5' end by n_bases.
    :rtype: coral.DNA

    """
    new_instance = dna.copy()
    if n_bases >= len(dna):
        new_instance._top.seq = ''.join(['-' for i in range(len(dna))])
    else:
        new_instance._top.seq = '-' * n_bases + str(dna)[n_bases:]

    new_instance = _remove_end_gaps(new_instance)

    return new_instance


def three_resect(dna, n_bases):
    """Remove bases from 3' end of top strand.

    :param dna: Sequence to resect.
    :type dna: coral.DNA
    :param n_bases: Number of bases cut back.
    :type n_bases: int
    :returns: DNA sequence resected at the 3' end by n_bases.
    :rtype: coral.DNA

    """
    new_instance = dna.copy()
    if n_bases >= len(dna):
        new_instance._top.seq = ''.join(['-' for i in range(len(dna))])
    else:
        new_instance._top.seq = str(dna)[:-n_bases] + '-' * n_bases

    new_instance = _remove_end_gaps(new_instance)

    return new_instance


def _remove_end_gaps(sequence):
    """Removes double-stranded gaps from ends of the sequence.

    :returns: The current sequence with terminal double-strand gaps ('-')
              removed.
    :rtype: coral.DNA

    """
    # Count terminal blank sequences
    def count_end_gaps(seq):
        gap = coral.DNA('-')
        count = 0
        for base in seq:
            if base == gap:
                count += 1
            else:
                break

        return count

    top_left = count_end_gaps(sequence.top())
    top_right = count_end_gaps(reversed(sequence.top()))
    bottom_left = count_end_gaps(reversed(sequence.bottom()))
    bottom_right = count_end_gaps(sequence.bottom())

    # Trim sequence
    left_index = min(top_left, bottom_left)
    right_index = len(sequence) - min(top_right, bottom_right)

    return sequence[left_index:right_index]
