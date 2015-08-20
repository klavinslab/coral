'''Check sequences for repeats that may impact cloning efficiency.'''
from collections import Counter


def repeats(seq, size):
    '''Count times that a sequence of a certain size is repeated.

    :param seq: Input sequence.
    :type seq: pymbt.DNA or pymbt.RNA
    :param size: Size of the repeat to count.
    :type size: int
    :returns: Occurrences of repeats and how many
    :rtype: tuple of the matched sequence and how many times it occurs

    '''
    seq = str(seq)
    n_mers = [seq[i:i + size] for i in range(len(seq) - size + 1)]
    counted = Counter(n_mers)
    # No one cares about patterns that appear once, so exclude them
    found_repeats = [(key, value) for key, value in counted.iteritems() if
                     value > 1]
    return found_repeats
