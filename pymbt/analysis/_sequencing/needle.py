'''Needleman-Wunsch alignment functions.'''
import multiprocessing
import pymbt
try:
    from .calign import aligner, score_alignment
except ImportError:
    from .align import aligner, score_alignment


def needle(reference, query, gap_open=-15, gap_extend=0, matrix='DNA_simple'):
    '''Do a Needleman-Wunsch alignment.

    :param reference: Reference sequence.
    :type reference: pymbt.DNA
    :param query: Sequence to align against the reference.
    :type query: pymbt.DNA
    :param gapopen: Penalty for opening a gap.
    :type gapopen: float
    :param gapextend: Penalty for extending a gap.
    :type gapextend: float
    :param matrix: Matrix to use for alignment - options are DNA_simple (for
                   DNA) and BLOSUM62 (for proteins).
    :type matrix: str
    :returns: (aligned reference, aligned query, score)
    :rtype: tuple of two pymbt.DNA instances and a float

    '''
    # Align using cython Needleman-Wunsch
    aligned_ref, aligned_res = aligner(str(reference),
                                       str(query),
                                       gap_open=gap_open,
                                       gap_extend=gap_extend,
                                       method='global_cfe',
                                       matrix=matrix)

    # Score the alignment
    score = score_alignment(aligned_ref, aligned_res, gap_open, gap_extend,
                            'DNA_simple')

    return pymbt.DNA(aligned_ref), pymbt.DNA(aligned_res), score


def run_needle(args):
    '''Run needle command using 4-tuple of the arguments (in the same order)
    as is used for needle. Necessary to make picklable function for
    multiprocessing.'''
    return needle(*args)


def needle_multi(references, queries, gap_open=-15, gap_extend=0,
                 matrix='DNA_simple'):
    '''Batch process of sequencing split over several cores. Acts just like
    needle but sequence inputs are lists.

    :param references: References sequence.
    :type references: pymbt.DNA list
    :param queries: Sequences to align against the reference.
    :type queries: pymbt.DNA list
    :param gap_open: Penalty for opening a gap.
    :type gap_open: float
    :param gap_extend: Penalty for extending a gap.
    :type gap_extend: float
    :param matrix: Matrix to use for alignment - options are DNA_simple (for
                   DNA) and BLOSUM62 (for proteins).
    :type matrix: str
    :returns: a list of the same output as pymbt.sequence.needle
    :rtype: list

    '''
    pool = multiprocessing.Pool()
    try:
        args_list = [[ref, que, gap_open, gap_extend, matrix] for ref, que in
                     zip(references, queries)]
        aligned = pool.map(run_needle, args_list)
    except KeyboardInterrupt:
        pool.terminate()
        raise KeyboardInterrupt

    return aligned
