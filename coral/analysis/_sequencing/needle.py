'''Needleman-Wunsch alignment functions.'''
import coral as cr
from . import substitution_matrices as submat
import multiprocessing
import warnings
try:
    from .calign import aligner, score_alignment
except ImportError:
    message = ('NW alignment extension could not be imported, falling back'
               'on native Python version (~100 times slower).')
    warnings.warn(message)
    from .align import aligner, score_alignment


def needle(reference, query, gap_open=-15, gap_extend=0,
           matrix=submat.DNA_SIMPLE):
    '''Do a Needleman-Wunsch alignment.

    :param reference: Reference sequence.
    :type reference: coral.DNA
    :param query: Sequence to align against the reference.
    :type query: coral.DNA
    :param gapopen: Penalty for opening a gap.
    :type gapopen: float
    :param gapextend: Penalty for extending a gap.
    :type gapextend: float
    :param matrix: Matrix to use for alignment - options are DNA_simple (for
                   DNA) and BLOSUM62 (for proteins).
    :type matrix: str
    :returns: (aligned reference, aligned query, score)
    :rtype: tuple of two coral.DNA instances and a float

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
                            matrix)

    return cr.DNA(aligned_ref), cr.DNA(aligned_res), score


def run_needle(args):
    '''Run needle command using 4-tuple of the arguments (in the same order)
    as is used for needle. Necessary to make picklable function for
    multiprocessing.'''
    return needle(*args)


def needle_msa(reference, results, gap_open=-15, gap_extend=0,
               matrix=submat.DNA_SIMPLE):
    '''Create a multiple sequence alignment based on aligning every result
    sequence against the reference, then inserting gaps until every aligned
    reference is identical

    '''
    gap = '-'
    # Convert alignments to list of strings
    alignments = []
    for result in results:
        ref_dna, res_dna, score = needle(reference, result, gap_open=gap_open,
                                         gap_extend=gap_extend, matrix=matrix)
        alignments.append([str(ref_dna), str(res_dna), score])

    def insert_gap(sequence, position):
        return sequence[:position] + gap + sequence[position:]

    i = 0
    while True:
        # Iterate over 'columns' in every reference
        refs = [alignment[0][i] for alignment in alignments]

        # If there's a non-unanimous gap, insert gap into alignments
        gaps = [ref == gap for ref in refs]
        if any(gaps) and not all(gaps):
            for alignment in alignments:
                if alignment[0][i] != gap:
                    alignment[0] = insert_gap(alignment[0], i)
                    alignment[1] = insert_gap(alignment[1], i)

        # If all references match, we're all done
        alignment_set = set(alignment[0] for alignment in alignments)
        if len(alignment_set) == 1:
            break

        # If we've reach the end of some, but not all sequences, add end gap
        lens = [len(alignment[0]) for alignment in alignments]
        if i + 1 in lens:
            for alignment in alignments:
                if len(alignment[0]) == i + 1:
                    alignment[0] = alignment[0] + gap
                    alignment[1] = alignment[1] + gap

        i += 1

        if i > 20:
            break

    # Convert into MSA format
    output_alignment = [cr.DNA(alignments[0][0])]
    for alignment in alignments:
        output_alignment.append(cr.DNA(alignment[1]))

    return output_alignment


def needle_multi(references, queries, gap_open=-15, gap_extend=0,
                 matrix=submat.DNA_SIMPLE):
    '''Batch process of sequencing split over several cores. Acts just like
    needle but sequence inputs are lists.

    :param references: References sequence.
    :type references: coral.DNA list
    :param queries: Sequences to align against the reference.
    :type queries: coral.DNA list
    :param gap_open: Penalty for opening a gap.
    :type gap_open: float
    :param gap_extend: Penalty for extending a gap.
    :type gap_extend: float
    :param matrix: Matrix to use for alignment - options are DNA_simple (for
                   DNA) and BLOSUM62 (for proteins).
    :type matrix: str
    :returns: a list of the same output as coral.sequence.needle
    :rtype: list

    '''
    pool = multiprocessing.Pool()
    try:
        args_list = [[ref, que, gap_open, gap_extend, matrix] for ref, que in
                     zip(references, queries)]
        aligned = pool.map(run_needle, args_list)
    except KeyboardInterrupt:
        print 'Caught KeyboardInterrupt, terminating workers'
        pool.terminate()
        pool.join()
        raise KeyboardInterrupt

    return aligned
