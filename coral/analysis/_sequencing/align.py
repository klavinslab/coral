'''Numpy implementation of Needlman-Wunsch algorithm'''
import numpy as np
from . import substitution_matrices as submat


def as_ord_matrix(matrix):
    '''Given the SubstitutionMatrix input, generate an equivalent matrix that
    is indexed by the ASCII number of each residue (e.g. A -> 65).'''
    ords = [ord(c) for c in matrix.alphabet]
    ord_matrix = np.zeros((max(ords) + 1, max(ords) + 1), dtype=np.integer)
    for i, row_ord in enumerate(ords):
        for j, col_ord in enumerate(ords):
            ord_matrix[row_ord, col_ord] = matrix[i, j]

    return ord_matrix


def max_index(array):
    '''Locate the index of the largest value in the array. If there are
    multiple, finds the earliest one in the row-flattened array.

    :param array: Any array.
    :type array: numpy.array

    '''
    return np.unravel_index(array.argmax(), array.shape)


def aligner(seqj, seqi, method='global', gap_open=-7, gap_extend=-7,
            gap_double=-7, matrix=submat.DNA_SIMPLE):
    '''Calculates the alignment of two sequences. The global method uses
    a global Needleman-Wunsh algorithm, local does a a local
    Smith-Waterman alignment, global_cfe does a global alignment with
    cost-free ends and glocal does an alignment which is global only with
    respect to the shorter sequence, also known as a semi-global
    alignment. Returns the aligned (sub)sequences as character arrays.

    Gotoh, O. (1982). J. Mol. Biol. 162, 705-708.
    Needleman, S. & Wunsch, C. (1970). J. Mol. Biol. 48(3), 443-53.
    Smith, T.F. & Waterman M.S. (1981). J. Mol. Biol. 147, 195-197.

    :param seqj: First sequence.
    :type seqj: str
    :param seqi: Second sequence.
    :type seqi: str
    :param method: Type of alignment: 'global', 'global_cfe', 'local', or
                   'glocal'.
    :type method: str
    :param gap_open: The cost of opening a gap (negative number).
    :type gap_open: float
    :param gap_extend: The cost of extending an open gap (negative number).
    :type gap_extend: float
    :param gap_double: The gap-opening cost if a gap is already open in the
                       other sequence (negative number).
    :type gap_double: float
    :param matrix: A score matrix dictionary name. Examples can be found in
                   the substitution_matrices module.
    :type matrix: str

    '''
    amatrix = as_ord_matrix(matrix)
    NONE, LEFT, UP, DIAG = range(4)  # NONE is 0
    max_j = len(seqj)
    max_i = len(seqi)

    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i
    else:
        flip = 0

    F = np.zeros((max_i + 1, max_j + 1), dtype=np.float32)
    I = np.ndarray((max_i + 1, max_j + 1), dtype=np.float32)
    I.fill(-np.inf)
    J = np.ndarray((max_i + 1, max_j + 1), dtype=np.float32)
    J.fill(-np.inf)
    pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)  # NONE

    if method == 'global':
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
        F[0, 1:] = gap_open + gap_extend * np.arange(0, max_j,
                                                     dtype=np.float32)
        F[1:, 0] = gap_open + gap_extend * np.arange(0, max_i,
                                                     dtype=np.float32)
    elif method == 'global_cfe':
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
    elif method == 'glocal':
        pointer[0, 1:] = LEFT
        F[0, 1:] = gap_open + gap_extend * np.arange(0, max_j,
                                                     dtype=np.float32)

    seqi_ord = [ord(base) for base in seqi]
    seqj_ord = [ord(base) for base in seqj]
    for i in range(1, max_i + 1):
        ci = seqi_ord[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj_ord[j - 1]
            # I
            I[i, j] = max(F[i, j - 1] + gap_open,
                          I[i, j - 1] + gap_extend,
                          J[i, j - 1] + gap_double)
            # J
            J[i, j] = max(F[i - 1, j] + gap_open,
                          J[i - 1, j] + gap_extend,
                          I[i - 1, j] + gap_double)
            # F
            diag_score = F[i - 1, j - 1] + amatrix[ci, cj]
            left_score = I[i, j]
            up_score = J[i, j]
            max_score = max(diag_score, up_score, left_score)

            F[i, j] = max(0, max_score) if method == 'local' else max_score

            if method == 'local':
                if F[i, j] == 0:
                    pass  # point[i,j] = NONE
                elif max_score == diag_score:
                    pointer[i, j] = DIAG
                elif max_score == up_score:
                    pointer[i, j] = UP
                elif max_score == left_score:
                    pointer[i, j] = LEFT
            elif method == 'glocal':
                # In a semi-global alignment we want to consume as much as
                # possible of the longer sequence.
                if max_score == up_score:
                    pointer[i, j] = UP
                elif max_score == diag_score:
                    pointer[i, j] = DIAG
                elif max_score == left_score:
                    pointer[i, j] = LEFT
            else:
                # global
                if max_score == up_score:
                    pointer[i, j] = UP
                elif max_score == left_score:
                    pointer[i, j] = LEFT
                else:
                    pointer[i, j] = DIAG

    align_j = []
    align_i = []
    if method == 'local':
        # max anywhere
        i, j = max_index(F)
    elif method == 'glocal':
        # max in last col
        i, j = (F[:, -1].argmax(), max_j)
    elif method == 'global_cfe':
        # from i,j to max(max(last row), max(last col)) for free
        row_max, col_idx = F[-1].max(), F[-1].argmax()
        col_max, row_idx = F[:, -1].max(), F[:, -1].argmax()
        if row_max > col_max:
            pointer[-1, col_idx + 1:] = LEFT
        else:
            pointer[row_idx + 1:, -1] = UP

    p = pointer[i, j]
    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1
            align_j.append(seqj[j])
            align_i.append(seqi[i])
        elif p == LEFT:
            j -= 1
            align_j.append(seqj[j])
            align_i.append('-')
        elif p == UP:
            i -= 1
            align_j.append('-')
            align_i.append(seqi[i])
        else:
            raise Exception('wtf!')
        p = pointer[i, j]
    align_i = ''.join(align_i[::-1])
    align_j = ''.join(align_j[::-1])
    # np.array(align_i.reverse())
    return ((align_i, align_j) if flip else (align_j, align_i))


def score_alignment(a, b, gap_open, gap_extend, matrix):
    '''Calculate the alignment score from two aligned sequences.

    :param a: The first aligned sequence.
    :type a: str
    :param b: The second aligned sequence.
    :type b: str
    :param gap_open: The cost of opening a gap (negative number).
    :type gap_open: int
    :param gap_extend: The cost of extending an open gap (negative number).
    :type gap_extend: int.
    :param matrix: A score matrix dictionary name. Examples can be found in
                   the substitution_matrices module.

    '''
    al = a
    bl = b
    l = len(al)
    score = 0
    assert len(bl) == l, 'Alignment lengths must be the same'
    mat = as_ord_matrix(matrix)

    gap_started = 0

    for i in range(l):
        if al[i] == '-' or bl[i] == '-':
            score += gap_extend if gap_started else gap_open
            gap_started = 1
        else:
            score += mat[ord(al[i]), ord(bl[i])]
            gap_started = 0
    return score
