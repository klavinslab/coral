import numpy as np
cimport numpy as np
from . import substitution_matrices as submat
from libc.string cimport strlen


# Access to the Python/C API
cdef extern from 'Python.h':
    ctypedef void PyObject
    PyObject *PyString_FromStringAndSize(char *, size_t)
    int _PyString_Resize(PyObject **, size_t)
    char * PyString_AS_STRING(PyObject *)


# Declaring numpy data types speeds things up massively
ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.float32_t DTYPE_FLOAT


cdef inline DTYPE_FLOAT max3(DTYPE_FLOAT a, DTYPE_FLOAT b, DTYPE_FLOAT c):
    '''Find largest of 3 floats. Much faster than using built-in max.

    :param a: First number.
    :type a: DTYPE_FLOAT
    :param b: Second number.
    :type b: DTYPE_FLOAT
    :param c: Third number.
    :type c: DTYPE_FLOAT

    '''
    if c > b:
        return c if c > a else a
    return b if b > a else a


cdef inline DTYPE_FLOAT max2(DTYPE_FLOAT a, DTYPE_FLOAT b):
    '''Find largest of 2 floats. Much faster than using built-in max and max3.

    :param a: First number.
    :type a: DTYPE_FLOAT
    :param b: Second number.
    :type b: DTYPE_FLOAT

    '''
    return b if b > a else a


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


def aligner(_seqj, _seqi, DTYPE_FLOAT gap_open=-7, DTYPE_FLOAT gap_extend=-7,
            DTYPE_FLOAT gap_double=-7, method='global',
            matrix=submat.DNA_SIMPLE):
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
    :type matrix: SubstitutionMatrix

    '''
    cdef int NONE = 0,  LEFT = 1, UP = 2,  DIAG = 3
    cdef bint flip = 0
    cdef char* seqj = _seqj
    cdef char* seqi = _seqi
    cdef size_t align_counter = 0

    cdef int imethod

    if method == 'global':
        imethod = 0
    elif method == 'local':
        imethod = 1
    elif method == 'glocal':
        imethod = 2
    elif method == 'global_cfe':
        imethod = 3

    cdef size_t max_j = strlen(seqj)
    cdef size_t max_i = strlen(seqi)
    if max_i == max_j == 0:
        return '', ''

    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i

    cdef char *align_j, *align_i
    cdef int i, j
    cdef char ci, cj
    cdef PyObject *ai, *aj

    assert gap_extend <= 0, 'gap_extend penalty must be <= 0'
    assert gap_open <= 0, 'gap_open must be <= 0'

    cdef np.ndarray[DTYPE_FLOAT, ndim=2] agap_i = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
    cdef np.ndarray[DTYPE_FLOAT, ndim=2] agap_j = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
    agap_i.fill(-np.inf)
    agap_j.fill(-np.inf)

    cdef np.ndarray[DTYPE_FLOAT, ndim=2] score = np.zeros((max_i + 1, max_j + 1), dtype=np.float32)

    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)
    cdef np.ndarray[DTYPE_INT, ndim=2] amatrix = matrix
    amatrix = as_ord_matrix(matrix)

    # START HERE:
    if imethod == 0:
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
        score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)
        score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.float32)
    elif imethod == 3:
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
    elif imethod == 2:
        pointer[0, 1:] = LEFT
        score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)

    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]
            # agap_i
            agap_i[i,j] = max3(
                         score[i, j - 1] + gap_open,
                         agap_i[i, j - 1] + gap_extend,
                         agap_j[i, j - 1] + gap_double)
            # agap_j
            agap_j[i,j] = max3(
                         score[i - 1, j] + gap_open,
                         agap_j[i - 1, j] + gap_extend,
                         agap_i[i - 1, j] + gap_double)
            # score
            diag_score = score[i - 1, j - 1] + amatrix[ci, cj]
            left_score = agap_i[i, j]
            up_score   = agap_j[i, j]
            max_score = max3(diag_score, up_score, left_score)

            score[i, j] = max_score

            # global
            if max_score == up_score:
                pointer[i,j] = UP
            elif max_score == left_score:
                pointer[i,j] = LEFT
            else:
                pointer[i,j] = DIAG


    if imethod == 0:
        # max anywhere
        i, j = max_index(score)
    elif imethod == 2:
        # max in last col
        i, j = (score[:,-1].argmax(), max_j)
    elif imethod == 3:
        # from i,j to max(max(last row), max(last col)) for free
        row_max, col_idx = score[-1].max(), score[-1].argmax()
        col_max, row_idx = score[:, -1].max(), score[:, -1].argmax()
        if row_max > col_max:
            pointer[-1,col_idx+1:] = LEFT
        else:
            pointer[row_idx+1:,-1] = UP

    seqlen = max_i + max_j
    ai = PyString_FromStringAndSize(NULL, seqlen)
    aj = PyString_FromStringAndSize(NULL, seqlen)

    # use this and PyObject instead of assigning directly...
    align_j = PyString_AS_STRING(aj)
    align_i = PyString_AS_STRING(ai)

    p = pointer[i, j]
    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = seqi[i]
        elif p == LEFT:
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = c'-'
        elif p == UP:
            i -= 1
            align_j[align_counter] = c'-'
            align_i[align_counter] = seqi[i]
        else:
            raise Exception('wtf!:pointer: %i', p)
        align_counter += 1
        p = pointer[i, j]

    _PyString_Resize(&aj, align_counter)
    _PyString_Resize(&ai, align_counter)

    if flip:
        return (<object>ai)[::-1], (<object>aj)[::-1]
    else:
        return (<object>aj)[::-1], (<object>ai)[::-1]


def score_alignment(a, b, int gap_open, int gap_extend, matrix):
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
    :type matrix: SubstitutionMatrix

    '''
    cdef char *al = a
    cdef char *bl = b
    cdef size_t l = strlen(al), i
    cdef int score = 0, this_score
    assert strlen(bl) == l, 'Alignment lengths must be the same'
    cdef np.ndarray[DTYPE_INT, ndim=2] mat
    mat = as_ord_matrix(matrix)

    cdef bint gap_started = 0

    for i in range(l):
        if al[i] == c'-' or bl[i] == c'-':
            score += gap_extend if gap_started else gap_open
            gap_started = 1
        else:
            score += mat[al[i], bl[i]]
            gap_started = 0
    return score
