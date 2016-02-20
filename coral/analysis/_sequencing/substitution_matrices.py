import numpy as np


class SubstitutionMatrix(np.ndarray):
    '''Generates and contains a substitution matrix indexed by the ASCII number
    of each character (e.g. 'A' is 65 in ASCII and 'T' is 84, so the
    substitution log odds of replacing an 'A' by a 'T' would be located at
    mat[65, 84].

    The (human-readable) inputs can be accessed with SubtitutionMatrix.alphabet
    and SubstitutionMatrix.matrix The ASCII number-indexed matrix is available
    at SubstitutionMatrix.ord_matrix.'''
    def __new__(cls, matrix, alphabet):
        '''
        :param matrix: A square 2D array-like (e.g. list of lists of numbers or
                       numpy array) corresponding to substitution matrix values
                       (often log-odds, as in BLOSUM62). Note: the
                       SubstitutionMatrix will use use integer values, so any
                       floats will be rounded (standard substitution matrices
                       use integers to speed up calculations).
        :type matrix: 2D array-like: np.ndarray or list of lists of integers
        :param alphabet: A string of the characters to use, in order, for the
                         square subsitution matrix (e.g. 'ATGCN' for a simple
                         5X5 DNA substitution matrix).
        :type alphabet: str

        '''
        ords = [ord(char) for char in alphabet]
        mat_size = max(ords) + 1

        template_mat = np.zeros((mat_size, mat_size), dtype=np.int)
        obj = np.asarray(template_mat).view(cls)

        for i, row_ord in enumerate(ords):
            for j, col_ord in enumerate(ords):
                obj[row_ord, col_ord] = matrix[i, j]

        obj.input_alphabet = alphabet
        obj.input_matrix = matrix

        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.input_alphabet = getattr(obj, 'input_alphabet', None)
        self.input_matrix = getattr(obj, 'input_matrix', None)


DNA_SIMPLE = SubstitutionMatrix(np.array([[1, -1, -1, -1, -1],
                                          [-1, 1, -1, -1, -1],
                                          [-1, -1, 1, -1, -1],
                                          [-1, -1, -1, 1, -1],
                                          [-1, -1, -1, -1, -1]]),
                                'ATGCN')

DNA = SubstitutionMatrix(
    np.array([[5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2, -4],
              [-4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, 5],
              [-4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2, -4],
              [-4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2, -4],
              [-4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1,
               -4],
              [1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1, 1],
              [1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1,
               -4],
              [-4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1, 1],
              [-4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1, 1],
              [1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1,
               -4],
              [-4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1,
               -1],
              [-1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1,
               -4],
              [-1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1,
               -1],
              [-1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1,
               -1],
              [-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
               -2],
              [-4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, 5]]),
    'ATGCSWRYKMBVHDNU')


BLOSUM62 = SubstitutionMatrix(
    np.array([[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1,
               1, 0, -3, -2, 0, -2, -1, -1, -1, -4],
              [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2,
               -1, -1, -3, -2, -3, -1, -2, 0, -1, -4],
              [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1,
               0, -4, -2, -3, 4, -3, 0, -1, -4],
              [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1,
               0, -1, -4, -3, -3, 4, -3, 1, -1, -4],
              [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2,
               -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4],
              [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0,
               -1, -2, -1, -2, 0, -2, 4, -1, -4],
              [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0,
               -1, -3, -2, -2, 1, -3, 4, -1, -4],
              [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2,
               0, -2, -2, -3, -3, -1, -4, -2, -1, -4],
              [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2,
               -1, -2, -2, 2, -3, 0, -3, 0, -1, -4],
              [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3,
               -2, -1, -3, -1, 3, -3, 3, -3, -1, -4],
              [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3,
               -2, -1, -2, -1, 1, -4, 3, -3, -1, -4],
              [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1,
               0, -1, -3, -2, -2, 0, -3, 1, -1, -4],
              [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2,
               -1, -1, -1, -1, 1, -3, 2, -1, -1, -4],
              [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4,
               -2, -2, 1, 3, -1, -3, 0, -3, -1, -4],
              [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,
               7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4],
              [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4,
               1, -3, -2, -2, 0, -2, 0, -1, -4],
              [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2,
               -1, 1, 5, -2, -2, 0, -1, -1, -1, -1, -4],
              [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1,
               -4, -3, -2, 11, 2, -3, -4, -2, -2, -1, -4],
              [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3,
               -3, -2, -2, 2, 7, -1, -3, -1, -2, -1, -4],
              [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2,
               -2, 0, -3, -1, 4, -3, 2, -2, -1, -4],
              [-2, -1, 4, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0,
               -1, -4, -3, -3, 4, -3, 0, -1, -4],
              [-1, -2, -3, -3, -1, -2, -3, -4, -3, 3, 3, -3, 2, 0, -3,
               -2, -1, -2, -1, 2, -3, 3, -3, -1, -4],
              [-1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0,
               -1, -2, -2, -2, 0, -3, 4, -1, -4],
              [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4],
              [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
               -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]]),
    'ARNDCQEGHILKMFPSTWYVBJZX*')
