import numpy as np


class SubstitutionMatrix(object):
    '''Container for a subsitution matrix - is simply a numpy ndarray with an
    additional .alphabet attribute indicating the labels for each row and
    column.'''

    def __init__(self, matrix, alphabet):
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
        self.matrix = matrix
        self.alphabet = alphabet
