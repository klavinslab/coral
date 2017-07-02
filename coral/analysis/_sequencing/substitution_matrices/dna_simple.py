import numpy as np
from .substitution_matrix import SubstitutionMatrix


DNA_SIMPLE = SubstitutionMatrix(
  np.array([[1, -1, -1, -1, -1],
            [-1, 1, -1, -1, -1],
            [-1, -1, 1, -1, -1],
            [-1, -1, -1, 1, -1],
            [-1, -1, -1, -1, -1]]),
  'ATGCN')
