import sys
import numpy as np

# Threshold the adjacency matrix
def thresh_mat(adjacency_matrix):
    # threshold at r > 0.25
    adjacency_matrix[adjacency_matrix < 0.25] = 0
    # fill diagonal with zeroes
    np.fill_diagonal(adjacency_matrix, 0)
    return adjacency_matrix

if __name__ == "__main__":
    cadjacency_matrix = sys.argv[1]
    thresh_mat(adjacency_matrix)


