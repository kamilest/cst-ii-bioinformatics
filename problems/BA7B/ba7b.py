import numpy as np

def limb_length(dist_matrix, leaf):
    leaf_row = dist_matrix[leaf, :]
    
    length = np.inf
    for i in range(len(leaf_row)):
        for j in range(i+1, len(leaf_row)):
            if i != leaf and j != leaf:
                length = min(length,\
                    (dist_matrix[leaf, i] +\
                     dist_matrix[leaf, j] -\
                     dist_matrix[i, j]) / 2.0)

    return int(length)


f = open('rosalind_ba7b.txt', 'r')
n = int(f.readline())
j = int(f.readline())

dist_matrix = np.zeros((n, n))
for i, line in enumerate(f):
    values = [int(m) for m in line.strip().split()]
    dist_matrix[i, :] = values

length = limb_length(dist_matrix, j)
print(length)
