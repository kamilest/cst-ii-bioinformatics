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


def additive_phylogeny(dist_matrix, n):
    if n == 2:
        return ({0: [1], 1: [0]}, {(1,0): dist_matrix[1, 0], (0, 1): dist_matrix[0, 1]})
    limb = limb_length(dist_matrix, n-1)
    dist_matrix[:-1, -1] -= limb
    dist_matrix[-1, :-1] -= limb
    
    # Find three leaves such that Di,k = Di,n + Dn,k
    i_leaf, k_leaf = None, None
    for i in range(n-2):
        for k in range(i+1, n-1):
            if dist_matrix[i, k] == \
                    dist_matrix[i, n-1] + dist_matrix[k, n-1]:
                i_leaf, k_leaf = i, k
                break
    x = dist_matrix[i_leaf, n-1]

    # Recursive call
    (T_neighbours, T_weights) = \
        additive_phylogeny(dist_matrix, n-1)

    # Potentially new node in T at distance x from i
    # on the path between i_leaf and k_leaf

    # new node v
    v = min(len(dist_matrix), max(T_neighbours.keys())+1)
    T_neighbours[v] = []
    T_neighbours[i_leaf].append(v)
    T_neighbours[k_leaf].append(v)
    T_neighbours[v].append(i_leaf)
    T_neighbours[v].append(k_leaf)
    T_weights[(v, i_leaf)] = T_weights[(i_leaf, v)] = x
    T_weights[(v, k_leaf)] = T_weights[(k_leaf, v)] = dist_matrix[i_leaf, k_leaf] - x

    # Add leaf n back to T by creating limb (v, n)
    # of length limb
    T_neighbours[v].append(n-1)
    T_neighbours[n-1] = [v]
    T_weights[(n, v)] = T_weights[(v, n)] = limb
    return (T_neighbours, T_weights)


f = open('ba7c.txt', 'r')
n = int(f.readline())

dist_matrix = np.zeros((n, n))
for i, line in enumerate(f):
    values = [int(m) for m in line.strip().split()]
    dist_matrix[i, :] = values

(T_neighbours, T_weights) = additive_phylogeny(dist_matrix, n)

# Print tree
for (v, w) in sorted(T_weights.keys()):
    print("{}->{}:{}".format(v, w, T_weights[(v, w)]))
