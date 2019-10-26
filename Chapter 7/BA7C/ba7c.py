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
        return ({0: [1], 1: [0]}, \
            {(1,0): dist_matrix[1, 0], (0, 1): dist_matrix[0, 1]})
    
    limb = limb_length(dist_matrix, n-1)
    dist_matrix[:-1, -1] -= limb
    dist_matrix[-1, :-1] -= limb
    
    # Find three leaves such that Di,k = Di,n + Dn,k
    i_leaf, k_leaf = None, None
    for k in range(len(dist_matrix)-1):
        # Removing Di,n from all Di,k (Di,k - Di,n)
        remaining_length = dist_matrix[k] - dist_matrix[-1]
        # Find indices where remaining_length = Dn,k
        # Di,k - Di,n = Dn,k
        ix = np.where(remaining_length == dist_matrix[k, -1])
        if len(ix[0] > 0):
            i_leaf, k_leaf = ix[0][0], k
            break

    x = dist_matrix[i_leaf, -1]

    # Recursive call for the distance matrix without outer row/col
    (T_neighbours, T_weights) = \
        additive_phylogeny(dist_matrix[:-1, :-1], n-1)

    # Potentially new node v in T at distance x from i_leaf on the path 
    # between i_leaf and k_leaf

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
