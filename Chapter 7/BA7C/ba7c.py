import numpy as np
from queue import Queue

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

# Inspired by samuelklee/ucsd-bioinformatics at the point I had given up on bugs
def additive_phylogeny(dist_matrix, n):
    def nodes_on_path(T_neighbours, i_leaf, k_leaf, total_n, visited):
        if k_leaf in T_neighbours[i_leaf]:
            return [i_leaf, k_leaf]
        to_visit = list(filter(lambda i: i >= total_n and i not in visited, T_neighbours[i_leaf]))
        if len(to_visit) == 0:
            return []

        visited.add(i_leaf)
        path = [i_leaf]

        for node in to_visit:
            path.extend(nodes_on_path(T_neighbours, node, k_leaf, total_n, visited))
            if path[-1] == k_leaf:
                break
            else:
                path = [i_leaf]

        return path
    
    def add_edge(T_neighbours, T_weights, origin, destination, weight):
        if origin not in T_neighbours.keys():
            T_neighbours[origin] = []
        if destination not in T_neighbours.keys():
            T_neighbours[destination] = []
        
        T_neighbours[origin].append(destination)
        T_neighbours[destination].append(origin)
        T_weights[(origin, destination)] = weight
        T_weights[(destination, origin)] = weight

    # Add leaf with weight limb at distance x from leaf i_keaf on path between leaves i_leaf and k_leaf, potentially creating new v_node
    def add_leaf(T_neighbours, T_weights, i_leaf, k_leaf, x, leaf, limb, total_n):
        path = nodes_on_path(T_neighbours, i_leaf, k_leaf, total_n, set([]))

        distance_remaining = x
        origin = i_leaf

        for destination in path[1:]:
            distance_remaining -= T_weights[(origin, destination)]
            if distance_remaining == 0:
                add_edge(T_neighbours, T_weights, leaf, destination, limb)
                break
            if distance_remaining < 0:
                weight = T_weights[(origin, destination)]
                origin_to_new_node_distance = distance_remaining + weight
                destination_to_new_node_distance = -distance_remaining

                v_ix = np.amax([k for k, v in T_neighbours.items()])
                v_node = total_n if v_ix < total_n else v_ix+1

                # print("Adding ", leaf, " with limb ", limb, " at node "
                # v_node, " between ", origin, " and ", destination)
                T_neighbours[origin].remove(destination)
                T_neighbours[destination].remove(origin)
                del T_weights[(origin, destination)]
                del T_weights[(destination, origin)]

                add_edge(T_neighbours, T_weights, origin, v_node, origin_to_new_node_distance)
                add_edge(T_neighbours, T_weights, destination, v_node, destination_to_new_node_distance)
                add_edge(T_neighbours, T_weights, v_node, leaf, limb)
                break

            origin = destination
        pass

    T_neighbours = {}
    T_weights = {}

    if n == 2:
        origin = 0
        destination = 1
        weight = dist_matrix[(origin, destination)]
        add_edge(T_neighbours, T_weights, origin, destination, weight)
        return (T_neighbours, T_weights)

    limb = limb_length(dist_matrix[:n, :n], n-1)
    dist_matrix[:n-1, n-1] -= limb
    dist_matrix[n-1, :n-1] -= limb
    
    # Find three leaves such that Di,k = Di,n + Dn,k
    i_leaf, k_leaf = None, None
    for k in range(n-1):
        # Dn,k = Di,k - Di,n
        remaining_length = dist_matrix[:n-1, k] - dist_matrix[:n-1, n-1]
        # Find indices where remaining_length = Dn,k
        ix = np.where(remaining_length == dist_matrix[n-1, k])
        for j in ix[0]:
            if j != k:
                i_leaf, k_leaf = j, k
            break

    x = dist_matrix[i_leaf, n-1]
    T_neighbours, T_weights = additive_phylogeny(dist_matrix, n-1)

    # print("Adding node at distance", x, "between", i_leaf, "and", k_leaf)
    add_leaf(T_neighbours, T_weights, i_leaf, k_leaf, x, n-1, limb,\
         len(dist_matrix))
    return (T_neighbours, T_weights)

f = open('rosalind_ba7c.txt', 'r')
n = int(f.readline())

dist_matrix = np.zeros((n, n))
for i, line in enumerate(f):
    values = [int(m) for m in line.strip().split()]
    dist_matrix[i, :] = values

(T_neighbours, T_weights) = additive_phylogeny(dist_matrix, n)

# Print tree
for (v, w) in sorted(T_weights.keys()):
    print("{}->{}:{}".format(v, w, int(T_weights[(v, w)])))
