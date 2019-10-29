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


def additive_phylogeny(dist_matrix, n):
    def nodes_on_path(T_neighbours, i_leaf, k_leaf, n, visited):
        if k_leaf in T_neighbours[i_leaf]:
            return [i_leaf, k_leaf]
        unvisited = np.where(T_neighbours[i_leaf] >= n and T_neighbours[i_leaf]not in visited)
        if len(unvisited) == 0:
            return []

        visited.add(i_leaf)
        path = [i_leaf]

        for node in unvisited:
            node.append(nodes_on_path(T_neighbours, node, k_leaf, n, visited))
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

    # Add leaf with limbLength at distance x from leaf i on path between leaves i and k, creating new internal node if needed
    def add_leaf(T_neighbours, T_weights, i_leaf, k_leaf, x, leaf, limb, n):
        path = nodes_on_path(T_neighbours, i_leaf, k_leaf, n, set([]))
        print("Path: ", path)

        distanceRemaining = x
        origin = i_leaf

        for destination in path[1:]:
            distanceRemaining -= T_weights[(origin, destination)]
            if distanceRemaining == 0:
                add_edge(T_neighbours, T_weights, leaf, destination, limb)
                break;
            if distanceRemaining < 0:
                weight = T_weights[(origin, destination)]
                origin_to_new_node_distance = distanceRemaining + weight
                destination_to_new_node_distance = -distanceRemaining
                v_node = min(n, \
                    np.amax([k for k, v in T_neighbours.items()])+1)

                print("Adding " + leaf + " with limb " + limb + " at node " + v_node + " between " + origin + " and " + destination)
                print()

                T_neighbours[origin].remove(destination)
                T_neighbours[destination].remove(origin)
                del T_weights[(origin, destination)]
                del T_weights[(destination, origin)]

                add_edge(T_neighbours, T_weights, origin, v_node, origin_to_new_node_distance)
                add_edge(T_neighbours, T_weights, destination, v_node, destination_to_new_node_distance)
                add_edge(T_neighbours, T_weights, v_node, leaf, limb)
                break

            origin = destination

    if n == 2:
        return ({0: [1], 1: [0]}, \
            {(1,0): dist_matrix[1, 0], (0, 1): dist_matrix[0, 1]})
    
    limb = limb_length(dist_matrix[:n, :n], n-1)
    dist_matrix[:n-1, n-1] -= limb
    dist_matrix[n-1, :n-1] -= limb
    
    # Find three leaves such that Di,k = Di,n + Dn,k
    i_leaf, k_leaf = None, None
    for k in range(n-1):
        # Dn,k = Di,k - Di,n
        remaining_length = dist_matrix[:n-1, k] - dist_matrix[:n-1, n-1]
        np.delete(remaining_length, k, 0)
        # Find indices where remaining_length = Dn,k
        ix = np.where(remaining_length == dist_matrix[n-1, k])
        if len(ix[0] > 0):
            i_leaf, k_leaf = ix[0][0], k
            break

    print(i_leaf, k_leaf)
    x = dist_matrix[i_leaf, n-1]

    # Recursive call for the distance matrix without outer row/col
    (T_neighbours, T_weights) = additive_phylogeny(dist_matrix, n-1)
    print(T_neighbours)
    print(T_weights)

    # Potentially new node v in T at distance x from i_leaf on the path 
    # between i_leaf and k_leaf
    s = {} # node -> (min_weight, prev_node)
    s[i_leaf] = (0, None)

    next_nodes = Queue()
    next_nodes.put(i_leaf)
    visited = set([])

    while not next_nodes.empty():
        v = next_nodes.get()
        visited.add(v)

        for w in T_neighbours[v]:
            if w != k_leaf and w not in visited:
                next_nodes.put(w)
            
            s[w] = (s[v][0] + T_weights[(v, w)], v)

    assert(k_leaf in s)
    # Backtrack through path from i_leaf to k_leaf
    path = [k_leaf]
    current = k_leaf
    while current is not i_leaf:
        print(current)
        current = s[current][1]
        path.append(current)
    path = path[::-1]
    print(path)


    # Search where is the v_node or where to insert it
    v_node_ix = None
    v_node = None
    for i in range(len(path)-1):
        if s[path[i]][0] < x and s[path[i+1]][0] >= x:
            v_node_ix = i+1
            break
    
    # Create new node or set v_node to existing node
    if s[path[v_node_ix]][0] == x:
        v_node = path[v_node_ix]
    else:
        v_node = min(len(dist_matrix), \
            np.amax([k for k, v in T_neighbours.items()])+1)
        i_rem = path[v_node_ix-1]
        k_rem = path[v_node_ix]
        T_neighbours[i_rem].remove(k_rem)
        T_neighbours[k_rem].remove(i_rem)
        T_neighbours[i_rem].append(v_node)
        T_neighbours[k_rem].append(v_node)
        
        T_weights[(v_node, i_rem)] = x - s[i_rem][0]
        T_weights[(i_rem, v_node)] = x - s[i_rem][0]
        T_weights[(v_node, k_rem)] = s[k_rem][0] - x
        T_weights[(k_rem, v_node)] = s[k_rem][0] - x
        del T_weights[(i_rem, k_rem)]
        del T_weights[(k_rem, i_rem)]

    # new node v
    # T_neighbours[v] = []
    # T_neighbours[i_leaf].append(v)
    # T_neighbours[k_leaf].append(v)
    # T_neighbours[v].append(i_leaf)
    # T_neighbours[v].append(k_leaf)
    # T_weights[(v, i_leaf)] = T_weights[(i_leaf, v)] = x
    # T_weights[(v, k_leaf)] = T_weights[(k_leaf, v)] = dist_matrix[i_leaf, k_leaf] - x

    # Add leaf n back to T by creating limb (v, n)
    # of length limb
    T_neighbours[v].append(n-1)
    T_neighbours[n-1] = [v]
    T_weights[(n-1, v)] = T_weights[(v, n-1)] = limb
    return (T_neighbours, T_weights)


# f = open('ba7c.txt', 'r')
# n = int(f.readline())
n = 4

f = ["0   13  21  22",
"13  0   12  13",
"21  12  0   13",
"22  13  13  0"]
dist_matrix = np.zeros((n, n))
for i, line in enumerate(f):
    values = [int(m) for m in line.strip().split()]
    dist_matrix[i, :] = values

(T_neighbours, T_weights) = additive_phylogeny(dist_matrix, n)

# Print tree
for (v, w) in sorted(T_weights.keys()):
    print("{}->{}:{}".format(v, w, T_weights[(v, w)]))
