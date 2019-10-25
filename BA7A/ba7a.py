import numpy as np
import re
from queue import Queue

def find_leaves(neighbours):
    leaves = []
    for v in neighbours:
        if len(neighbours[v]) == 1:
            leaves.append(v)
    return leaves

def distances_from_leaf(neighbours, weights, leaf, leaves): 
    distances = {leaf: 0}
    next_nodes = Queue()
    next_nodes.put(leaf)

    while not next_nodes.empty():
        v = next_nodes.get()
        for w in neighbours[v]:
            if w not in distances.keys():
                next_nodes.put(w)
                distances[w] = distances[v] + weights[(v, w)]
    
    return [distances[leaf] for leaf in leaves]

def distance_between_leaves(neighbours, weigths):
    leaves = find_leaves(neighbours)
    dist_matrix = np.zeros((len(leaves), len(leaves)))
    for i, leaf in enumerate(leaves):
        dist_matrix[:, i] = distances_from_leaf(neighbours, weights, leaf, leaves)
    
    return dist_matrix.astype(int)

f = open("rosalind_ba7a.txt", "r")
n = int(f.readline())

neighbours = {}
weights = {}

for line in f:
    [v, w, weight] = [int(s) for s in re.split('->|:', line)]
    if v not in neighbours:
        neighbours[v] = []
    if w not in neighbours:
        neighbours[w] = []
    neighbours[v].append(w)
    # Adjacency list contains weights in both directions
    weights[(v, w)] = weight
f.close()

dist_matrix = distance_between_leaves(neighbours, weights)
for row in dist_matrix:
    print(' '.join(map(str, row)))