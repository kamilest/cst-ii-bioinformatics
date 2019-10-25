# import numpy as np
from queue import Queue
import re

def longest_path_dag(source, sink, neighbours, weights):
    s = {} # node -> (max_weight, prev_node)
    s[source] = (0, None)
    # visited = set()

    next_nodes = Queue()
    next_nodes.put(source)

    while not next_nodes.empty():
        v = next_nodes.get()
        # visited.add(v)

        for w in neighbours[v]:
            if w != sink:
                next_nodes.put(w)
            if w in s:
                weight = max(s[w][0], s[v][0] + weights[(v, w)])
                if weight == s[v][0] + weights[(v, w)]:
                    s[w] = (weight, v)
            else:
                s[w] = (s[v][0] + weights[(v, w)], v)
    
    weight = s[sink][0]

    # backtrack for the path
    path = [sink]
    current = sink
    while current is not source:
        current = s[current][1]
        path.append(current)
    
    return (weight, path)


f = open("rosalind_ba5d.txt", "r")
source = int(f.readline())
sink = int(f.readline())

neighbours = {}
weights = {}

for line in f:
    [v, w, weight] = [int(s) for s in re.split('->|:', line)]
    if v not in neighbours:
        neighbours[v] = []
    if w not in neighbours:
        neighbours[w] = []
    neighbours[v].append(w)
    weights[(v, w)] = weight

weight, path = longest_path_dag(source, sink, neighbours, weights)
print(weight)
print('->'.join(map(str, path[::-1])))