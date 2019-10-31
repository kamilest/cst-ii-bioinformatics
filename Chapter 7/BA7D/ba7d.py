import numpy as np

def upgma(dist_matrix, n):
    def cluster_distance(dist_matrix, c1, c2):
        distance = 0
        for i in c1:
            for j in c2:
                distance += dist_matrix[(c1, c2)]
        return distance / (len(c1) * len(c2))
    pass

f = open('ba7d.txt', 'r')
n = int(f.readline())

dist_matrix = np.zeros((n, n))
for i, line in enumerate(f):
    values = [int(m) for m in line.strip().split()]
    dist_matrix[i, :] = values

(T_neighbours, T_weights) = upgma(dist_matrix, n)

# Print tree
for (v, w) in sorted(T_weights.keys()):
    print("{}->{}:{}".format(v, w, int(T_weights[(v, w)])))