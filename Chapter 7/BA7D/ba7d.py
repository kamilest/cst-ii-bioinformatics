import numpy as np

def upgma(dist_matrix, n):
    def cluster_distance(dist_matrix, c1, c2):
        distance = 0
        for i in c1:
            for j in c2:
                distance += dist_matrix[(i, j)]
        return distance / (len(c1) * len(c2))
    
    def closest_clusters(dist_matrix, clusters):
        ci, cj = None, None
        min_distance = np.inf
        for i in range(len(clusters)):
            for j in range(i+1, len(clusters)):
                distance = cluster_distance(dist_matrix, \
                    clusters[i], clusters[j])
                if distance < min_distance:
                    ci = i
                    cj = j
                    min_distance = distance

        return ci, cj, distance
    
    clusters = {}
    C_distances = {}
    T_neighbours = {}
    T_ages = {}
    for i in range(len(dist_matrix)):
        clusters[i] = set([i]) # (root_node, cluster)
        T_neighbours[i] = []
        T_ages[i] = 0
    
    for i in range(len(dist_matrix)):
        C_distances[i] = {}
        for j in range(len(dist_matrix)):
            C_distances[i][j] = dist_matrix[i, j]

    while len(clusters) > 1:
        ci, cj, distance = closest_clusters(C_distances, clusters)
        ci_node, ci_cluster = clusters[ci]
        cj_node, cj_cluster = clusters[cj]

        # Merge clusters
        cnew = np.amax([k for k, _ in T_neighbours.items()]) + 1
        cnew_cluster = ci_cluster.union(cj_cluster)
        T_ages[cnew] = distance / 2

        # Connect node Cnew to Ci and Cj
        T_neighbours[cnew] = [ci_node, cj_node]
        T_neighbours[ci_node].append[cnew]
        T_neighbours[cj_node].append[cnew]

        # Remove rows and columns of dist_matrix for Ci and Cj
        del C_distances[(ci, cj)]
        del C_distances[(cj, ci)]
        dist_matrix = np.delete(dist_matrix, [ci, cj], 0)
        dist_matrix = np.delete(dist_matrix, [ci, cj], 1)
        clusters.remove(clusters[ci])
        clusters.remove(clusters[cj])

        dist_to_cnew = [lambda c: cluster_distance(dist_matrix, c, cnew_cluster) for c in clusters]

        clusters.append((cnew, cnew_cluster))
        
        
    root = clusters[0][0]




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