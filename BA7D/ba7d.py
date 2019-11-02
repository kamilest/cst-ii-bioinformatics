import numpy as np

def upgma(dist_matrix, n):
    class Cluster:
        def __init__(self, label, son=None, daughter=None, size=1):
            self.label = label
            self.son = son
            self.daughter = daughter
            self.size = size
    
    def closest_clusters(cluster_distances, clusters):
        ci, cj = None, None
        min_distance = np.inf
        for i in clusters:
            for j in clusters:
                if (i != j):
                    distance = cluster_distances[i][j]
                    if distance < min_distance:
                        ci = i
                        cj = j
                        min_distance = distance

        return ci, cj, min_distance

    def remove_cluster_distances(cluster_distances, ci_label, cj_label):
        del cluster_distances[ci_label]
        del cluster_distances[cj_label]

        for c in cluster_distances:
            del cluster_distances[c][ci_label]
            del cluster_distances[c][cj_label]

    def compute_cnew_distances(clusters, cluster_distances, ci, cj):
        cnew_distances = {}
        for cm in clusters:
            cnew_distances[cm] = \
                (cluster_distances[ci.label][clusters[cm].label] * ci.size + \
                 cluster_distances[cj.label][clusters[cm].label] * cj.size) / \
                    (ci.size + cj.size)
        
        return cnew_distances
    
    def add_cnew_distances(clusters, cluster_distances, cnew_label, cm_label, distance):
        cluster_distances[cnew_label][cm_label] = distance
        cluster_distances[cm_label][cnew_label] = distance
    
    def add_weights(T_weights, T_ages, cluster):
        if cluster.son is None or cluster.daughter is None:
            return T_weights
        
        else:
            label = cluster.label
            son_label = cluster.son.label
            daughter_label = cluster.daughter.label

            label_age = T_ages[label]
            son_age = T_ages[son_label]
            daughter_age = T_ages[daughter_label]

            T_weights[(label, son_label)] = abs(label_age - son_age)
            T_weights[(son_label, label)] = abs(label_age - son_age)

            T_weights[(label, daughter_label)] = abs(label_age - daughter_age)
            T_weights[(daughter_label, label)] = abs(label_age - daughter_age)

            T_weights = add_weights(T_weights, T_ages, cluster.son)
            T_weights = add_weights(T_weights, T_ages, cluster.daughter)
        
        return T_weights

    clusters = {}
    cluster_distances = {}
    cnew = len(dist_matrix) - 1
    T_ages = {}
    for i in range(len(dist_matrix)):
        clusters[i] = Cluster(i)
        T_ages[i] = 0
    
    for i in range(len(dist_matrix)):
        cluster_distances[i] = {}
        for j in range(len(dist_matrix)):
            cluster_distances[i][j] = dist_matrix[i, j]

    while len(clusters) > 1:
        ci, cj, distance = closest_clusters(cluster_distances, clusters)
        ci_cluster = clusters[ci]
        cj_cluster = clusters[cj]

        # Merge clusters
        cnew += 1
        cnew_cluster = Cluster(cnew, ci_cluster, cj_cluster, \
            ci_cluster.size + cj_cluster.size)

        T_ages[cnew] = distance / 2.

        # Remove Ci and Cj from clusters
        del clusters[ci]
        del clusters[cj]

        # Computing D(Cnew, C) for each C in clusters
        cnew_distances = \
            compute_cnew_distances(clusters, cluster_distances, ci_cluster, cj_cluster)

        remove_cluster_distances(cluster_distances, ci, cj)
        cluster_distances[cnew] = {}
        for cm in cnew_distances:
            add_cnew_distances(clusters, cluster_distances, cnew, cm, cnew_distances[cm])

        clusters[cnew] = cnew_cluster

    # Binary tree traversal
    for cluster in clusters.values():
        T_weights = add_weights({}, T_ages, cluster)

    return T_weights

f = open('rosalind_ba7d.txt', 'r')
n = int(f.readline())

dist_matrix = np.zeros((n, n))
for i, line in enumerate(f):
    values = [int(m) for m in line.strip().split()]
    dist_matrix[i, :] = values

T_weights = upgma(dist_matrix, n)

# Print tree
for (v, w) in sorted(T_weights.keys()):
    print("{}->{}:{}".format(v, w, "{0:.3f}".format(T_weights[(v, w)])))