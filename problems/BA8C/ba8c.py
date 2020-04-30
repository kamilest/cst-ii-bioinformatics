import numpy as np

def distance(p, q):
    return np.linalg.norm(np.array(p)-np.array(q))

def centroid(ps):
    length = len(ps)
    sum_x = np.sum([p[0] for p in ps])
    sum_y = np.sum([p[1] for p in ps])
    return [sum_x/length, sum_y/length]

def centers_to_clusters(centers, data):
    clusters = [[] for _ in range(len(centers))]
    for p in data:
        distances = [distance(p, c) for c in centers]
        clusters[np.argmin(distances)].append(p)
    return clusters
            
def clusters_to_centers(clusters):
    return [centroid(c) for c in clusters]

def lloyd_algorithm(data, k):
    # Generate arbitrary centers.
    data_ = data.copy()
    [xs, ys] = [list(t) for t in zip(*data_)]

    xcenter = np.random.uniform(min(xs), max(xs), size=k)
    ycenter = np.random.uniform(min(ys), max(ys), size=k)
    
    new_centers = np.array([list(t) for t in zip(xcenter, ycenter)])
    centers = np.zeros(new_centers.shape)

    while(not np.allclose(centers, new_centers)):
        centers = new_centers
        clusters = centers_to_clusters(centers, data)
        new_centers = clusters_to_centers(clusters)
    
    return new_centers

# Parsing input
f = open('ba8c.txt', 'r')
[k, m] = [int(m) for m in f.readline().strip().split()]

data = []
for line in f:
  p = [float(m) for m in line.strip().split()]
  data.append(p)

kmeans = lloyd_algorithm(data, k)
for p in kmeans:
    print('{0:.3f} {1:.3f}'.format(p[0], p[1]))

