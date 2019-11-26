import numpy as np

def centers_to_clusters(centers, data):
    pass

def clusters_to_centers(clusters):
    pass

def lloyd_algorithm(data, k):
    # Generate arbitrary centers.
    data_ = data.copy()
    [xs, ys] = [list(t) for t in zip(*data_)]

    xcenter = np.random.uniform(min(xs), max(xs), size=k)
    ycenter = np.random.uniform(min(ys), max(ys), size=k)
    centers = list(zip(xcenter, ycenter))

    
    

    pass

# Parsing input
f = open('ba8c.txt', 'r')
[k, m] = [int(m) for m in f.readline().strip().split()]

data = []
for line in f:
  [x, y] = [float(m) for m in line.strip().split()]
  data.append((x, y))

