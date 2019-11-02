import numpy as np

def add_edge(T_neighbours, T_weights, origin, destination, weight):
  if origin not in T_neighbours.keys():
    T_neighbours[origin] = []
  if destination not in T_neighbours.keys():
    T_neighbours[destination] = []
  
  T_neighbours[origin].append(destination)
  T_neighbours[destination].append(origin)
  T_weights[(origin, destination)] = weight
  T_weights[(destination, origin)] = weight

def limb_length(dist_matrix, leaf):
  leaf_row = dist_matrix[leaf, :]
  
  length = np.inf
  for i in range(len(leaf_row)):
    for j in range(i+1, len(leaf_row)):
      if i != leaf and j != leaf:
        length = min(length,\
            (dist_matrix[leaf, i] + dist_matrix[leaf, j] - dist_matrix[i, j]) \
              / 2.0)

  return length

def total_distance(dist_matrix, i):
  return np.sum(dist_matrix[i, :])

def neighbour_joining_matrix(dist_matrix):
  dist_nj = np.zeros(dist_matrix.shape)

  n = len(dist_matrix)
  for i in range(n):
    for j in range(n):
      if i != j:
        dist_nj[i, j] = (n - 2) * dist_matrix[i, j] \
          - total_distance(dist_matrix, i) \
          - total_distance(dist_matrix, j)
    
  return dist_nj

def neighbour_joining(dist_matrix, labels, n, total_n):
  T_weights = {}
  print(labels)

  if n == 2:
    T_weights[(m, m+1)] = dist_matrix[0, 1]
    T_weights[(m+1, m)] = dist_matrix[0, 1]
    return T_weights

  dist_nj = neighbour_joining_matrix(dist_matrix)

  # Find i and j such that dist_nj[i, j] is a minimum non-diagonal element of D'
  mask = np.ones(dist_nj.shape, dtype=bool)
  np.fill_diagonal(mask, np.inf)
  (i, j) = np.unravel_index(np.argmin(dist_nj[mask], axis=None), dist_nj.shape)
  i_label, j_label = labels[i], labels[j]

  delta = (total_distance(dist_matrix, i) - total_distance(dist_matrix, j)) / (n - 2)
  limb_i = (dist_matrix[i, j] + delta) / 2
  limb_j = (dist_matrix[i, j] - delta) / 2

  # Add a new row/column m to D so that Dk,m = Dm,k ... for any k
  m_col = (dist_matrix[:, i] + dist_matrix[:, j] - dist_matrix[i, j]) / 2

  dist_matrix = np.r_[dist_matrix, [m_col]]
  dist_matrix = np.c_[dist_matrix, np.zeros(dist_matrix.shape[0])]
  dist_matrix[:-1, -1] = m_col

  dist_matrix = np.delete(dist_matrix, [i, j], 0)
  dist_matrix = np.delete(dist_matrix, [i, j], 1)
  labels.remove(i_label)
  labels.remove(j_label)

  T_weights = neighbour_joining(dist_matrix, labels, n-1, m)
  labels.append(m)

  T_weights[(i_label, m)] = T_weights[(m, i_label)] = limb_i
  T_weights[(j_label, m)] = T_weights[(m, j_label)] = limb_j

  return T_weights

f = open('ba7e.txt', 'r')
n = int(f.readline())

dist_matrix = np.zeros((n, n))
for i, line in enumerate(f):
  values = [int(m) for m in line.strip().split()]
  dist_matrix[i, :] = values

# lengths to track node labels as columns in dist_matrix get removed
T_weights = neighbour_joining(dist_matrix, list(range(len(dist_matrix))), n, n)

# Print tree
for (v, w) in sorted(T_weights.keys()):
  print("{}->{}:{}".format(v, w, "{0:.3f}".format(T_weights[(v, w)])))