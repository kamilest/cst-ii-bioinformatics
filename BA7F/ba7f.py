import numpy as np
import re

class ParsimonyNode:
  def __init__(self, label, tag=0, characters=''):
    self.label = label
    self.tag = tag
    self.characters = characters

    self.min_score = 0

    # Stores scores for all the characters rather than running the 
    # algorithm for each position in the string
    self.scores = []

    for c in characters:
      self.scores.append({'A': np.inf, 'C': np.inf, 'T': np.inf, 'G': np.inf})
      self.scores[-1][c] = 0

    self.left = None
    self.right = None
  
  def ripe(self):
    if self.left is None and self.right is None:
      return True
    return self.left.tag == 1 and self.right.tag == 1

def small_parsimony(node):
  small_parsimony_bottom_up(node)
  small_parsimony_top_down(node)

def small_parsimony_bottom_up(node):
  def delta(i, k):
    return 0 if i == k else 1
  
  if node is None or node.left is None or node.right is None:
    return

  # Postorder traversal
  small_parsimony_bottom_up(node.left)
  small_parsimony_bottom_up(node.right)
  node.tag = 1

  for i in range(len(node.left.scores)):
    node.scores.append({})

  for ix, s in enumerate(node.scores):
    for k in 'ACTG':
      si = [node.left.scores[ix][i] + delta(i, k) for i in 'ACTG']
      sj = [node.right.scores[ix][j] + delta(j, k) for j in 'ACTG']
      
      s[k] = min(si) + min(sj)
  

def small_parsimony_top_down(node):
  if node is None or node.left is None or node.right is None:
    return
  
  min_scores = []
  for s in node.scores:
    sc = [s[k] for k in 'ACTG']
    min_scores.append(min(sc))

    sc_ix = np.argmin(sc)
    # node.characters += 'ACTG'[sc_ix]
  
  node.min_score = np.sum(min_scores)

  small_parsimony_top_down(node.left)
  small_parsimony_top_down(node.right)

def print_small_parsimony(T, T_adj):
  def hamming_distance(characters_i, characters_j):
    sum = 0
    for i, j in zip(characters_i, characters_j):
      sum += 0 if i == j else 1
    return sum 

  print(T[max(T)].min_score)
  for (v, w) in sorted(T_adj):
    print(v, w)
    print("{}->{}:{}".format(T[v].characters, T[w].characters, hamming_distance(T[v].characters, T[w].characters)))

T = {}
T_adj = []
label = 0

f = open('ba7f.txt', 'r')
n = int(f.readline())

for line in f:
  [v, w] = re.split('->', line.strip())
  v = int(v)
  if v not in T:
      T[v] = ParsimonyNode(v)

  if w.isdigit():
    u = int(w)
    if u not in T:
      T[u] = ParsimonyNode(u)
  else:
    u = label
    T[u] = ParsimonyNode(u, tag=1, characters=w)
    label+=1
  
  if T[v].left is None:
    T[v].left = T[u]
  else:
    T[v].right = T[u]
  
  T_adj.append((u, v))
  T_adj.append((v, u))

f.close()

# Assuming maximum key is the root
root = T[max(T)]
small_parsimony(root)
print_small_parsimony(T, T_adj)