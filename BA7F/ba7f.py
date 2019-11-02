import numpy as np
import re

class ParsimonyNode:
  def __init__(self, label, tag=0, characters=''):
    self.label = label
    self.tag = tag
    self.characters = characters

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



def small_parsimony(root):
  def delta(i, k):
    return 0 if i == k else 1
  
  if root is None or root.left is None or root.right is None:
    return

  small_parsimony(root.left)
  small_parsimony(root.right)
  root.tag = 1

  for i in range(len(root.left.scores)):
    root.scores.append({})

  for ix, s in enumerate(root.scores):
    for k in 'ACTG':
      si = [root.left.scores[ix][i] + delta(i, k) for i in 'ACTG']
      sj = [root.right.scores[ix][j] + delta(j, k) for j in 'ACTG']
      print(si, sj)

      s[k] = min(si) + min(sj)

T = {}
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

f.close()

# Assuming maximum key is the root
small_parsimony(T[max(T)])
scores = []
for s in T[max(T)].scores:
  scores.append(min([s[k] for k in 'ACTG']))

print(np.sum(scores))