---
layout: default
title: Nussinov RNA folding algorithm
parent: Alignment
nav_order: 6
---

## Nussinov RNA folding algorithm

* DP matrix tracking the maximum number of bonds if the inner string between the two indices in the matrix folds optimally 
* if $i < j$:
  1. $x_i$ paired with $x_j$:
     1. $F(i, j) = s(x_i, x_j) + F(i+1, j-1)$
     2. score of folding of $x_i$ and $x_j$ as substring from $i+1$ to $j-1$ folds optimally
  2. $x_i$ is not paired with $x_j$:
     1. $F(i, j) = \mathrm{max}_{k: i \leq k < j}\ F(i, k) + F(k+1, j)$
* Final $F(i, j)$ is the maximum of previous two cases
* Dot-bracket representation for *pseudoknot-free* structure; extended pseudoknot representation
* Secondary structure:
  * set of paired positions on interval [i, j]
  * optimal structure is built by extending optimal substructures
  * if optimal substructures known, the structure is formed in four ways
    * $i, j$ paired
    * $i$ unpaired
    * $j$ unpaired
    * combining two substructures
* Initialise DP matrix to $F(i, {i-1}) = 0$ and $F(i, i) = 0$
* Then $F(i, j)$ is the maximum of (a more straightforward illustration of the four cases)
  * $F(i+1, j)$
  * $F(i, j-1)$
  * $F(i+1, j-1) + s(i, j)$ where $s(i, j)$ is indicator of complementary base pairs
  * $\mathrm{max} \{k: i < k < j\}\ F(i, k) + F(k+1, j)$

### Pseudocode
```
for i in range(n-1):
  s[i,i] = 0
  s[i+1,i] = 0

s[n,n] = 0
for k in range(n):
  for i in range(n-k):
    j = i + k
    s[i,j] = np.amax([s[i+1,j], 
                      s[i,j-1], 
                      s[i+1,j-1] + score[v[i],v[j]],
                      max_i<k<j(s[i,k] + s[k+1,j])])
```

### Complexity

* $O(n^2)$ terms
* each taking $O(n)$ in case of bifurcation
* Time: $O(n^3)$
* Space: $O(n^2)$