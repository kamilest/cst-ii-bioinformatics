---
layout: default
title: Global alignment problem
parent: Alignment
nav_order: 2
---

## Global alignment problem
*Find the highest-scoring alignment between two strings using a scoring matrix.*

**Given**: Two amino acid strings and a scoring matrix.

**Return**: The maximum alignment score of these strings followed by an alignment achieving this maximum score. 

### Method (Needleman-Wunsch algorithm)

[Implementation](https://github.com/kamilest/cst-ii-bioinformatics/blob/master/problems/BA5E/ba5e.py)

1. Initialise first row and column to the `(cells to source) * penalty`
2. Dynamic programming: highest score between indel and match/mismatch
```
score = scoring_matrix[(v[i-1], w[j-1])]
s[i, j] = max(s[i-1, j] - indel_penalty,
              s[i, j-1] - indel_penalty
              s[i-1, j-1] + score)
```
3. Backtrack as appropriate
  * to the left (if insertion in column string)
  * up (if insertion in row string)
  * diagonal (match/mismatch) 
4. Going back through the coordinates from sink to source returns the original path.

### Complexity
* Space: $O(mn)$
* Time: $O(mn)$
* Backtrace: $O(m+n)$
