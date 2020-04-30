---
layout: default
title: Local alignment problem
parent: Alignment
nav_order: 3
---

## Local alignment problem
When the global alignment is not as important and we care about the locally aligned sequences. 

*Find the highest-scoring local alignment between two strings.*

**Given**: Two amino acid strings and a scoring matrix.

**Return**: The maximum score of a local alignment of the strings, followed by a local alignment of these strings achieving the maximum score.

### Method (Smith-Waterman algorithm)

1. Initialise first row and column of the alignment matrix to 0 (do not add penalty for initial indels).
2. Dynamic programming 
```
score = scoring_matrix[(v[i-1], w[j-1])]
s[i, j] = np.amax([0, \
            s[i-1, j] - indel_penalty, \
            s[i, j-1] - indel_penalty, \
            s[i-1, j-1] + score])
```
3. Changes in backtracking: now backtrack to source if the maximum score is 0, and start backtracking from any cell in alignment matrix.

The main change from the Needleman-Wunsch is that the *score cannot be negative*.

### Complexity
* Space: $O(mn)$
* Time: $O(mn)$
* Backtrace: $O(m+n)$
