---
layout: default
title: Longest common subsequence
parent: Alignment
nav_order: 1
---

## Longest common subsequence

**Given**: Two strings.

**Return**: A longest common subsequence of these strings.

### Method

Dynamic programming algorithm:

1. Initialise scoring matrix to 0
2. At each step allow only insertions or deletions (no mismatches)
3. Positive score for matches, no score for indels
4. Aim for the maximal number of symbols for strings `v` and `w` such that `v_t = w_t`

Edit distance over Hamming distance–DP method of solving through edit graphs.

### Overlap detection variant

* Can have unlimited number of gaps at the beginning and the end of the sequence.
* Initialialisation: set first row and column to 0. (arbitrary gap length from the beginning).
* Termination: Optimal Score: $\max(\max_{i}F[M, i]), \max_{j}(F[j, N]))$ (arbitrary gap length from the end).

### Complexity
* Time: $O(mn)$ 
* Space: $O(mn)$
