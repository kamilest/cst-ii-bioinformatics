# Overlap Alignment Problem

http://rosalind.info/problems/ba5i/

*Construct a highest-scoring overlap alignment between two strings.*

**Given**: Two protein strings $v$ and $w$, and a scoring matrix.

**Return**: The score of an optimal overlap alignment of $v$ and $w$, followed by an alignment of a suffix $v'$ of $v$ and a prefix $w'$ of $w$ achieving this maximum score. 

Use an alignment score in which matches count +1 and both the mismatch and indel penalties are 2. (If multiple overlap alignments achieving the maximum score exist, you may return any one.)

Assume each $v$ and $w$ of length at most 1000.

Sample Dataset
```
PAWHEAE
HEAGAWGHEE
```

Sample Output
```
1
HEAE
HEAG
```

## Background
When we assembled genomes, we discussed how to use overlapping reads to assemble a genome, a problem that was complicated by errors in reads. We would like to find overlaps between error-prone reads as well.

An **overlap alignment** of strings $v = v_1 \dots v_n$ and $w = w_1 \dots w_m$ is a global alignment of a suffix of $v$ with a prefix of $w$. An optimal overlap alignment of strings $v$ and $w$ maximizes the global alignment score between an $i$-suffix of $v$ and a $j$-prefix of $w$ (i.e., between $v_i \dots v_n$ and $w_1 \dots w_j$) among all $i$ and $j$.