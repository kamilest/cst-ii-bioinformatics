# Alignment with Affine Gap Penalties Problem
http://rosalind.info/problems/ba5j/

*Construct a highest-scoring global alignment of two strings (with affine gap penalties).*

**Given**: Two amino acid strings $v$ and $w$, and a scoring matrix.

**Return**: The maximum alignment score between $v$ and $w$, followed by an alignment of $v$ and $w$ achieving this maximum score. 

Use the [BLOSUM62](http://rosalind.info/glossary/blosum62/) scoring matrix, a gap opening penalty of 11, and a gap extension penalty of 1.

Assume each $v$ and $w$ of length at most 100.

Sample Dataset
```
PRTEINS
PRTWPSEIN
```
Sample Output
```
8
PRT---EINS
PRTWPSEIN-
```

## Background
A **gap** is a contiguous sequence of spaces in a row of an alignment. One way to score gaps more appropriately is to define an **affine penalty** for a gap of length $k$ as $\sigma + \varepsilon(k − 1)$, where $\sigma$ is the **gap opening penalty**, assessed to the first symbol in the gap, and $\varepsilon$ is the **gap extension penalty**, assessed to each additional symbol in the gap. We typically select $\varepsilon$ to be smaller than $\sigma$ so that the affine penalty for a gap of length $k$ is smaller than the penalty for $k$ independent single-nucleotide indels ($\sigma k$).