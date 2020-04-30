---
layout: default
title: Affine gap variant
parent: Alignment
nav_order: 4
---

## Affine gap variant

*Construct a highest-scoring global alignment of two strings (with affine gap penalties).*

**Input:** Two strings, a scoring matrix Score, and numbers $\sigma$ and $\varepsilon$.

**Output:** A highest scoring global alignment between these strings, as defined by Score and by the gap opening and extension penalties $\sigma$ and $\varepsilon$.

### Method
[Implementation](https://github.com/kamilest/cst-ii-bioinformatics/blob/master/problems/BA5J/BA5J.md)

* Fixed penalty $\sigma$ might be a problem
* Affine gap penalty for gap length $k$: $\sigma + \varepsilon(k-1)$
  * $\sigma$ gap opening penalty
  * $\varepsilon$ gap extension penalty
  * $\sigma > \varepsilon$ opening gap is more expensive than extending
* Use three levels of backtracking matrices
  * Lower level for insertions (down arrows)
  * Middle for matches/mismatches (diagonal arrows)
  * Upper level for deletions (right arrows)
* Lower and upper level penalties are $\varepsilon$ as they account for extensions
* Switching from middle level to lower/upper costs $\sigma$ gap opening penalty
* Going back from gap extension levels back to the same node in middle level does not cost anything
  * middle level updates to max of previous middle + score/same node in lower/same node in upper
  * lower updates to max of lower + $\varepsilon$ and middle level + $\sigma$
  * upper updates to max of upper + $\varepsilon$ and middle level + $\sigma$

Alternatively,

* $F$ scoring matrix would compute matches/gap opening/gap close (score from $G$ plus score of match)
* $G$ scoring matrix would compute gap extensions
* to remember overall optimal alignment need to know the best alignment if gap is still open or if gap is closed

### Complexity

* Time $O(nm)$ for $n > m$
* Space $O(nm)$