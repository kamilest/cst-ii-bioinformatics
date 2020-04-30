---
title: Alignment
author: Kamile Stankeviciute
category: Alignment
layout: post
mathjax: true
---

* Align genome and protein sequences
* Detect differences at the single base to block of bases
* RNA: aligning molecule with itself
* *Dynamic programming algorithms*
* BA Chapter 5


# Longest common subsequence

https://github.com/kamilest/ii-bioinformatics/tree/master/BA5C

**Given**: Two strings.

**Return**: A longest common subsequence of these strings.

**Complexity**: $O(nm)$ for sequences of length $n$ and $m$ (DP algorithm)

* Allow only insertions and deletions but no mismatches
* Scoring matrix initialise to 0
* +1 for matches, 0 for indels
* Aim for the maximal number of symbols for strings `v` and `w` such that `v_t = w_t`
* Edit distance over Hamming distance–DP method of solving through edit graphs


# Needleman-Wunsch

https://github.com/kamilest/ii-bioinformatics/tree/master/BA5E

## Global alignment problem
*Find the highest-scoring alignment between two strings using a scoring matrix.*

**Given**: Two amino acid strings and a scoring matrix.

**Return**: The maximum alignment score of these strings followed by an alignment achieving this maximum score. 

**Complexity**

* Space: $O(mn)$
* Time: $O(mn)$
* Backtrace: $O(m+n)$

## Method

* Initialise first row and column to the `cells to source * penalty`
* Dynamic programming: highest score between
  * indel
  * match or mismatch
* Backtrack as appropriate: left, right, diagonal step. Going back through the coordinates from sing to source returns the original path.

## Overlap detection variant
* Can have unlimited number of gaps at the beginning and the end of the sequence 
* Initialise first row and column to 0

## Smith-Waterman

https://github.com/kamilest/ii-bioinformatics/tree/master/BA5F

### Local alignment problem
*Find the highest-scoring local alignment between two strings.*

**Given**: Two amino acid strings and a scoring matrix.

**Return**: The maximum score of a local alignment of the strings, followed by a local alignment of these strings achieving the maximum score.

**Complexity**:

* Space: $O(mn)$
* Time: $O(mn)$
* Backtrace: $O(m+n)$

### Method

* Initialise first row and column to 0
* Dynamic programming: highest score between
  * 0
  * indel
  * match or mismatch

## Affine gap

* Fixed penalty $\sigma$ might be a problem
* Affine gap penalty for gap lenth $k$: $\sigma + \varepsilon(k-1)$
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
  * middle level updates to max of previous middle + score/same node in lower/same node in uppter
  * lower and uppter update to max of lower/upper + $\varepsilon$/middle level + $\sigma$

### Simulating affine gaps by long edges

* $F$ scoring matrix would compute matches/gap opening/gap close (score from $G$ plus score of match)
* $G$ scoring matrix would compute gap extensions
* to remember overall optimal alignment need to know the best alignment if gap is still open or if gap is closed

**Complexity** (alignment with gaps):

* Time $O(N^2M)$ for $N > M$
* Space $O(NM)$

## Banded dynamic programming
* assume strings are already similar and have few insertions/deletions
* update scorings within some distance from the diagonal ($|i-j| < k(N)$)

**Complexity**:

* Time $O(N \cdot k(N)) \ll O(N^2)$

## Hirschberg algorithm
Computing alignment with linear memory

* Use two columns at a time and throw away others
* Detect middle edge based on where the score is highest in the middle column (forward and reverse)
* Recursively compute other midpoints, throwing the rest of the scoring matrix away

### Middle edge in linear space problem.
*Find a middle edge in the alignment gap in linear space.*

**Input**: Two strings and scoring matrix.

**Output**: Middle edge in the alignment graph of these strings

**Complexity**:

* Time $O(nm)$
* Space $O(n)$

## Nussinov RNA folding algorithm
* DP matrix tracking the maximum number of bonds if the inner string between the two indices in the matrix folds optimally 
* if $i < j$:
  1. $x_i$ paired with $x_j$:
     1. $F(i, j) = s(x_i, x_j) + F(i+1, j-1)$
     2. score of folding of $x_i$ and $x_j$ as substring from $i+1$ to $j-1$ folds optimally
  2. $x_i$ is not paired with $x_j$:
     1. $F(i, j) = \mathrm{max}\{k: i \leq k < j\}\ F(i, k) + F(k+1, j)$
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

**Complexity**:

* $O(n^2)$ terms
* each taking $O(n)$ in case of bifurcation
* Time: $O(n^3)$
* Space: $O(n^2)$