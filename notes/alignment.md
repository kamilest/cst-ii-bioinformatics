# Alignment lectures

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
* Initialise first row and column to 0.

### Complexity
* Time: $O(mn)$ 
* Space: $O(mn)$

for sequences of length $n$ and $m$.

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

## Banded dynamic programming
Assume strings are already similar and have few insertions/deletions, updating scorings within some distance from the diagonal ($|i-j| < k(N)$)

### Complexity

* Time $O(N \cdot k(N)) \ll O(N^2)$

## Alignment in linear memory

Based on the following subproblem:

*Find a middle edge in the alignment gap in linear space.*

**Input**: Two strings and scoring matrix.

**Output**: Middle edge in the alignment graph of these strings

### Method (Hirschberg algorithm)
Computing alignment with linear memory

* Use two columns at a time and throw away others
* Detect middle edge based on where the score is highest in the middle column (forward and reverse)
* Recursively compute other midpoints, throwing the rest of the scoring matrix away

### Complexity
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

### Complexity

* $O(n^2)$ terms
* each taking $O(n)$ in case of bifurcation
* Time: $O(n^3)$
* Space: $O(n^2)$