# Middle Edge in Linear Space Problem

http://rosalind.info/problems/ba5k/

*Find a middle edge in the alignment graph in linear space.*

**Given**: Two amino acid strings.

**Return**: A middle edge in the alignment graph of these strings, where the optimal path is defined by the [BLOSUM62](http://rosalind.info/glossary/blosum62/) scoring matrix and a linear indel penalty equal to 5. Return the middle edge in the form `(i, j) (k, l)`, where `(i, j)` connects to `(k, l)`.

Sample Dataset
```
PLEASANTLY
MEASNLY
```
Sample Output
```
(4, 3) (5, 4)
```