# Global Alignment Problem

[http://rosalind.info/problems/ba5e/]()

*Find the highest-scoring alignment between two strings using a scoring matrix.*

**Given**: Two amino acid strings and a scoring matrix.

**Return**: The maximum alignment score of these strings followed by an alignment achieving this maximum score. 

Use the [BLOSUM62](http://rosalind.info/glossary/blosum62/) scoring matrix and indel penalty σ = 5. (If multiple alignments achieving the maximum score exist, you may return any one.)

Sample dataset
```
PLEASANTLY
MEANLY
```

Sample output
```
8
PLEASANTLY
-MEA--N-LY
```

## Needleman-Wunsch algorithm
* Initialise first row and column to the `cells to source * penalty`
* Dynamic programming: highest score between
  * indel
  * match or mismatch
* Backtrack as appropriate: left, right, diagonal step. Going back through the coordinates from sing to source returns the original path.