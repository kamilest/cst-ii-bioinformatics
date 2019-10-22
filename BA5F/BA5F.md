# Local Alignment Problem
http://rosalind.info/problems/ba5f/

*Find the highest-scoring local alignment between two strings.*

**Given**: Two amino acid strings and a scoring matrix.

**Return**: The maximum score of a local alignment of the strings, followed by a local alignment of these strings achieving the maximum score. 

Use the [PAM250](http://rosalind.info/glossary/pam250/) scoring matrix and indel penalty σ = 5. (If multiple local alignments achieving the maximum score exist, you may return any one.)

Sample dataset
```
MEANLY
PENALTY
```

Sample output
```
15
EANL-Y
ENALTY
```


## Smith-Waterman algorithm
Modified **Needleman-Wunsch** (global alignment algorithm)

* Ignore badly aligning regions
* *Best* local alignment: `s_opt = max_{i,j} s[i, j]`
* The local alignment algorithm introduces zero-weight edges connecting the source (0, 0) to every other node in the alignment graph, as well as zero-weight edges connecting every node to the sink node.
