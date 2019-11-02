# Edit Distance Problem

http://rosalind.info/problems/ba5g/

*Find the edit distance between two strings.*

**Given**: Two amino acid strings.

**Return**: The edit distance between these strings.

Sample dataset
```
PLEASANTLY
MEANLY
```

Sample output
```
5
```

## Wagner-Fischer algorithm
* Dynamic programming
* Store the edit distances of the prefixes
* Just like **Needleman-Wunsch** but:
  * the goal is to find the *minimum* length path
  * indel penalty +1 (1 edit)
  * mismatch penalty +1 (1 substitution edit)
  * match score 0
  