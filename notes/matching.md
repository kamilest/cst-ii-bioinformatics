# Genomics pattern matching

* mapping low-divergent sequences against a large reference genome
* using a reference genome to help reconstruct the new sample genome faster by recognising the parts matching the reference rather than assembling from scratch
* *read mapping*: determining where each read has high similarity to the reference genome
* *fitting alignment*: aligning each pattern to the best substring in genome
  * runtime $O(|\text{Pattern}| \times |\text{Genome}|)$ for a single pattern
  * $O(|\text{Patterns}| \times |\text{Genome}|)$ for a collection of patterns

## Single pattern matching problem
*Find where reads match the reference genome exactly.*

**Input:** string *Pattern* and *Genome*

**Output:** positions in *Genome* where *Pattern* appears as a substring

## Multiple pattern matching problem
**Input:** collection of *Patterns* and *Genome*

**Output:** all positions in *Genome* where a string from *Patterns* appears as a substring

#### Brute force approach complexity
  * $O(|\text{Pattern}| \times |\text{Genome}|)$ for a single pattern
  * $O(|\text{Patterns}| \times |\text{Genome}|)$ for a collection of patterns

## Using tries for pattern matching

* sliding the trie down the genome
* walk down the trie to see if the substring starting at that genome position matches any pattern in the trie.

**Complexity**



## Suffix tree compression

## Burrows-Wheeler transform

### Method

* form a $N\times N$ matrix by cyclically rotating the text to form the rows
* lexicographically sort the rotated rows
* the last row is the Burrows-Wheeler transform
* *last-to-first* property: the $i$-th occurrence of a character in last column corresponds to the $i$-th occurrence of that character in the first column