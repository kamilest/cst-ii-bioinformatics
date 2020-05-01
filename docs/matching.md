---
layout: default
title: Hidden Markov models
nav_order: 5
---

# Genomics pattern matching

* mapping low-divergent sequences against a large reference genome
* using a reference genome to help reconstruct the new sample genome faster by recognising the parts matching the reference rather than assembling from scratch
* *read mapping*: determining where each read has high similarity to the reference genome
* *reference genome*: template based on multiple individual genomes
* *fitting alignment*: aligning each pattern to the best substring in genome
  * runtime $O(\vert \text{Pattern}\vert  \times \vert \text{Genome}\vert )$ for a single pattern
  * $O(\vert \text{Patterns}\vert  \times \vert \text{Genome}\vert )$ for a collection of patterns

## Single pattern matching problem
*Find where reads match the reference genome exactly.*

**Input:** string *Pattern* and *Genome*

**Output:** positions in *Genome* where *Pattern* appears as a substring


## Trie data structure
* each trie has a single root node with indegree 0
* each edge of the trie is labelled with a letter of alphabet
* edges leaving a given node have distinct labels
* every string in *Patterns* is spelled out by concatenating letters along some path from root downward
* every path from the root to a *leaf*, or node with outdegree 0, spells a string in Patterns.


### Trie construction problem

*Construct a trie from a collection of patterns*

**Input:** collection of strings *Patterns*.

**Output:** the trie corresponding to *Patterns*.

* sliding the trie down the genome
* walk down the trie to see if the substring starting at that genome position matches any pattern in the trie.

### Method (trie construction problem)

```python
def construct_trie(patterns):
    trie = []
    trie.append({})

    for pattern in patterns:
        current = 0
        for i in range(len(pattern)):
            if pattern[i] in trie[current]:
                current = trie[current][pattern[i]]
            else:
                trie[current][pattern[i]] = len(trie)
                current = len(trie)
                trie.append({})

    return trie
```

### Complexity
* Time: $O(\vert \text{Patterns} \vert * n)$ where $n$ is the length of a pattern.


## Multiple pattern matching problem
**Input:** collection of *Patterns* and *Genome*

**Output:** all positions in *Genome* where a string from *Patterns* appears as a substring

### Complexity (brute force)
  * $O(\vert \text{Pattern}\vert  \times \vert \text{Genome}\vert )$ for a single pattern
  * $O(\vert \text{Patterns}\vert  \times \vert \text{Genome}\vert )$ for a collection of patterns (applying the Single pattern matching algorithm for every pattern)

### Method (prefix trie matching)

1. Construct a trie from the *prefixes* of Patterns (therefore prefix trie).
2. Run each suffix of Genome through the trie.

### Complexity (prefix trie matching)

* Time: $O(|\text{Patterns}| + |\text{Genome}| \times |\text{Pattern}|$
* Space: $O(|\text{Patterns}|$

## Suffix tree compression

Attempts to avoid storing the pattern trie in memory (reads (pattenrs) from human genome could be upwards of 1TB)

### Method
1. Construct a trie out of the *suffixes* of Genome
  * Append dollar signs to the end of the Genome text
2. Compress the non-branching paths into a single path to avoid storing unnecessary pointers to leaves and the Genome itself—
   * store the starting position of the string
   * length of the substring until the next branching point

### Complexity (suffix tries)
* Space: $O(|\text{Genome}|$ (but this is still impractical because of huge constant factors which make the genome not fit into RAM).


## Burrows-Wheeler transform

### Method

1. form a $N\times N$ matrix by cyclically rotating the text to form the rows
2. lexicographically sort the rotated rows
3. the last row is the Burrows-Wheeler transform

### Last-to-first property
The $i$-th occurrence of a character in last column corresponds to the $i$-th occurrence of that character in the first column

## Pattern matching with Burrows-Wheeler transform
* Each row of the BWT of the text corresponds to a suffix of text, so the indices of the characters in the first column correspond to a characters in the suffix array
* However, undesirable to store the entire BWT matrix other than the first and last columns because of the high memory requirement. How to reconstruct unknown second, third etc. symbol when they are not stored

### Method
1. Move backward through the pattern
2. Match the $n$-th symbol of the Pattern to the FirstColumn (easy by binary search)
3. Look up the LastColumn of the matched suffixes to get $(n-1)$-th symbol
4. Look up which of those match to $(n-1)$-th symbol of the Patten
5. Repeat