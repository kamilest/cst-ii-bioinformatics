---
layout: default
title: Genome pattern matching
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

* Time: $O(\vert \text{Patterns}\vert  + \vert \text{Genome}\vert  \times \vert \text{Pattern}\vert)$
* Space: $O(\vert \text{Patterns}\vert)$

## Suffix tree compression

Attempts to avoid storing the pattern trie in memory (reads (pattenrs) from human genome could be upwards of 1TB)

### Method
1. Construct a trie out of the *suffixes* of Genome
  * Append dollar signs to the end of the Genome text
2. Compress the non-branching paths into a single path to avoid storing unnecessary pointers to leaves and the Genome itself—
   * store the starting position of the string
   * length of the substring until the next branching point

### Complexity (suffix tries)
* Space: $O(\vert \text{Genome}\vert $ (but this is still impractical because of huge constant factors which make the genome not fit into RAM).


## Burrows-Wheeler transform

### Method

1. add an end-of-string character $
2. form a $N\times N$ matrix by cyclically rotating the text to form the rows
3. lexicographically sort the rotated rows
4. the last column is the Burrows-Wheeler transform, return it run-length coded
5. additionally return index to the original unrotated string in the BWT

### Complexity (BWT)

* $O(N^2 \log N)$ for a string of length $N$. 


## Inverse Burrows-Wheeler transform

### Method
1. get the last column in run-length coding
2. compute the first column by sorting the last column (making use of last-to-first property)
3. starting from the index of the original unrotated string

```
start from L[I] (set i=I):

while not back at L[I] again:
  add L[i] to output stack
  the char c = L[i] is j-th occurrence of c in array L
  let k be index of j-th occurrence of c in array F (last-to-first)
  next entry is then L[k] => set i=k

pop the characters from the output stack
```

### Complexity (inverse-BWT)
* $O(N \log N)$ for string of length $N$

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


## Multiple approximate pattern matching problem

*Find all approximate occurrences of a collection of patterns in a text.*

**Input:** a string `Text`, a collection of shorter strings `Patterns`, integer *d*.

**Output:** all starting positions in `Text` where a string from `Patterns` appears as a substring with at most *d* mismatches.

### Method 

1. *Seeding*: if at most $d$ mismatches are permitted and the `Text` is split into $d+1$ parts, then by pigeonhole principle at least one part will be matching completely
2. Sample the slices and see if there is one perfect slice match (seed detection).
3. Extend matching seed in both directions to verify further if `Pattern` occurs with at most $d$ mismatches (seed extension).

## BLAST: comparing a sequence against a database

Basic local alignment search tool.

* looks for similarities between proteins
* tolerant to extreme mutations of particular amino acids (looks for similarities in overall structure)
* very fast and works at large scale (compares against all proteins in the database)