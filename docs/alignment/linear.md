---
layout: default
title: Alignment in linear memory
parent: Alignment
nav_order: 5
---

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