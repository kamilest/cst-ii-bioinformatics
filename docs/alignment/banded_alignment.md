---
layout: default
title: Banded alignment variant
parent: Alignment
nav_order: 4
---

## Banded dynamic programming
Assume strings are already similar and have few insertions/deletions, updating scorings within some distance from the diagonal ($|i-j| < k(N)$)

### Complexity

* Time $O(N \cdot k(N)) \ll O(N^2)$