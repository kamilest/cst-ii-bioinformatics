# Limb Length Problem

http://rosalind.info/problems/ba7b/

*Find the limb length for a leaf in a tree.*

**Given**: An integer $n$, followed by an integer $j$ between $0$ and $n - 1$, followed by a space-separated additive distance matrix $D$ (whose elements are integers).

**Return**: The limb length of the leaf in `Tree(`$D$`)` corresponding to row $j$ of this distance matrix (use 0-based indexing).

Sample Dataset
```
4
1
0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0
```
Sample Output
```
2
```