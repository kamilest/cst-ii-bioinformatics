# Trie Construction Problem
http://rosalind.info/problems/ba9a/

Construct a trie on a collection of patterns.

**Given**: A collection of strings *Patterns*.

**Return**: The adjacency list corresponding to `Trie(Patterns)`, in the following format. 

If `Trie(Patterns)` has `n` nodes, first label the root with `1` and then label the remaining nodes with the integers `2` through `n` in any order you like. Each edge of the adjacency list of `Trie(Patterns)` will be encoded by a triple: the first two members of the triple must be the integers labeling the initial and terminal nodes of the edge, respectively; the third member of the triple must be the symbol labeling the edge.

Sample Dataset
```
ATAGA
ATC
GAT
```

Sample Output
```
0->1:A
1->2:T
2->3:A
3->4:G
4->5:A
2->6:C
0->7:G
7->8:A
8->9:T
```

## Multiple Pattern Matching Problem

Find all occurrences of a collection of patterns in a text.

* **Input**: A string *Text* and a collection *Patterns* containing (shorter) strings.
* **Output**: All starting positions in *Text* where a string from *Patterns* appears as a substring.

