## Longest Path in a DAG Problem
http://rosalind.info/problems/ba5d/

*Find a longest path between two nodes in an edge-weighted DAG.*

**Given**: An integer representing the source node of a graph, followed by an integer representing the sink node of the graph, followed by an edge-weighted graph. The graph is represented by a modified adjacency list in which the notation "0->1:7" indicates that an edge connects node 0 to node 1 with weight 7.

**Return**: The length of a longest path in the graph, followed by a longest path. (If multiple longest paths exist, you may return any one.)

Sample dataset
```
0
4
0->1:7
0->2:4
2->3:2
1->4:1
3->4:3
```

Sample output
```
9
0->2->3->4
```

Extra dataset output
```
62
0->14->29->44
```
