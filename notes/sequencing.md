# Genome sequencing
* modern sequencing machines cannot read entire genome one nucleotide at a time
* generate short *reads*
* reconstructing the genome from *overlapping reads*

## Genome sequencing problem

*Reconstruct a genome from reads.*

**Input:** collection of reads

**Output:** string corresponding to the genome

## String reconstruction
*Reconstruct a string from its $k$-mer composition*

**Input:** collection of $k$-mers

**Output:** genome such that $\text{Composition}_k(\text{Genome})$ is equal to the collection of $k$-mers.

Leads to the reconstruction of genome as a *path* between the different reads so that every $k$-mer is visited (perhaps every *type* of $k$-mer, accounting for multiplicity as several copies of the same DNA are used at a time)

* however repeats would be more frequent than the general multiplicity

### Hamiltonian Path
* Interpret $k$-mers as nodes, add an edge between nodes $v \leadsto w$ if $\text{Suffix}(v) = \text{Prefix}(w)$
* Creates a graph of all possible connections between $k$-mers.
* Find the path such that *all $k$-mers are visited once*.

#### Complexity
NP-Complete.

### Eulerian Path

* Interpret $k$-mers as edges, nodes on every side of the edge correspond to prefixes and suffixes of neighbouring edges.
* Merge the *nodes* which have the *same label* (prefix/suffix) on them for the de Bruijn graph.
* To reconstruct the genome find the *Eulerian path*, visiting each *edge* exactly once.

### De Bruijn graphs

* represent every $k$-mer as an edge between the prefix and suffix
* all nodes with identical labels merged together
* find the *Eulerian cycle* which can be easily converted to a path.

#### Complexity (graph construction)

* for each $k$-mer, 1 edge and up to 2 nodes
* $O(1)$
* assume hash map encodes nodes+edges
* assume $(k-1)$-mers fit in $O(1)$ machine words, hashing is $O(1)$
* querying/adding key is $O(1)$
* so $O(1)$ per $k$-mer, $O(n)$ for $n$ $k$-mers.

### Euler's theorem
*A graph is Eulerian if it contains an Eulerian cycle. The Eulerian cycle exists if the graph is balanced, i.e. for every node the indegree is equal to the outdegree.*

* Every Eulerian graph is balanced
* Every balanced and strongly connected graph is Eulerian
* A node is *semi-balanced* if indegree differs from outdegree by 1
  * can be the case with genomes in perfect coverage case where the first node does not match the last node
  * add a pseudo-edge to make the graph balanced, run Eulerian cycle

### Method (Hierholzer)
1. Start at a random node and keep walking until a cycle is constructed.
2. If there are nodes with unused edges, start at that node, use the unexplored edges and extend the previously constructed cycle with a second loop. The new cycle will start at that node which will be visited twice to traverse both of the loops.
3. Keep going until all edges are used.

Alternatively,
1. Form a *Cycle* by randomly walking in balanced graph
2. While *Cycle* is not Eulerian
   1. select the node *newStart* in *Cycle* with still unexplored outgoing edges
   2. form a *Cycle'* by traversing *Cycle* from *newStart* by randomly walking
   3. assign *Cycle'* to current *Cycle*


#### Complexity
Time: $O(|E|)$ for $E$ edges if the implementation is efficient.

## Read-pair assembly
A *paired $k$-mer* is a pair of $k$-mers at a *fixed* distance $d$ apart in the genome.
* (in real world), the $k$-mers are sampled from different anti-symmetric strands
* (in real world), the distance between $k$-mers is not exact and can include errors

### String reconstruction from read-pairs

*Reconstruct a string from its paired $k$-mers.*

**Input:** collection of paired $k$-mers.

**Output:** string such that its $\text{PairedComposition}$ is equal to the collection of paired $k$-mers.

### Method (paired de Bruijn)

* Construct de Bruijn graph as always but edges have paired labels and nodes are paired
* Still glue nodes but now such that *paired* labels are the same. 
* Decreases the number of possibilities for a path—simple unpaired de Bruijn usually gives multiple possible paths so many potential genome assemblies. 
* Paired de Bruijn narrows the assembly choice down and accounts for the placement of various *repeats* in the genome.

## Constraints
Assumptions include:

* *perfect coverage* of genome by reads
  * breaking reads into shorter $k$-mers which are more likely to cover the entire genome
  * long reads for perfect coverage should contain every possible $k$-mer in the composition which is unlikely
* *error-free* reads
  * in reality *bubbles* of de Bruijn graph are possible
  * bubble detection (could be due to errors or due to actual mutations—learn to distinguish those cases)
* *known $k$-mer multiplicities*
  * actually unknown but can estimate relative to the multiplicities of everything else (or average multiplicity)
* *distances* between reads within read-pairs are *exact*
