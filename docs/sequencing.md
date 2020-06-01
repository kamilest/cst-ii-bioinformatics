---
layout: default
title: Genome sequencing
nav_order: 4
---
- [Genome sequencing](#genome-sequencing)
  - [Genome sequencing problem](#genome-sequencing-problem)
  - [String reconstruction](#string-reconstruction)
    - [Hamiltonian path](#hamiltonian-path)
    - [Complexity](#complexity)
    - [Eulerian path](#eulerian-path)
    - [De Bruijn graphs](#de-bruijn-graphs)
    - [Complexity (De Bruijn graph construction)](#complexity-de-bruijn-graph-construction)
  - [Finding Eulerian cycle in De Bruijn graph](#finding-eulerian-cycle-in-de-bruijn-graph)
    - [Euler's theorem](#eulers-theorem)
    - [Method (Hierholzer)](#method-hierholzer)
    - [Complexity](#complexity-1)
  - [Read-pair assembly](#read-pair-assembly)
    - [String reconstruction from read-pairs](#string-reconstruction-from-read-pairs)
    - [Method (paired de Bruijn)](#method-paired-de-bruijn)
    - [Constraints](#constraints)

# Genome sequencing
* modern sequencing machines cannot read entire genome one nucleotide at a time
* generate short *reads*
* reconstructing the genome from *overlapping reads*

## Genome sequencing problem

*Reconstruct a genome from reads.*

**Input:** collection of reads

**Output:** string corresponding to the genome

## String reconstruction
*Reconstruct a string from its $k$-mer composition. A $k$-mer is an arbitrary substring of length $k$.*

**Input:** collection of $k$-mers

**Output:** genome such that $\text{Composition}_k(\text{Genome})$ is equal to the collection of $k$-mers.

Leads to the reconstruction of genome as a *path* between the different reads so that every $k$-mer is visited (perhaps every *type* of $k$-mer, accounting for multiplicity as several copies of the same DNA are used at a time)

* however repeats would be more frequent than the general multiplicity

### Hamiltonian path
* Interpret $k$-mers as nodes, add an edge between nodes $v \leadsto w$ if $\text{Suffix}(v) = \text{Prefix}(w)$
* Creates a graph of all possible connections between $k$-mers.
* Find the path such that *all $k$-mers are visited once*.

### Complexity
NP-Complete.

### Eulerian path

* Interpret $k$-mers as edges, nodes on every side of the edge correspond to prefixes and suffixes of neighbouring edges.
* Merge the *nodes* which have the *same label* (prefix/suffix) on them for the de Bruijn graph.
* To reconstruct the genome find the *Eulerian path*, visiting each *edge* exactly once.

### De Bruijn graphs

* represent every $k$-mer $s_i, \dots, s_{i + k - 1}$ as an edge between a suffix node ($s_{i + 1}, \dots, s_{i + k - 1}$) and a prefix node ($s_i, \dots, s_{i + k - 2}$), with edge label $s_i, \dots, s_{i + k - 1}$.
* all nodes with identical labels merged together
* find the *Eulerian cycle* which can be easily converted to a path.

### Complexity (De Bruijn graph construction)

* for each $k$-mer, 1 edge and up to 2 nodes
* $O(1)$
* assume hash map encodes nodes+edges
* assume $(k-1)$-mers fit in $O(1)$ machine words, hashing is $O(1)$
* querying/adding key is $O(1)$
* so $O(1)$ per $k$-mer, $O(n)$ for $n$ $k$-mers.

## Finding Eulerian cycle in De Bruijn graph

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

Alternatively, ()
1. Form a *Cycle* by randomly walking in balanced graph
2. While *Cycle* is not Eulerian
   1. select the node *newStart* in *Cycle* with still unexplored outgoing edges
   2. form a *Cycle'* by traversing *Cycle* from *newStart* by randomly walking
   3. assign *Cycle'* to current *Cycle*


### Complexity
Time: $O(|E|)$ for $E$ edges if the implementation is efficient.

### Problems
With the above methods, we can arrive at multiple Eulerian paths, and some of them may not correspond to the k-mer composition of the genome. We may attempt to arrive at a unique reconstruction of the genome by increasing the read length. However, this is not always possible due to limitations of the sequencing machine. A more elegant method is to solve the sequencing problem on a de Bruijn graph constructed from paired k-mers (Below). 

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

**Problem:** Not every Eulerian path in the paired de Bruijn graph constructed from a (k, d)-mer composition spells out a
solution of the String Reconstruction from Read-Pairs Problem. To solve this, generate all Eulerian paths and output the path that spells out a string whose paired (k, d)-mer composition is equal to the initial set of (k, d)-mers.

## Constraints
Assumptions include:

* *perfect coverage* of genome by reads 
  * in reality, only a small fraction of k-mers are captured from the genome, thereby violating the key assumption of de Bruijn graphs.
  * breaking reads into shorter $k$-mers which are more likely to cover the entire genome.
  * long reads for perfect coverage should contain every possible $k$-mer in the composition which is unlikely
  * However, Even after read breaking, most assemblies still have gaps in k-mer coverage, causing the de Bruijn graph to have missing edges, and so the search for an Eulerian path fails. In this case, biologists often settle on assembling contigs (long, contiguous segments of the genome) rather than entire chromosomes. For example, a typical bacterial sequencing project may result in about a hundred contigs, ranging in length from a few thousand to a few hundred thousand nucleotides. For most genomes, the order of these contigs along the genome remains unknown. Needless to say, biologists would prefer to have the entire genomic sequence, but the cost of ordering the contigs into a final assembly and closing the gaps using more expensive experimental methods is often prohibitive. We can derive contigs from the de Bruijn graph, by picking out non-branching paths. In practice, biologists have no choice but to break genomes into contigs, even in the case of perfect coverage, since repeats prevent them from being able to infer a unique Eulerian path.
* *error-free* reads
  * in reality, we may perform an incorrect read. Since errors are independent, we may arrive at *bubbles* at the Bruijn graph. More formally, bubbles are short disjoint paths, shorter than some threshold length, that connect the same pair of nodes in the de Bruijn graph. The presence of bubbles clearly affects the correctness of our genome reconstruction.
  * A single rerror in read results in a bubble of length k in a de Bruijn graph.
  * Therefore, we may want to perform bubble detection (could be due to errors or due to actual mutations—learn to distinguish those cases). 
* *known $k$-mer multiplicities*
  * actually unknown but can estimate relative to the multiplicities of everything else (or average multiplicity). The multiplicity of a k-mer in a genome can often be estimated using its coverage. Indeed, k-mers that appear t times in a genome
are expected to have approximately t times higher coverage than k-mers that appear just once. Needless to say, coverage varies across the genome, and this condition is often violated. As a result, existing assemblers often assemble repetitive regions in genomes without knowing the exact number of times each k-mer from this region occurs in the
genome.
* *distances* between reads within read-pairs are *exact*

### Impact of k-mer length
The optimal $k$-mer size depends on the read length and the read depth and sequence complexity. 

* longer reads and/or higher read depth $\rightarrow$  use larger $k$-mers which could resolve complex areas of the graph 
* short reads are short and/or low read depth $\rightarrow$ use shorter $k$-mers
* many separate disconnected subgraphs in assembly graph (i.e. there are many small groups of contigs that have no connections to the rest of the graph) $\rightarrow$  $k$-mer size may be too large
* connected (i.e. all contigs are tied together in a single graph structure) but is very dense/tangled assembly graph $\rightarrow$ k-mer size may be too small
