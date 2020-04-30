---
layout: default
title: Hidden Markov models
nav_order: 5
---

# Hidden Markov Models
- [Hidden Markov Models](#hidden-markov-models)
  - [Three problems of HMMs](#three-problems-of-hmms)
  - [Probability of a hidden path](#probability-of-a-hidden-path)
  - [Probability of an outcome given a hidden path](#probability-of-an-outcome-given-a-hidden-path)
  - [Decoding problem](#decoding-problem)
    - [Method (Viterbi algorithm)](#method-viterbi-algorithm)
    - [Complexity](#complexity)
  - [Evaluation problem (outcome likelihood)](#evaluation-problem-outcome-likelihood)
    - [Method (Forward algorithm)](#method-forward-algorithm)
    - [Complexity](#complexity-1)
  - [Learning problem](#learning-problem)
    - [Method (Backward algorithm)](#method-backward-algorithm)
  - [Baum-Welch learning](#baum-welch-learning)

**Definition.** HMM:

* alphabet $\Sigma$ of emitted symbols
* set of states $S$
* transition matrix $T \in \mathbb{R}^{\vert S\vert  \times \vert S\vert }$ representing the probability of going from one state to another, $\forall j. \sum_{i}T_{ij} = 1$
* emission matrix $E \in \mathbb{R}^{\vert S\vert  \times \vert \Sigma\vert }$, representing the probability of emitting a given symbol from a given state, $\forall j. \sum_{i} E_{ij} = 1$


**Definition.** Hidden path (or *parse*) $\pi = \pi_1 \dots \pi_n$: sequence of states HMM passes through.

**Definition.** Sequence $x = x_1\dots x_n$ of emitted symbols.

## Three problems of HMMs
1. **Evaluation.** Given a HMM $M$ and a sequence $x$, find probability of the *observable sequence* $\mathrm{Pr}(x\vert M)$ over all possible paths (i.e. $\sum_\pi \mathrm{Pr}(x, \pi\vert M)$)—solved by *Forward algorithm*
2. **Decoding.** Given a HMM $M$ and a sequence $x$, find the *hidden state sequence* $\pi$ that maximises $\mathrm{Pr}(x, \pi\vert M)$ (i.e. $\mathrm{argmax}_\pi \mathrm{Pr}(x, \pi \vert M)$)—solved by *Viterbi algorithm*
3. **Learning.** Given a HMM $M$ with unspecified transition/emission probabilities and a sequence $x$, find parameters $\theta = (E, T)$ that maximise $\mathrm{Pr}(x\vert \theta)$—solved by *Backward algorithm*/*Baum-Welch learning*

where model $M$ is defined by architecture – alphabet $\Sigma$, states $S$ (or $Q$) – and parameters $\theta = (E, T)$ ($T$ can also be denoted as $A$ or $a_{ij}$, $E$ can be denoted as $e_i[\cdot]$).

## Probability of a hidden path

*Compute the probability of an HMM's hidden path.*

**Input:** hidden path $\pi$ in an HMM $(\Sigma, S, T, E)$

**Output:** probability of the path, $\mathrm{Pr}(\pi)$

$$\mathrm{Pr}(\pi) = \prod_{i=1}^n T_{\pi_{i-1}\pi_i}$$

## Probability of an outcome given a hidden path

*Compute the probability that an HMM will emit a string given its hidden path.*

**Input:** string (or *sequence*) $x = x_1\dots x_n$ emitted by an HMM $(\Sigma, S, T, E)$ and a hidden path $\pi = \pi_1\dots \pi_n$.

**Ouput:** conditional probability $\mathrm{Pr}(x\vert \pi)$ that $x$ will be emitted given that the HMM follows hidden path $\pi$.

\[\mathrm{Pr}(x\vert \pi) = \prod_{i=1}^n \mathrm{Pr}(x_i \vert  \pi_i) = \prod_{i=1}^n E_{\pi_i x_i}\]


We also have

\[\mathrm{Pr}(x, \pi) = \mathrm{Pr}(x\vert \pi)\mathrm{Pr}(\pi) = \prod_{i=1}^n \mathrm{Pr}(x_\vert \pi_i)\mathrm{Pr}(\pi_i) = \prod_{i=1}^n T_{\pi_{i-1} \pi_i} E_{\pi_i x_i}\]


## Decoding problem

*Find the most optimal hidden path in an HMM given a string of its emitted symbols.*

**Input:** string $x = x_1\dots x_n$ emitted by an HMM $(\Sigma, S, T, E)$

**Output:** path $\pi$ that maximises probability $\mathrm{Pr}(x, \pi)$ over all possible paths through this HMM.


### Method (Viterbi algorithm)
 dynamic programming: the solution from the source to sink corresponds to the maximum of the solutions from the first node to sink multiplied by probabilities of transitioning to that node.

1. define $s_{k,i}$ as the weight of the optimal path from *source* to node $(k,i)$
2. *dynamic programming*: first $i-1$ transitions from source to $(k,i)$ must be from an optimal path from source to $(l, i-1)$ for some state $l$.

```python
s[k,i] = max(s[l, i-1] * weight[(l, i-1), (k, i)])
       = max(s[l, i-1] * T[path[i-1], path[i]] * E[path[i], symbol[i]])
```
3. source is connected to all first nodes with equal weight so the first step depends only on emissions

```
s[0,0] = 1
s[j,0] = 0 for all j > 0

iterate:
  s[k,i] = max(s[l, i-1] * T[path[i-1], path[i]] * E[path[i], symbol[i]])
  ptr[k,i] = argmax_j(T[j, k] * s[j, i-1])

terminate:
  P(x, path*) = max(s[k, N])

traceback:
  path[N]* = argmax_j(s[j, N])
  path[i-1]* = ptr[path[i], i]
```

### Complexity

Runtime is linear in terms of the *number of edges in the Viterbi graph*, which is $O(\vert S\vert ^2 n)$ where $n$ is the number of emitted symbols (length of string).

* Time: $O(\vert S\vert ^2 n)$
* Space: $O(\vert S\vert n)$


## Evaluation problem (outcome likelihood)

*Find the probability that an HMM emits a given string.*

**Input:** string $x=x_1\dots x_n$ emitted by an HMM $(\Sigma, S, T, E)$

**Output:** probability $\mathrm{Pr}(x)$ that the HMM emits $x$.

### Method (Forward algorithm)
Slight change to the Viterbi algorithm, replacing maximisation with the sum.

```python
s[k,i] = sum(s[l, i-1] * weight[(l, i-1), (k, i)])
       = sum(s[l, i-1] * T[path[i-1], path[i]] * E[path[i], symbol[i]])
```
*useful technique:* rescale at each position by multiplying by a constant

### Complexity

* Time: $O(\vert S\vert ^2n)$
* Space: $O(\vert S\vert n)$

## Learning problem

*Compute the probability distribution of hidden states given the observed states.*

**Input:** an HMM $M$, sequence of observed states $x$.

**Output:** $\mathrm{Pr}(\pi_i =k \vert  x)$

### Method (Backward algorithm)

\[\mathrm{Pr}(\pi_i=k\vert x) = \mathrm{Pr}(x_1\dots x_i, \pi=k, x_{i+1}\dots x_N)\]

\[... = \mathrm{Pr}(x_1\dots x_i, \pi_i = k) \mathrm{Pr}(x_{i+1} \dots x_N \vert  x_1\dots x_i, \pi_i = k)\]

\[... = \mathrm{Pr}(x_1\dots x_i, \pi_i = k)\mathrm{Pr}(x_{i+1}\dots x_N\vert  \pi_i = k)\]

*useful technique:* rescale at each position by multiplying by a constant

## Baum-Welch learning

*Expectation maximisation (EM) algorithm for parameter estimation.*

Start with estimators for parameters:

\[T_{lk}' = \frac{T_{lk}}{\sum_j T_{lj}}\]

\[E_{kb}' = \frac{E_{kb}}{\sum_{c\in \Sigma} E_{kc}}\]