# Hidden Markov Models

**Definition.** HMM:

* alphabet $\Sigma$ of emitted symbols
* set of states $S$
* transition matrix $T \in \mathbb{R}^{|S| \times |S|}$ representing the probability of going from one state to another, $\forall j. \sum_{i}T_{ij} = 1$
* emission matrix $E \in \mathbb{R}^{|S| \times |\Sigma|}$, representing the probability of emitting a given symbol from a given state, $\forall j. \sum_{i} E_{ij} = 1$


**Definition.** Hidden path (or *parse*) $\pi = \pi_1 \dots \pi_n$: sequence of states HMM passes through.

**Definition.** Sequence $x = x_1\dots x_n$ of emitted symbols.

## Three problems of HMMs
1. **Evaluation.** Given a HMM $M$ and a sequence $x$, find probability of the *observable sequence* $\mathrm{Pr}(x|M)$ over all possible paths (i.e. $\sum_\pi \mathrm{Pr}(x, \pi|M)$)—solved by *Forward algorithm*
2. **Decoding.** Given a HMM $M$ and a sequence $x$, find the *hidden state sequence* $\pi$ that maximises $\mathrm{Pr}(x, \pi|M)$ (i.e. $\mathrm{argmax}_\pi \mathrm{Pr}(x, \pi |M)$)—solved by *Viterbi algorithm*
3. **Learning.** Given a HMM $M$ with unspecified transition/emission probabilities and a sequence $x$, find parameters $\theta = (E, T)$ that maximise $\mathrm{Pr}(x|\theta)$—solved by *Baum-Welch learning*

where model $M$ is defined by architecture – alphabet $\Sigma$, states $S$ (or $Q$) – and parameters $\theta = (E, T)$ ($T$ can also be denoted as $A$ or $a_{ij}$, $E$ can be denoted as $e_i[\cdot]$).

### Probability of a hidden path

*Compute the probability of an HMM's hidden path.*

**Input:** hidden path $\pi$ in an HMM $(\Sigma, S, T, E)$

**Output:** probability of the path, $\mathrm{Pr}(\pi)$

$$\mathrm{Pr}(\pi) = \prod_{i=1}^n T_{\pi_{i-1}\pi_i}$$

### Probability of an outcome given a hidden path

*Compute the probability that an HMM will emit a string given its hidden path.*

**Input:** string (or *sequence*) $x = x_1\dots x_n$ emitted by an HMM $(\Sigma, S, T, E)$ and a hidden path $\pi = \pi_1\dots \pi_n$.

**Ouput:** conditional probability $\mathrm{Pr}(x|\pi)$ that $x$ will be emitted given that the HMM follows hidden path $\pi$.

$$\mathrm{Pr}(x|\pi) = \prod_{i=1}^n \mathrm{Pr}(x_i | \pi_i) = \prod_{i=1}^n E_{\pi_i x_i}$$


We also have

$$\mathrm{Pr}(x, \pi) = \mathrm{Pr}(x|\pi)\mathrm{Pr}(\pi) = \prod_{i=1}^n \mathrm{Pr}(x_|\pi_i)\mathrm{Pr}(\pi_i) = \prod_{i=1}^n T_{\pi_{i-1} \pi_i} E_{\pi_i x_i}$$

<!-- ### Performance
### Classifying proteins with profile HMMs

### Soft decoding problem -->

## Viterbi algorithm

### Decoding problem

*Find the most optimal hidden path in an HMM given a string of its emitted symbols.*

**Input:** string $x = x_1\dots x_n$ emitted by an HMM $(\Sigma, S, T, E)$

**Output:** path $\pi$ that maximises probability $\mathrm{Pr}(x, \pi)$ over all possible paths through this HMM.

**Method:** dynamic programming: the solution from the source to sink corresponds to the maximum of the solutions from the first node to sink multiplied by probabilities of transitioning to that node.

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

**Complexity:** runtime is linear in terms of the *number of edges in the Viterbi graph*, which is $O(|S|^2 n)$ where $n$ is the number of emitted symbols (length of string).

* Time: $O(|S|^2 n)$
* Space: $O(|S|n)$


## Forward-Backward algorithm

### Outcome likelihood problem

*Find the probability that an HMM emits a given string.*

**Input:** string $x=x_1\dots x_n$ emitted by an HMM $(\Sigma, S, T, E)$

**Output:** probability $\mathrm{Pr}(x)$ that the HMM emits $x$.

**Method.** *Forward algorithm*: Slight change to the Viterbi algorithm, replacing maximisation with the sum.

```python
s[k,i] = sum(s[l, i-1] * weight[(l, i-1), (k, i)])
       = sum(s[l, i-1] * T[path[i-1], path[i]] * E[path[i], symbol[i]])
```
*useful technique:* rescale at each position by multiplying by a constant

**Complexity**

* Time: $O(|S|^2n)$
* Space: $O(|S|n)$

## Backward algorithm

*Compute the probability distribution of hidden states given the observed states.*

**Input:** an HMM $M$, sequence of observed states $x$.

**Output:** $\mathrm{Pr}(\pi_i =k | x)$

**Method.**

$$\mathrm{Pr}(\pi_i=k|x) = \mathrm{Pr}(x_1\dots x_i, \pi=k, x_{i+1}\dots x_N)$$
$$ ... = \mathrm{Pr}(x_1\dots x_i, \pi_i = k) \mathrm{Pr}(x_{i+1} \dots x_N | x_1\dots x_i, \pi_i = k)$$
$$... = \mathrm{Pr}(x_1\dots x_i, \pi_i = k)\mathrm{Pr}(x_{i+1}\dots x_N| \pi_i = k)$$

*useful technique:* rescale at each position by multiplying by a constant

## Baum-Welch learning

*Expectation maximisation (EM) algorithm for parameter estimation.*

Start with estimators for parameters:

$$
T'_{lk} = \frac{T_{lk}}{\sum_j T_{lj}}
$$

$$
E'_{kb} = \frac{E_{kb}}{\sum_{c\in \Sigma} E_{kc}}
$$