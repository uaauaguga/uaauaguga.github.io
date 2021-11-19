---
layout: post
title:  "Probablistic Model for Sequence Analysis"
date:   2021-11-17 16:15:22 +0800
usemathjax: true
categories: jekyll update
---

## HMM

- 看了不知道多少遍总是记不清楚，详细记录一下，希望别再忘了:)

- Hidden state at time $$i$$: $$\pi_{i}$$
- Emitted symbol at time $$i$$: $$x_{i}$$
- Transition probability from state $$s$$ to $$t$$
 
  $$a_{st} = P(x_{i}=t \mid x_{i-1}=s)$$

- Probability of emitting symbol $$b$$ at state $$k$$

$$e_{k}(b)=P(x_{i}=b \mid \pi_{i}=k)$$

### Viterbi algorithm (decoding)
- Find the most probable path of hidden states, given the model and observation

### Forward algorithm (inference)
- Determine the likelihood of an observed series, marginalized for all possible hidden states

### Impose sparsity to the transition matrix

### Model training
- Estimation HMM parameter from observations
- Direct optimize for the likelihood by gradient descent
- Baum–Welch algorithm: forward-backward algorithm based, expectation maximization style training
- Viterbi training: viterbi algorithm based model training


## Profile HMM

- A spacial HMM primarily used for protein homolog modeling

- The hidden states in a profile HMM for local alignment, with N key positions
  - start state
  - insertion state 0, modeling left flanking background sequence
  - match state at key position 1, insertion state at key position 1, deletion state at key position 1
  - ...
  - match state at key position N, insertion state at key position N (modeling right flanking background sequence), deletion state at key position N
  - end state


- Possible transitions
  - begin $$\Rightarrow$$ insertion state 0
  - begin $$\Rightarrow$$ deletion state at key position 1 (allow condition when first residue is missing)
  - begin $$\Rightarrow$$  match state at key position 1
  
  - insertion state 0 $$\Rightarrow$$ insertion state 0 
      - self looping for left flanking background sequence modeling
      - The looping probability on the flanking states should be close to 1, since they must account for long stretches of sequence
  - insertion state 0 $$\Rightarrow$$ deletion state at key position 1
  - insertion state 0 $$\Rightarrow$$ match state at key position 1

  - match state at key position i-1 $$\Rightarrow$$ match state at key position i
  - insertion state at key position i-1 $$\Rightarrow$$ match state at key position i
  - deletion state at key position i-1 $$\Rightarrow$$ match state at key position i

  - match state at key position i-1 $$\Rightarrow$$ deletion state at key position i
  - insertion state at key position i-1 $$\Rightarrow$$ deletion state at key position i
  - deletion state at key position i-1 $$\Rightarrow$$ deletion state at key position i

  - match state at key position i $$\Rightarrow$$ insertion state at key position i
  - insertion state at key position i $$\Rightarrow$$ insertion state at key position i
  - deletion state at key position i $$\Rightarrow$$ insertion state at key position i 
    - "We have added transitions between insert and delete states, as they did, although these are usually very improbable. **Leaving them out has negligible effect on scoring a match, but can create problems when building the model from MSA**"
    - In most cases, profile HMM is constructed from MSA. Consider Figure 5.7 in [Biological Sequence Analysis](https://www.cambridge.org/core/books/biological-sequence-analysis/921BB7B78B745198829EF96BC7E0F29D). We can see given a column assignment for MSA, the $$D_i \Rightarrow I_i$$ transition could occur

{:refdef: style="text-align: center;"}
![]({{site.baseurl}}/images/bsa-figure-5.7.png){: height="70%" width="70%"}
{:refdef}

  - match state at key position i $$\Rightarrow$$ end
  - insertion state at key position i $$\Rightarrow$$ end
    - self looping for right flanking background sequence modeling, should be close to 1
  - deletion state at key position i $$\Rightarrow$$ end

- HMM Profile searching
  - Viterbi algorithm and forward algorithm can both be used for HMM profiling searching
  - Viterbi: assign a mostly likely hidden state for each nucleotide, and report the corresponding likelihood
  - Viterbi variable for matching, insertion and deletion states at key position j and sequence position i

$$
\begin{equation}
V_j^M(i) = \log{\frac{e_{M_j}(x_i)}{q_{x_i}}} + 
max \left\{
\begin{aligned}
  V_{j-1}^M(i-1) + \log a_{M_{j-1}M_j} \\
  V_{j-1}^I(i-1) + \log a_{I_{j-1}M_j}   \\
  V_{j-1}^D(i-1) + \log a_{D_{j-1}M_j}  \\
\end{aligned}
\right.
\end{equation}
$$

$$
\begin{equation}
V_j^I(i) = \log{\frac{e_{I_j}(x_i)}{q_{x_i}}} + 
max \left\{
\begin{aligned}
  V_{j}^M(i-1) + \log a_{M_{j}I_j} \\
  V_{j}^I(i-1) + \log a_{I_{j}I_j}   \\
  V_{j}^D(i-1) + \log a_{D_{j}I_j}  \\
\end{aligned}
\right.
\end{equation}
$$

$$
\begin{equation}
V_j^I(i) = 
max \left\{
\begin{aligned}
  V_{j-1}^M(i-1) + \log a_{M_{j-1}D_j} \\
  V_{j-1}^I(i-1) + \log a_{I_{j-1}D_j}   \\
  V_{j-1}^D(i-1) + \log a_{D_{j-1}D_j}  \\
\end{aligned}
\right.
\end{equation}
$$

  - forward: assign a likelihood that the input sequence is generated from a HMM profile

- See <https://notebook.community/jmschrei/pomegranate/tutorials/B_Model_Tutorial_3_Hidden_Markov_Models>


### Model construction

- To construct profile HMM from MSA, we have to make a decision. Each column in the MSA should either be assigned to a key position (marked column), or an insertion (unmarked column).
- **For columns marked as key positions, residues are assigned to match state, gaps are assigned to deletion state**
- **For columns marked as insertions, residues are assigned to insertions, gaps are ignored**
- There are $$2^L$$ ways to mark a MSA with $$L$$ columns.
- Three ways for marking MSA columns
  - Manual
  - Heuristics: rule based assignment. For example, assigning all columns will more than a certain fraction of gap characters to insert states
  - MAP (maximum *a posteriori*)



