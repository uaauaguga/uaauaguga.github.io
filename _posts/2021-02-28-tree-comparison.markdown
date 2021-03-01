---
layout: post
title:  "Tree Comparisons"
date:   2021-02-28 14:54:30 +0800
usemathjax: true
categories: jekyll update
---

## Tree Comparisons


## Tree edit model
- Finding a series of edit operations that
transforms one input tree into the other with minimum overall
cost, defined as the accumulated cost of the basic edit operations.
- Basic edit operations
  - replacements
  - insertions
  - deletion
- Tree alignment
  - give two input trees, create a super tree to minimize the cost
- Tree edit
  - give two input trees, find a common subtree to minimize the cost

## Implementations
- Zhang-Shasha algorithm 
  -  <http://www.grantjenks.com/wiki/_media/ideas:simple_fast_algorithms_for_the_editing_distance_between_tree_and_related_problems.pdf>
  - <https://zhang-shasha.readthedocs.io/en/latest/>
  - A tutorial <https://arxiv.org/pdf/1805.06869.pdf>

## Example applications 

### RNA structural comparison
- RNAdistance

  - Calculate distances between RNA secondary structures
  -  <https://www.tbi.univie.ac.at/RNA/RNAdistance>

- RNAforester

  -  compare RNA secondary structures via forest alignment
  - <https://www.tbi.univie.ac.at/RNA/RNAforester>
- Gardenia
  - Multiple alignment of RNA secondary structure
  - <https://bioinfo.lifl.fr/RNA/gardenia/gardenia.php>


### Tracking clonal evolution of tumor cells
- Detecting evolutionary patterns of cancers using consensus trees <https://doi.org/10.1093/bioinformatics/btaa801>





