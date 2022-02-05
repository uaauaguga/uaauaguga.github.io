---
layout: post
title:  "Notes for deep learning in biological sequence analysis"
date:   2021-06-25 13:51:23 +0800
usemathjax: true
categories: jekyll update
---

- A collection of studies that applies neural network to learn representations for biological sequence

### Related Resources

- 2017, Cell Systems, [Enhancing Evolutionary Couplings with DeepConvolutional Neural Networks](https://linkinghub.elsevier.com/retrieve/pii/S2405-4712(17)30542-2)

- 2019, [Pre-Training of Deep Bidirectional Protein Sequence Representations with Structural Information](https://arxiv.org/abs/1912.05625)

- 2019, Nature Method, [Unified rational protein engineering with sequence-based deep representation learning](https://www.nature.com/articles/s41592-019-0598-1)
  - MLSTM, minimize next amino-acid prediction cross entropy loss, use fixed size hidden states as feature representation
  - One hidden state activity for each amino acid 
  - Average the activity cross all AAs in the full length sequence to get representation for the protein 
  - Another choice is use the activity of the last hidden state
  - doc2vec: <https://github.com/fhalab/embeddings_reproduction>
  - <https://github.com/churchlab/UniRep>
  

- 2018, NeurIPS, [Neural Edit Operations for Biological Sequences](https://proceedings.neurips.cc/paper/2018/file/d0921d442ee91b896ad95059d13df618-Paper.pdf)
  - Replace argmax with softmax in sequence alignment, to make the sequence alignment loss differentiable 
  - Related works:
    - Differentiable DTW loss for time series: <https://arxiv.org/pdf/1703.01541.pdf> 
    - Sequence alignment kernel: 2004, *Bioinformatics*, [Protein homology detection using string alignment kernels](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bth141)

- 2021, ICLR, [Bertology Meets Biology Interpreting Attention Protein Language Moldes](https://arxiv.org/abs/2006.15222)

- 2021, Bioinformatics, [DNABERT: pre-trained Bidirectional Encoder Representations from Transformers model for DNA-language in genome](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab083/6128680)
  - <https://github.com/jerryji1993/DNABERT>

- 2021, PNAS, [Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences](https://www.pnas.org/content/118/15/e2016239118)
  - <https://github.com/facebookresearch/esm>

- 2021, Current Potocols, [Learned Embeddings from Deep Learning to Visualize and Predict Protein Sets](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.113)
  - <https://github.com/sacdallago/bio_embeddings>