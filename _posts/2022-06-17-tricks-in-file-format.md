---
layout: post
title:  "Tricks on file formating"
date:   2022-06-17 16:21:48 +0800
usemathjax: true
categories: jekyll update
---

## Manipulate fasta file 

- <https://github.com/lh3/seqtk> is extremely useful

- Filter contigs shorter than 1000nt:

```bash
seqtk seq -L 1000 contigs.fa > contigs.gt1000.fa
# other functions include sampling sequence, etc.
```

- Convert one line fasta to multiline fasta

```bash
# reformat.sh provided in bbtools
reformat.sh fastawrap=70 in=input.fa out=output.fa
```


## Download refseq sequences

```python
#!/usr/bin/env python
from Bio import Entrez
import time
import argparse
import os
import logging

def main():
    parser = argparse.ArgumentParser(description='Download genebank sequence')
    parser.add_argument('--query', '-q',help="Genebank id for download",required=True)
    parser.add_argument('--fasta','-f',help="Path to write the downloaded data",required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
    logger = logging.getLogger('Retrieve sequence')

    Entrez.email = "xxx"
    fout = open(args.fasta,"w")

    try:
        logger.info("Start retriving {} from ncbi nucleotide in fasta format ...".format(args.query))
        handle = Entrez.efetch(db="nucleotide", id=args.query, rettype="fasta", retmode="text")
        content = handle.read()
        fout.write(content)
        logger.info("Done.")
    except:
        logger.info("Error retriving {}, skip ...".format(args.query))
    fout.close()

if __name__ == "__main__":
    main()

```