---
layout: post
title:  "Notes on rRNA annotation"
date:   2022-01-25  15:37:11 +0800
usemathjax: true
categories: jekyll update
---


- <https://en.wikipedia.org/wiki/Ribosomal_RNA>
- <https://en.wikipedia.org/wiki/Ribosome_biogenesis>


## SSU rRNA and LSU rRNA

- [SSU rRNA](https://en.wikipedia.org/wiki/SSU_rRNA) (small subunit rRNA) 
  - 16S rRNA in bacteria, archea and plastid
  - 18S rRNA in eukaryotes
  - 12S rRNA in mitochondria

- [LSU rRNA](https://en.wikipedia.org/wiki/LSU_rRNA) (large subunit rRNA) 
  - 23S rRNA in bacteria, archea and plastid
  - 5.8S & 28S rRNA in eukaryotes. The 5.8S rRNA is homologous to the 5' end of non-eukaryotic LSU rRNA. In eukaryotes, the insertion of ITS2 breaks LSU rRNA into 5.8S and 28S rRNAs.
  - 16S in Mitochondria
  - 5S rRNA

- Biogenesis
  - In bacteria, 16S rRNA, 23S rRNA, and 5S rRNA are typically organized as a co-transcribed operon.
  - In eukaryotes, the 28S, 5.8S, and 18S rRNAs are encoded by a single transcription unit (45S) separated by 2 internally transcribed spacers, and transcribed by RNA polymerase I. 5S rRNA is transcribed independetly by RNA polymerase III. 
  

## Where to get rRNA data

- [barrnap](https://github.com/tseemann/barrnap) is widely used tool for rRNA annotation, [SortMeRNA](https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases) is a tool to sort out rRNA from NGS data. Here list data resources used by these two tools 
  - [barrnap](https://github.com/tseemann/barrnap)
  - [SortMeRNA](https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases)


- We can see two major sources of rRNA information come from:
  - The [silva](https://www.arb-silva.de/) database. Provide sequence alignments, focus on metagenomics (16S rRNA is the most wildely used marker gene in metagenomics).
  - Also see 
    - [greengenes](https://greengenes.secondgenome.com/) for 16S rRNA
    - [RDP](https://rdp.cme.msu.edu/) Bacterial and Archaeal 16S rRNA sequences, and Fungal 28S rRNA
  - The [Rfam](https://rfam.xfam.org/) database. Provide sequence alignments, reference structure and covariance model.
    - 5S and 5.8S rRNA presents in ribosome large subunit, but not included in SILVA LSU. [barrnap](https://github.com/tseemann/barrnap) and [SortMeRNA](https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases) use Rfam for 5S and 5.8S rRNA annotation.
    - 5S rRNA: RF00001
    - 5.8S rRNA: RF00002



## rRNA sequence in reference genome assembly

- rRNA is generally recognized as repeat sequence in eukaryotic genomes, and abscent from some genome assemblies.

- Human genome
  - 5S, 5.8S, 18S, 28S rRNAs are all presents in [hg38 assembly](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz).
  - [gencode v38 annotation](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz) contains two types of genomic rRNA (5S, 5.8S)and two types of mitochondrial rRNA (16S, 12S), but not genomic 16S and 28S rRNA.
  - refseq provides annotations for all 4 types of rRNAs, see [this link](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz).
  - See discussion in here: <http://seqanswers.com/forums/showthread.php?t=41868> 

- rRNAs have considerable variation among populations
  - See 2018, *Science Advances*, [Variant ribosomal RNA alleles are conserved and exhibit tissue-specific expression](https://www.science.org/doi/10.1126/sciadv.aao0665?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
    - Mouse rRNA prototypes
      - NR_003278.3
      - NR_003279.1
      - NR_003280.2
      - NR_030686.1
    - Mouse rDNA prototypes
      - BK000964.3
      - chromosome 8 [123538334,123539354]
    - Human rRNA prototypes
      - NR_003285.2
      - NR_003286.2
      - NR_003287.2
      - NR_023379.1
    - Human rDNA prototypes
      - U13369.1
      - X12811.1

GenBank IDs of the rDNA prototypes used were U13369.1 and X12811.1
(human) and  (mouse).

