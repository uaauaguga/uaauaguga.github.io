#!/bin/bash
for rfam in $(cat data/Rfam/info/rfam-id-st50-lt10.txt);do
echo $rfam
wget -O data/Rfam/full-alignments/fasta/${rfam}.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/${rfam}.fa.gz
done
