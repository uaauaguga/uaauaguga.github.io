#!/bin/bash
for rfam in $(cat data/Rfam/info/rfam-id-st50-lt10.txt);do
bsub -q TEST-1U "make -f makefiles/prepare-data.makefile rfamid=${rfam} source=seed-alignments > log/prepare-data/${rfam}-seed-alignments.log 2>&1"
done
