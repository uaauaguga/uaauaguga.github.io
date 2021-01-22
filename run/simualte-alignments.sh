#!/bin/bash
for rfam in $(cat data/Rfam/info/rfam-id-st50-lt10.txt);do
echo $rfam
bsub -q TEST-1U "cmemit -N 1000 -a -o benchmark/alignments/simulated-alignments/${rfam}.stk data/Rfam/cm-models/${rfam}.cm "
done

