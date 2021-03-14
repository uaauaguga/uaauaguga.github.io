#!/bin/bash
for path in $(ls benchmark/alignments/seed-alignments);do
for L in 200 400 800;do
#L=100
bsub -q TEST-1U "make -f  makefiles/prepare-data.makefile rfamid=${path%.*} flanking=${L} > log/prepare-data/${path%.*}-fl-${L}.log"
done
done
