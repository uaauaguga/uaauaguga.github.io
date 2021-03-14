#!/bin/bash
#while IFS='\t' read -r rfamid length
#do
rfamid=RF00037
for alpha in 0.1 0.3 0.5 0.7 0.9;do
for size in 4 8 16;do
for init in sequence structure;do
#bsub -q TEST-1U "scripts/em-discovery-zoos-sequence-structure.py -i benchmark/datasets/seed-alignments/${rfamid}/ViennaRNA-400/pairing-probability-without-shape.txt --gamma 0.99 --alpha ${alpha} --seed-size ${size} --initialize ${init} --locations performance/tune/${rfamid}-${alpha}-${size}-${init}.bed -w 45 "
echo "${rfamid}-${alpha}-${size}-${init}"
scripts/evaluate-motif-location.py -i performance/tune/${rfamid}-${alpha}-${size}-${init}.bed -r benchmark/datasets/seed-alignments/${rfamid}/ViennaRNA-400/locations.bed --overlaps performance/tune/${rfamid}-${alpha}-${size}-${init}-overlap.txt --performance performance/tune/${rfamid}-${alpha}-${size}-${init}.performance.txt
done
done
done
#done < rfam-0.8-gt-10.txt
