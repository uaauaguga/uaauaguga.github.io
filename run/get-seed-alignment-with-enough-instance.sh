#!/bin/bash
for rfamid in $(ls benchmark/datasets/seed-alignments);do
n=$(cat benchmark/datasets/seed-alignments/${rfamid}/ViennaRNA-100/locations.bed | wc -l | sed 's/\n//g')
L=$(head -n 1  benchmark/datasets/seed-alignments/${rfamid}/ViennaRNA-100/locations.bed | awk '{print $3-$2}')
if [[ "$n" -ge "10" ]];then
echo -e "${rfamid}\t${L}"
fi
done
