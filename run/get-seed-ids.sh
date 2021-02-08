#!/bin/bash
(for stk in $(ls data/Rfam/seed-alignments);do
esl-reformat fasta data/Rfam/seed-alignments/${stk} | grep '>' | sed 's/>//' | awk -F '.' '{print $1}'
done) | sort | uniq > seed-seq-ids.txt
