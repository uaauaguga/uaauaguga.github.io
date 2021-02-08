## Notes for benchmark data preparation
- Split Rfam seed alignments
```{bash}
scripts/separate-seed.py -i data/Rfam/seed-alignments-concatenated/Rfam.seed -o data/Rfam/seed-alignments
scripts/get-stk-length.py -i data/Rfam/seed-alignments -o data/Rfam/info/average-length.txt
cat data/Rfam/info/average-length.txt | awk 'NR>1&&$2<50&&$3>10{print $1}' > data/Rfam/info/rfam-id-st50-lt10.txt
cat data/Rfam/info/average-length.txt | awk 'NR>1&&$2<50&&$3>5{print $1}' > data/Rfam/info/rfam-id-st50-lt5.txt
```

## Download full alignment
- Currently only consider motif with average length shorter than 50, and more than 10 seed sequences

```{bash}
cripts/stk2dbn.py -i data/Rfam/seed-alignments/RF00032.stk -o test/RF00032.dot
scripts/add-flanking-sequence.py -i RF00032.dot --output RF00032.RNAstructure.const --flanking-length 100 -of RNAstructure --fasta RF00032.RNAstructure.fa
scripts/add-flanking-sequence.py -i test/RF00032.dot --output test/RF00032.ViennaRNA.const --flanking-length 100 -of ViennaRNA --fasta test/RF00032.ViennaRNA.fa

scripts/fold.py --fasta  test/RF00032.ViennaRNA.fa -m ViennaRNA --constraint test/RF00032.ViennaRNA.const --tmp-dir test/RF00032.ViennaRNA.tmp --enforce --output test/RF00032.ViennaRNA.dot
```


## Curation of reference structure of some structural noncoding RNA 
- Get sequence id in seed alignemnts `run/get-seed-ids.sh`
- Get reference structure
```{bash}
run/get-reference-structure-of-model-organism.py --output data/reference-structure/model-organism/homo-sapien-reference-structure.dot --taxid 9606
run/get-reference-structure-of-model-organism.py --output data/reference-structure/model-organism/mus-musculus-structure.dot --taxid 10090
run/get-reference-structure-of-model-organism.py --output data/reference-structure/model-organism/saccharomyces-cerevisiae.dot --taxid 4932
run/get-reference-structure-of-model-organism.py --output data/reference-structure/model-organism/oryza-sativa.dot --taxid 4530
run/get-reference-structure-of-model-organism.py --output data/reference-structure/model-organism/arabidopsis-thaliana.dot --taxid 3702
for file in $(ls dbn);do cat dbn/${file} | grep -v '^[\.\(]' > fasta/${file%%.*}.fa ;done
mkdir -p fasta-filtered
for file in $(ls fasta);do cd-hit -c 0.7 -i fasta/${file} -o fasta-filtered/${file};done
```


