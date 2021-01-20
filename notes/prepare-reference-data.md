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
