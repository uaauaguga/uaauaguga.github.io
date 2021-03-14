### Demo motif finder
```bash
benchmark/datasets/seed-alignments/RF00037/ViennaRNA-100/sequence.fa
scripts/pairing-probability.py -f benchmark/datasets/seed-alignments/RF00037/ViennaRNA-100/sequence.fa --output RF00037-100-ViennaRNA.txt -m ViennaRNA --tmp-dir tmp/RF00037-100
scripts/simulate-genextreme.py --input test/motif-finder-demo/data/RF00037-100-hardconstraint-ViennaRNA.dot --output test/motif-finder-demo/data/RF00037-100-hardconstraint-ViennaRNA-genextreme.shape --model patteRNA-weeks-genextreme.pkl
```
