```{bash}
scripts/structure-accuracy.py -r test/data/patteRNA-weeks.dot -p test/predicted/RNAstructure-real-shape.dot --performance test/performance/RNAstructure-real-shape.txt
```

scripts/fold.py -m RNAstructure -f test/data/patteRNA-weeks.fa  -s test/reactivity-simulation/shaker-loo.shape -o test/predicted/RNAstructure-shaker-shape.dot -t tmp/RNAstructure-shaker-shape

scripts/structure-accuracy.py -r test/data/patteRNA-weeks.dot -p test/predicted/RNAstructure-shaker-shape.dot --performance test/performance/RNAstructure-shaker-shape.txt
