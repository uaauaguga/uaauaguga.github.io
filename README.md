## Scripts for RNA Analysis

scripts/add-flanking-sequence.py --input RF00032.dot --in-format dbn --output RF00032.CON --out-format ViennaRNA --flanking-length 100 --fasta RF00032.fa

scripts/build-index.py --fasta ${path} --aligner star --prefix data/reference-structure/model-organism/${data}/star-index/${organism}/ --tmp-dir tmp/star/${data}-${organism} > log/build-star-index/${data}-${organism}.log
