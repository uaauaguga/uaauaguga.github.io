source := seed-alignments
method := ViennaRNA
flanking := 400
outdir := benchmark/datasets/$(source)/$(rfamid)/$(method)-$(flanking)
#PY27 := /BioII/lulab_b/jinyunfan/anaconda3/envs/shaker-env-py27/bin/python 
identity_filter := -f 0.8
sampling := -s 50

all: $(outdir)/genextreme.shape $(outdir)/folded-soft-constraint.dot $(outdir)/pairing-probability-with-shape.txt $(outdir)/pairing-probability-without-shape.txt


benchmark/motif-structure/$(source)/$(rfamid).dot: benchmark/alignments/$(source)/$(rfamid).stk
	mkdir -p benchmark/alignments/$(source)
	@echo -e "Convert stk to dot bracket notation ..."
	scripts/stk2dbn.py -i $^ -o $@ $(identity_filter) $(sampling)
	@echo -e "Done ."

$(outdir)/constraint $(outdir)/sequence.fa $(outdir)/locations.bed : benchmark/motif-structure/$(source)/$(rfamid).dot
	@echo -e "Add flanking sequence ..."
	mkdir -p $(outdir)
	scripts/add-flanking-sequence.py -i $^ --output $(outdir)/constraint --flanking-length $(flanking) -of $(method) --fasta $(outdir)/sequence.fa --locations $(outdir)/locations.bed
	@echo -e "Done ."

$(outdir)/folded-hard-constraint.dot: $(outdir)/constraint $(outdir)/sequence.fa
	@echo -e "Fold sequence ..."
	scripts/fold.py --fasta $(outdir)/sequence.fa -m $(method) --constraint $(outdir)/constraint --tmp-dir tmp/$(source)-$(rfamid)-$(method)-$(flanking)  --output $(outdir)/folded-hard-constraint.dot
	@echo -e "Done ."


$(outdir)/genextreme.shape: $(outdir)/folded-hard-constraint.dot
	@echo  -e "Simulate SHAPE reacitivity with shaker ..."
	#$(PY27) scripts/simulate-reactivity.py -i $^ -o $@
	scripts/simulate-genextreme.py -i $^ -s $@ -m data/reactivity/patteRNA-weeks-genextreme-model.pkl
	@echo -e "Done ."

$(outdir)/folded-soft-constraint.dot: $(outdir)/genextreme.shape $(outdir)/sequence.fa
	@echo -e "Fold using simulated reactivity ..."
	scripts/fold.py --fasta $(outdir)/sequence.fa -m $(method) --shape $(outdir)/genextreme.shape  --tmp-dir tmp/$(source)-$(rfamid)-$(method)-$(flanking) --output $(outdir)/folded-soft-constraint.dot


$(outdir)/pairing-probability-with-shape.txt: $(outdir)/genextreme.shape $(outdir)/sequence.fa
	@echo -e "Get pairing probability with shape constraint"
	scripts/pairing-probability.py --fasta $(outdir)/sequence.fa -m $(method) --shape $(outdir)/genextreme.shape --tmp-dir $(outdir)/bppm-with-shape --jobs 4 --output $(outdir)/pairing-probability-with-shape.txt

$(outdir)/pairing-probability-without-shape.txt: $(outdir)/sequence.fa
	@echo -e "Get pairing probability without  shape constraint"
	scripts/pairing-probability.py --fasta $(outdir)/sequence.fa -m $(method) --shape $(outdir)/genextreme.shape --tmp-dir $(outdir)/bppm-without-shape --jobs 4 --output $(outdir)/pairing-probability-without-shape.txt

