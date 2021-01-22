source := seed-alignments
method := ViennaRNA
flanking := 100
outdir := benchmark/datasets/$(source)/$(rfamid)/$(method)-$(flanking)
PY27 := /BioII/lulab_b/jinyunfan/anaconda3/envs/shaker-env-py27/bin/python 


all: $(outdir)/shaker.shape


benchmark/motif-structure/$(source)/$(rfamid).dot: benchmark/alignments/$(source)/$(rfamid).stk
	@echo -e "Convert stk to dot bracket notation ..."
	scripts/stk2dbn.py -i $^ -o $@
	@echo -e "Done ."

$(outdir)/constraint $(outdir)/sequence.fa $(outdir)/locations.bed : benchmark/motif-structure/$(source)/$(rfamid).dot
	@echo -e "Add flanking sequence ..."
	mkdir -p $(outdir)
	scripts/add-flanking-sequence.py -i $^ --output $(outdir)/constraint --flanking-length $(flanking) -of $(method) --fasta $(outdir)/sequence.fa --locations $(outdir)/locations.bed
	@echo -e "Done ."

$(outdir)/folded.dot: $(outdir)/constraint $(outdir)/sequence.fa
	@echo -e "Fold sequence ..."
	scripts/fold.py --fasta $(outdir)/sequence.fa -m $(method) --constraint $(outdir)/constraint --tmp-dir tmp/$(source)-$(rfamid)-$(method)-$(flanking) --enforce --output $(outdir)/folded.dot
	@echo -e "Done ."


$(outdir)/shaker.shape: $(outdir)/folded.dot
	@echo  -e "Simulate SHAPE reacitivity with shaker ..."
	$(PY27) scripts/simulate-reactivity.py -i $^ -o $@
	@echo -e "Done ."
