---

layout: post
title:  "Play with Genomic Features "
date:   2021-03-20 1:11:20 +0800
usemathjax: false
categories: jekyll update
---

- Genomic data is often genomic ranges associate with some metadata
- gtf, bed, vcf, sam, wig, bedgraph ...
- Multiple tools is available for manipulate such genomic ranges

### A closer look at file formats

- The coordinate system
- gtf / gff file
- bed / bedgraph and bed12 format
  - bed: genome interval associate with some metadata
  - bed12: a specific bed, each line specify the structure of a transcript
- vcf format
- wiggle and bigwig
- sam

### Parse / Load different file format & file format conversion

- R / bioconductor, https://bioconductor.org/developers/how-to/commonImportsAndClasses/
  - [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) support data importing in multiple format
    - "In summary, for the typical use case of combining gene models with experimental data, GFF is preferred for gene models and `BigWig` is preferred for quantitative score vectors. "
  - https://kasperdanielhansen.github.io/genbioconductor/html/rtracklayer_Import.html
  - [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) is useful for manipulate genome annotation

```r
  library(rtracklayer)
  tx <- rtracklayer::import(file.bed12,"bed")
```

```r
  library(GenomicFeatures)
  txdb <- makeTxDbFromGFF(file=gtf.path,format="gtf")
  # Output:
  #Import genomic features from the file as a GRanges object ... OK
  #Prepare the 'metadata' data frame ... OK
  #Make the TxDb object ... The "phase" metadata column contains non-NA values for features of type stop_codon. 

  # Get exons by tx
  # Return a list of GRanges
  tx.exons <- exonsBy(txdb, by="tx")
  #  information was ignored.OK
```


### Get transcript from genome with bed12/gtf gene model

- Command line tools

  - <http://ccb.jhu.edu/software/stringtie/dl/gffread-0.11.4.Linux_x86_64.tar.gz>

  - Get transcripts fasta from genome fasta

```bash
    #gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] [-o <outfile>] [-t <trackname>] [-r [[<strand>]<chr>:]<start>..<end> [-R]][-CTVNJMKQAFPGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>][-i <maxintron>] [--bed] [--table <attrlist>] [--sort-by <refseq_list.txt>]
    gffread -g genome.fa -s genome.size -W -M -F -G -A -O -E -w transcriptome.fa -d transcriptome.collapsed.info genome.gtf
    #  -W  for -w and -x options, also write for each fasta record the exon coordinates projected onto the spliced sequence
    #  -M/--merge : cluster the input transcripts into loci, collapsing matching transcripts (those with the same exact introns and fully contained)
    # -F  preserve all GFF attributes (for non-exon features)
    # -G  do not keep exon attributes, move them to the transcript feature (for GFF3 output)
    # -A  use the description field from <seq_info.fsize> and add it as the value for a 'descr' attribute to the GFF record
    # -O  process also non-transcript GFF records (by default non-transcript records are ignored)
    # -E  expose (warn about) duplicate transcript IDs and other potential problems with the given GFF/GTF records
    #  -w  write a fasta file with spliced exons for each GFF transcript
```

  - Convert gff to bed12

```bash
    gffread -W -M -F -G -A -E --bed {gtf} > {bed}
```

  - Get fasta from bed12 file

    - The result is strange, need check (seems intron is not removed)
    
```bash
    bedtools getfasta -name+ -fi genome.fasta -bed {bed12} -split -s > tx.fa
```

- Bioconductor packages

  - This implementation is quite slow, seems can be improved

    
```R
  library(Biostrings)
  library(GenomicFeatures)
  library(pbapply)
  
  # Load gene model from gtf file
  txdb <- makeTxDbFromGFF(file=gtf.path,format="gtf")
  tx.exons <- exonsBy(txdb, by="tx")
  
  # Load genome sequence from fasta file
  fasta.path <- "Path to genomic fasta file"
  genome.bs <- Biostrings::readDNAStringSet(fasta.path)
  
  # Extract tx sequence
  getTxSequence <- function(tx.grs,genome.dna.set){
  chrom <- seqnames(tx.grs)[1]
  strand <- strand(tx.grs)[1]
  dna <- genome.dna.set[[chrom]]
  tx.name <- gsub(".exon.*$","",tx.grs$exon_name[1])
  tx.irs <- IRangesList(ranges(tx.grs))
  names(tx.irs) <- tx.name
  tx.seq <-   GenomicFeatures::extractTranscriptSeqs(dna,transcripts=tx.irs,strand=strand)
  tx.seq
  }
  
  tx.sequences <- unlist(pblapply(tx.exons,100,getTxSequence,genome.dna.set=genome.bs))
  names(tx.sequences) <- NULL
  tx.sequences <- do.call(c,tx.sequences)
  Biostrings::writeXStringSet(tx.sequences,"tx.fasta")
```

- Perl scripts
  
  - Check this <https://metacpan.org/release/Bio-ViennaNGS>, seems not work for spliced gene



### Liftover between different version of genome
- [CrossMap](http://crossmap.sourceforge.net/)
- [rtracklayer::liftOver](https://bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html)


### Projection between genome coordinate and transcript coordinate

- You have a list of genome interval, and gene model, you want to project the interval to the coordinate of transcript

- Bioconductor implementation

  - [GenomicFeatures](https://rdrr.io/bioc/GenomicFeatures/) package

    - `mapToTranscripts`

    - `mapFromTranscripts`

    - `transcriptLocs2refLocs `

      

```R
    library(GenomicFeatures)
    gtf.path <- "annotation/gene-models/gtf/Arabidopsis_thaliana.TAIR10.46.gtf"
    txdb <- makeTxDbFromGFF(file=gtf.path,format="gtf")
    tx.exons <- exonsBy(txdb, by="tx",use.names=TRUE)
    tx.used <- tx.exons["AT1G01010.1"]
    gr <- GRanges(c("1","1","2","2"),IRanges(c(4000,4160,2000,3000), width=100),c("+","+","-","+"))
    names(gr) <- rep("AT1G01010.1",length(gr))
    mapToTranscripts(gr,tx.used)
```


### Get CDS relative to transcript coordinate

```R
library(GenomicFeatures)
gtf.path <- "annotation/gene-models/gtf/Arabidopsis_thaliana.TAIR10.46.gtf"
txdb <- makeTxDbFromGFF(file=gtf.path,format="gtf")
tx.utr5p <- fiveUTRsByTranscript(txdb,use.names=TRUE)
tx.cds <- cdsBy(txdb,by="tx",use.names=TRUE)
utr5p.lengths <- sum(width(tx.utr5p))
cds.lengths <- sum(width(tx.cds))
utr5p.add <- rep(0,length(tx.no5putr))
tx.no5putr <- setdiff(names(cds.lengths),names(utr5p.lengths))
names(utr5p.add) <- tx.no5putr
utr5p.lengths <- c(utr5p.lengths,utr5p.add)
utr5p.lengths <- utr5p.lengths[names(cds.lengths)]
tx.coordinates <- data.frame(utr5p.lengths,cds.lengths)
tx.coordinates[["utr5p-start"]] <- 0
tx.coordinates[["utr5p-end"]] <- tx.coordinates[["utr5p.lengths"]]
tx.coordinates[["cds-end"]] <- tx.coordinates[["utr5p.lengths"]] + tx.coordinates[["cds.lengths"]]
utr5p <- tx.coordinates[,c("utr5p-start","utr5p-end")]
utr5p <- utr5p[utr5p[["utr5p-end"]]>0,]
write.table(utr5p,file="5putr.bed",sep="\t",quote = F, col.names = F)
write.table(tx.coordinates[,c("utr5p-end","cds-end")],file="cds.bed",sep="\t",quote = F, col.names = F)

```


#### Useful tools and resource

- Useful tools
  
  - [bedtools](https://bedtools.readthedocs.io/en/latest/)
  
- Useful links
  
  - Manipulate genomic interval
    - <https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GRanges_and_GRangesList_slides.pdf>
    - <http://bioconductor.org/help/course-materials/2010/EMBL2010/GenomicRanges.pdf>
    - <http://bioconductor.org/packages/release/bioc/vignettes/HelloRanges/inst/doc/tutorial.pdf>
    - <http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf>
    - <http://bioconductor.org/help/course-materials/2014/CSAMA2014/2_Tuesday/lectures/Ranges_Sequences_and_Alignments-Lawrence.pdf>
    - <https://www.bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops1/Obenchain/Lecture-sequences.html>

### Some task related to interval comparison

#### Take intersection of feature sets
- Annotate input features with known feature set
  - `bedtools annotate`
  - `bedtools intersect`

- Count coverage

```bash
bedtools coverage -a features.bed -b input.bam -counts > counts.txt
```

- Take some arithmetic operation across genomic interval

```bash
# -o can be arithmetic operation like sum, min, max, mean, median ...
bedtools map -a interval.tosummarize.bed -b input.value.bed  -c 4 -o sum > output.bed
```


#### Take union

- Merge bed file

```bash
# Consider strandness
sort -k1,1 -k2,2n ${input} | bedtools merge -s -c 6 -o distinct | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,".",".",$4}' > ${output}
```

- Merge / reduce exon from different isoforms from a gene
  
  - Get gene length for TPM or FPKM calculation (If you use tools like featureCount, gene length information is provided in the output)
```R
#!/usr/bin/env Rscript
# Refer to https://www.biostars.org/p/83901/
library(GenomicFeatures)
gtf.path <- "annotation/gene-models/gtf/Arabidopsis_thaliana.TAIR10.46.gtf"
txdb <- makeTxDbFromGFF(file=gtf.path,format="gtf")
gene.exons <- exonsBy(txdb, by="gene")
# The key operation is GenomicRanges::reduce
gene.lengths <- sum(GenomicRanges::width(GenomicRanges::reduce(gene.exons)))
write.table(gene.lengths,"gene.length.txt",col.name=F,quote=F,sep="\t")
```

#### Take difference

- Get intron location from gff/gtf file
- `bedtools subtract`

### Algorithm behind these tools

- Interval tree
- [Nested Containment List](https://academic.oup.com/bioinformatics/article/23/11/1386/199545)

### Implement similar function in other programming language

- perl
  - [Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree)

```perl
#!/usr/bin/env perl
use Set::IntervalTree;
use Getopt::Long;

GetOptions ("database=s" => \my $db_path,   
            "query=s"   => \my $query_path)      
or die("Error in parsing command line arguments\n");

open(FDB,"<",$db_path) or die("Failed opening database file\n");
open(FQUERY,"<",$query_path) or die("Failed opening query file\n");


my $treedb = Set::IntervalTree->new;


#print "Load database file ...\n";
my $chr,$start,$end; #,$name;
while(<FDB>){
  chomp;
  @fields = split "\t",$_;
  $chr = $fields[0];
  $start = $fields[1];
  $end = $fields[2];
  $treedb->insert($_."",$start,$end);
}
#print "Done .\n";


while(<FQUERY>){
  chomp;
  @fields = split "\t",$_;
  $start = $fields[1];
  $end = $fields[2];
  $entries = $treedb->fetch($start,$end);
  for my $entry (@$entries){
    print $entry,"\n";
  }
}

close(FDB);
close(FQUERY);
```

- python
  - There are also interval tree implementations in python [intervaltree](https://github.com/chaimleib/intervaltree), not try yet
  - Seems [pyranges](https://github.com/biocore-ntnu/pyranges) is a good alternative in python 
  - [bx-python](https://github.com/bxlab/bx-python)

- C++
  - [cgranges](https://github.com/lh3/cgranges)

