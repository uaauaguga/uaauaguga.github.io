#!/usr/bin/env python
import re 
import argparse
def main():
    parser = argparse.ArgumentParser(description='Convert transterm prediction to bed format')
    parser.add_argument('--input', '-i', type=str, required=True, help='input transterm prediction')
    parser.add_argument('--output','-o',type=str, required=True, help="output bed file")
    args = parser.parse_args()
    """
    SEQUENCE SRR5947820:k119_84915:0-6619  (length 6619)

    PF10090.12       1051 - 2040     + | 
      TERM 1         2060 - 2087     + F   100  -8.4 -5.32599 | opp_overlap 2060, overlap 2062
      AGTATCATAAAAAAA        GACGGACACTTC CAAA GAAATGTCCGCC         TTTTTCTATTTTATC

        TERM 19  15310 - 15327  -      F     99      -12.7 -4.0 |bidir
    (name)   (start - end)  (sense)(loc) (conf) (hp) (tail) (notes)
    """
    fout = open(args.output,"w")
    with open(args.input,encoding = "ISO-8859-1") as f:
        for line in f:
            if line.startswith("SEQUENCE"):
                line = line.strip()
                contig_id = re.split("\s+",line)[1]
            elif line.startswith("  TERM"):
                line = line.strip() 
                fields = re.split("\s+",line)
                start, end = int(fields[2]), int(fields[4])
                start -= 1
                if start > end:
                    end, start = start, end
                strand = fields[5]
                score = fields[7]
                print(contig_id,start,end,".",score,strand,sep="\t",file=fout)
    fout.close()
                    



if __name__ == "__main__":
    main()
