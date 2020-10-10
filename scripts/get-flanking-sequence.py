import argparse
from tqdm import tqdm
import numpy as np
import pyfaidx
parser = argparse.ArgumentParser(description='Get flanking sequence')
parser.add_argument('--query','-q',type=str,required=True,help="gbid,start,end")
parser.add_argument('--database', '-db', type=str,help='input full length sequence',default="Rfam/seeds/CMfinder-set-full-length.fasta")
parser.add_argument('--output','-o',type=str,required=True,help='samples sequence with flanking regions')
parser.add_argument('--length','-l',type=int,default=200)
parser.add_argument('--random-state','-rs',type=int,default=777)
args = parser.parse_args()


faObj = pyfaidx.Fasta(args.database)
np.random.seed(args.random_state)

queries = open(args.query)
f = open(args.output,"w")

for query in tqdm(queries):
    gbid,start0,end0 = query.strip().split("\t")
    ## Get position of seed alignment on 0-based coordinate of full length sequence
    start0 = int(start0)
    end0 = int(end0)
    if start0 < end0:
        strand = "+"
        start,end = start0,end0
    else:
        strand = "-"
        start,end = end0,start0
    start -= 1 # Convert to zero-based coordiante
    seed_length = end - start
    if gbid not in faObj.keys():
        print("{} is not present, skipping".format(gbid))
        continue
    length = len(faObj[gbid][:])
    # Ramdomly generate left flanking length
    left_length = np.random.randint(0,args.length) if args.length > 0 else 0
    # Left position, 0-based coordiante, relative to full length sequence
    left_pos = start - left_length if start >= left_length else start
    # Remaining sequence length on the right side
    remain = length - end
    # Right position, 0 based coordiante
    right_length = min(args.length - left_length,remain)
    right_pos = end + right_length
    seed_start = left_length
    seed_end = left_length+seed_length
    seqIdOut = gbid + ":" + "-".join([str(start0),str(end0)]) + ":" + "-".join([str(seed_start),str(seed_end)])
    seq = faObj[gbid][left_pos:right_pos]
    if strand == "-":
        seq = seq.reverse.complement
    print(">"+seqIdOut,file=f)
    print(seq.seq,file=f)
f.close()
    


