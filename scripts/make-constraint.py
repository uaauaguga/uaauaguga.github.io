import argparse
import tqdm
import re
import os 
import sys

def load_pairs(path,preserveSS=True):
    """Load base pairs and sequence from .ct file"""
    pairs = []
    unpaired = []
    firstLine = True
    sequence = ""
    pairedPos = set()
    with open(path) as f:
        for line in f:
            if firstLine:
                firstLine = False
                continue
            line = line.strip()
            if len(line)==0:
                continue
            fields = re.split(r"\s+",line)
            x,y = int(fields[0]),int(fields[4])
            if y==0:
                unpaired.append(x)
            else:
                pairedPos.add(y)
                if x not in pairedPos:
                    pairs.append((x,y))
            sequence += fields[1] 
    return sequence,pairs,unpaired

def adjustPairs(pairs,unpaired,offset):
    adjuated_pairs = []
    adjusted_unpaired = []
    for x,y in pairs:
        x_adjusted = int(x) + int(offset)
        y_adjusted = int(y) + int(offset)
        adjuated_pairs.append((x_adjusted,y_adjusted))
    for ss in unpaired:
        adjusted_unpaired.append(int(ss)+offset)
    return adjuated_pairs,adjusted_unpaired


def makeRNAstructureConstraint(pairs,unpaired):
    """
    Make constrain file in RNA structure" 
    """
    const = ""
    const += "DS:\n-1\n"
    const += "SS:\n"
    for ss in unpaired:
        const += "{}\n".format(ss)
    const += "-1\n"
    const += "Mod:\n-1\n"
    const += "Pairs:\n"
    for x,y in pairs:
        const += "{} {}\n".format(x,y)
    const += "-1 -1\n"
    const += "FMN:\n-1\n"
    const += "Forbids:\n-1 -1"
    return const
    
def makeViennaRNAConstraint(pairs,unpaired,length):
    """
    . No constraint
    x Force unpaired
    < Force left pair
    > Force right pair
    """
    const = length*"."
    for x,y in pairs:
        x_,y_ = x-1,y-1 if x < y else y-1,x-1
        const[x_] = "<"
        const[y_] = ">"
    for ss in unpaired:
        const[ss-1] = "x"
    return const


def main():
    #genebank-id:rfam 1based start-end:seed start relative to flanked sequence 0based-end
    parser = argparse.ArgumentParser(description='Prepare constraint file for Fold')
    parser.add_argument('--ct-file','-ct',required=True,help="Ct file of the seed alignment")
    parser.add_argument('--start','-s',required=True,type=int,help="Start of the seed, 0 based relative to sampled sequence")
    parser.add_argument('--end','-e',required=True,type=int,help="End of the seed, 0 based relative to sampled sequence")
    parser.add_argument('--output','-o',required=True,help="Output constraint file")
    parser.add_argument('--preserve-single','-ps',action="store_true",help="Keep both double strand and single strand in the structure",default=True)
    args = parser.parse_args()
    start,end = args.start,args.end
    output = args.output
    offset = start
    seed_length = end - start
    ctPath = args.ct_file
    if not os.path.exists(ctPath):
        print("Error, {} does not exist.".format(ctPath))
        sys.exit(1)
    else:
        sequence,pairs,unpaired = load_pairs(ctPath,preserveSS=args.preserve_single)
        if len(sequence)!=seed_length:
            print("Inconsistent seed length for ct file and seed length")
            sys.exit(2)
        pairs_adjusted,unpaired_adjusted = adjustPairs(pairs,unpaired,int(offset))
        f = open(output,"w")
        makeConstraint(pairs_adjusted,unpaired_adjusted,f)
        f.close()
