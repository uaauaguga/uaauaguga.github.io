#!/usr/bin/env python
import argparse
import tqdm
import re
import numpy as np
import os 
import sys



def parseDBN(s):
    """
    Get base pairs and unpaired positions from dot bracket notations
    """
    s = s.strip()
    L = len(s)
    position_stack = []
    pairs = []
    unpaired = []
    for j in range(L):
        if s[j] == "(":
            position_stack.append(j)
        elif s[j] == ")":
            i = position_stack.pop()
            pairs.append((i,j))
        else:
            unpaired.append(j)
    assert len(position_stack) == 0
    return pairs,unpaired
    

def loadPairsFromDBN(path):
    """
    Load base pairing information from dot bracket file
    Input: dot bracket file path, may contain multiple records
    Sequence and structure should take one line
    Output: A dict, sequence name as key, 
    tuple of (sequence, list of base pair (i,j), unpaired position k) as value
    """
    structures = {}
    i = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if i%3 == 0:
                name = line.split(" ")[0].replace(">","")
            elif i%3 == 1:
                sequence = line
            else:
                pairs,unpaired = parseDBN(line)
                structures[name] = (sequence,pairs,unpaired)
            i += 1
    return structures
            
            
def loadFasta(path):
    """
    Load fasta file into an sequence dict
    """
    sequences = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if len(line)==0:
                continue
            if line.startswith(">"):
                seqid = line.replace(">","").strip() 
                sequences[seqid] = ""
            else:
                sequences[seqid] += line
    return sequences


def adjustPairs(pairs,unpaired,offset):
    """
    Shift the base pairing with a given offset
    Input: List if base pairs, unpaired positions, and off set
    Output: Shifted positions
    """
    adjuated_pairs = []
    adjusted_unpaired = []
    for x,y in pairs:
        x_adjusted = int(x) + int(offset)
        y_adjusted = int(y) + int(offset)
        adjuated_pairs.append((x_adjusted,y_adjusted))
    for ss in unpaired:
        adjusted_unpaired.append(int(ss)+offset)
    return adjuated_pairs,adjusted_unpaired



def checkPairing(seq,pairs):
    """
    Removing non-canonical base pairing
    Input:
    RNA sequence, List of base pairing tuple
    Output:
    Filtered list of base pairing tuple
    """
    canonical = set([("A","U"),("A","T"),
                    ("C","G"),
                    ("G","U"),("G","T"),("G","C"),
                    ("T","A"),("T","G"),
                    ("U","A"),("U","G")])
    pairs_filtered = set()
    n_noncanonical = 0
    for x,y in pairs:
        pair = seq[int(x)],seq[int(y)]
        if pair  in canonical:
            pairs_filtered.add((x,y))
        else:
            print(x,y,"\t",pair[0],pair[1])
            n_noncanonical += 1
    if n_noncanonical > 0:
        print("{} non conanical base pairing in the seed structure was removed".format(n_noncanonical))
    return pairs_filtered
        



def makeRNAstructureConstraint(pairs,unpaired):
    """
    Make constrain file for Fold in RNA structure package
    """
    const = ""
    const += "DS:\n-1\n" # No double strand constraint
    const += "SS:\n"     # Single strand constraint
    for ss in unpaired:
        const += "{}\n".format(ss+1)
    const += "-1\n"
    const += "Mod:\n-1\n"
    const += "Pairs:\n"  # Pairing constraint
    for x,y in pairs:
        const += "{} {}\n".format(x+1,y+1)
    const += "-1 -1\n"
    const += "FMN:\n-1\n"
    const += "Forbids:\n-1 -1"
    return const
    
def makeViennaRNAConstraint(pairs,unpaired,length):
    """
    Make constraint for RNAfold in ViennaRNA package
    . No constraint
    x Force unpaired
    ( Force left pair
    ) Force right pair
    < Force to be paired, and at left side
    > Force to be paired, and at right side
    """
    const = list(length*".")
    for x,y in pairs:
        if x > y:
            x, y = y, x
        const[x] = "("
        const[y] = ")"
    for ss in unpaired:
        const[ss] = "x"
    return "".join(const)


def getRandomSequence(L):
    """
    Generate fully random sequence of given length
    """
    alphabet = "ACGT"
    seq = ""
    for i in range(L):
         seq += alphabet[np.random.randint(4)]
    return seq



def main():
    parser = argparse.ArgumentParser(description='Add flanking sequence. If the structure information is provided in dot bracket format or ct format, the output could be ')
    parser.add_argument('--input','-i',required=True, help="Input path. For fasta and dbn input, accept multiple record")
    parser.add_argument('--output', '-o', required=True, help = "Output path for fasta / ViennaRNA constraint, output directory for RNAstructure constraint ")
    parser.add_argument('--in-format','-if', type=str, default="dbn", choices = ["dbn","fasta"],help="Input format")
    parser.add_argument('--out-format','-of', type=str, default="ViennaRNA", choices = ["ViennaRNA","RNAstructure","fasta"], help="Output format")
    #parser.add_argument('--database','-db', type=str,help="Data base of full length sequence. If not none, first attemp to retrieve sequence from this fasta file instead of random generate flanking sequence")
    parser.add_argument('--flanking-length', '-fl', type=int, help = "Flanking length of simulated sequence")
    parser.add_argument('--seed','-s',type=int,default=777, help = "Seed for random flanking sequencing generation")
    parser.add_argument('--fasta', '-f', type=str, help = "Whether / Where to write fasta file (if the specified output format is not fasta format)" )
    parser.add_argument('--locations','-loc',type=str,required=True,help = "Location of the motif in the flanked sequence")
    args = parser.parse_args()

    np.random.seed(args.seed)
   

    if not os.path.exists(args.input):
        print("Error, {} does not exist.".format(args.input))
        sys.exit(1)
    
    if args.in_format == "dbn":
        print("Input in dot bracket format .")
    elif args.in_format == "fasta":
        print("Input is fasta file .") 


    if args.out_format == "RNAstructure" or args.out_format == "ViennaRNA":
        outdir = args.output
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        print("Write {} constraint to {}".format(args.out_format,outdir))

    if args.fasta is not None and args.out_format != "fasta":
        f_fasta = open(args.fasta,"w")
        print("Fasta file will also be generated.")
    if args.out_format == "fasta":
        f_fasta = open(args.output,"w")
        print("Output fasta file with flanking sequence.")

    fbed = open(args.locations,"w")


    if args.in_format == "dbn":
        structureDict = loadPairsFromDBN(args.input)
        for name,record in structureDict.items():
            print("Processing",name,"...")
            sequence, pairs, unpaired = record
            pairs = checkPairing(sequence, pairs)
            offset = np.random.randint(args.flanking_length)
            pairs, unpaired = adjustPairs(pairs,unpaired,offset)
            flanking_sequence = getRandomSequence(args.flanking_length)
            sequence_flanked = flanking_sequence[:offset] + sequence + flanking_sequence[offset:]
            start, end = str(offset),str(offset + len(sequence))
            name = name + ":" + start + "-" + end
            fbed.write("\t".join([name,start,end])+"\n")
            path = os.path.join(outdir,name + ".CON")
            if args.out_format == "RNAstructure":
                constraint = makeRNAstructureConstraint(pairs,unpaired)
                with open(path,"w") as fout:
                    fout.write(constraint)
            elif args.out_format == "ViennaRNA":
                constraint = makeViennaRNAConstraint(pairs,unpaired,len(sequence_flanked))
                with open(path,"w") as fout:
                    fout.write(">" + name + "\n")
                    fout.write(sequence_flanked.replace("U","T") + "\n")
                    fout.write(constraint + "\n")
            if args.fasta is not None or args.out_format == "fasta":
                f_fasta.write(">" + name + "\n")
                f_fasta.write(sequence_flanked.replace("U","T") + "\n")
    elif args.in_format == "fasta":
        sequenceDict = loadFasta(args.input)
        for name,sequence in sequenceDict.items():
            flanking_sequence = getRandomSequence(args.flanking_length)
            offset = np.random.randint(args.flanking_length)
            sequence_flanked = flanking_sequence[:offset] + sequence + flanking_sequence[offset:]
            start, end = str(offset),str(offset + len(sequence))
            name = name + ":" + start + "-" + end
            fbed.write("\t".join([name,start,end])+"\n")
            f_fasta.write(">" + name + "\n")
            f_fasta.write(sequence_flanked.replace("U","T") + "\n")
    if args.fasta is not None or args.out_format == "fasta":
        f_fasta.close()


if __name__ == "__main__":
    main()
