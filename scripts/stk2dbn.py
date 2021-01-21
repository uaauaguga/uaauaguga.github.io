#!/usr/bin/env python
import argparse
from Bio import AlignIO

def checkDBN(s):
    """
    Check whether a dot bracket notation represents a valid secondary structure
    """
    stack = []
    for c in s:
        if c == "(":
            stack.append(c)
        elif c == ")":
            _ = stack.pop()
        else:
            continue
    assert len(stack) == 0


def getIdentity(s1,s2):
    assert len(s1) == len(s2)
    identical = 0
    L = len(s1)
    for a,b in zip(s1,s2):
        if a == b:
            if a == "-":
                L -= 1
            else:
                identical += 1
    return identical/L


def getSimilarSequence(align,threshold = 0.9):
    """
    Input: an alignment object
    Output: sequence id to be removed
    """
    depleted = set()
    records = [(record.id,record.seq) for record in align]
    n = len(records)
    for i in range(n):
        name1, sequence1 = records[i]
        if name1 in depleted:
            continue
        for j in range(i+1,n):
            name2, sequence2 = records[j]
            if name2 in depleted:
                continue
            similarity = getIdentity(sequence1,sequence2)
            if similarity > threshold:
                depleted.add(name2)
    return depleted,n
    
            
                


def main():
    parser = argparse.ArgumentParser(description='Convert stockholm format to dot bracket format.')
    parser.add_argument('--input', '-i',required=True,help="Input path, read from stdin by default.")
    parser.add_argument('--output','-o',required=True,help="Output path, print to stdout by default.")
    parser.add_argument('--filter','-f',type=float,help="Should be a float between 0 and 1, sequence with identity higher than this value will be filtered")
    args = parser.parse_args()

    align = AlignIO.read(args.input, "stockholm")
    deWUSS = dict(zip("([{<-._,:~>}])","((((......))))"))
    structure_WUSS = align.column_annotations["secondary_structure"]
    L = len(structure_WUSS)
    structure_dbn = ""

    ## Get base pairing in consensus secondary structure
    position_stack = []
    left_pairing = {}
    right_pairing = {}
    for j in range(L):
        c = deWUSS[structure_WUSS[j]]
        structure_dbn += c
        if c == "(":
            position_stack.append(j)
        elif c == ")":
            i = position_stack.pop()
            left_pairing[i] = j
            right_pairing[j] = i
    fout = open(args.output,"w")

    if args.filter is not None:
        depleted,n = getSimilarSequence(align,threshold = args.filter)
        print("Among {} input sequences, {} with high similarity would be removed ".format(n,len(depleted)))
        print("{} sequences were retained".format(n - len(depleted)))

    ## Assign structure to each sequence
    for record in align:
        if args.filter is not None and record.id in depleted:
            #print(record.annotations["accession"],"was filtered.")
            continue
        sequence = record.seq
        current_structure = list(structure_dbn)
        ## Delete base pairing corresponding to insertions
        for i in range(L):
            if sequence[i] == "-":
                if i in left_pairing.keys():
                    current_structure[i],current_structure[left_pairing[i]] = ".","."
                elif i in right_pairing.keys():
                    current_structure[i],current_structure[right_pairing[i]] = ".","."
        ungapped_sequence = ""
        ungapped_structure = ""
        ## Remove insertions in the sequence
        for i in range(L):
            if sequence[i] != "-":
                ungapped_structure += current_structure[i]
                ungapped_sequence += sequence[i]
        checkDBN(ungapped_structure)
        name  = record.annotations["accession"] + ":" + str(record.annotations["start"]) + ":" + str(record.annotations["end"])
        
        fout.write(">" + name + "\n")
        fout.write(ungapped_sequence + "\n")
        fout.write(ungapped_structure + "\n")
    fout.close()
    
if __name__ == "__main__":
    main()
