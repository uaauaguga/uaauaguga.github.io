#!/usr/bin/env python
import argparse
from Bio import AlignIO
import numpy as np
from utilities import checkDBN 

def getIdentity(s1,s2):
    assert len(s1) == len(s2)
    identical = 0
    L = len(s1)
    for a,b in zip(s1,s2):
        if a == b:
            if a == "-" or a == ".":
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
    return depleted
    
            
                


def main():
    parser = argparse.ArgumentParser(description='Convert stockholm format to dot bracket format.')
    parser.add_argument('--input', '-i',required=True,help="Input path, read from stdin by default.")
    parser.add_argument('--output','-o',required=True,help="Output path, print to stdout by default.")
    parser.add_argument('--filter','-f',type=float,help="Should be a float between 0 and 1, sequence with identity higher than this value will be filtered")
    parser.add_argument("--sampling",'-s',type=int,help="Subsampling a given number of sequences")
    parser.add_argument("--number",'-n',type=int,default=1000,help="Subsampling a given number of sequences")
    args = parser.parse_args()

    align = AlignIO.read(args.input, "stockholm")
    deWUSS = dict(zip("([{<-._,:~Aa>}])","((((........))))"))
    structure_WUSS = align.column_annotations["secondary_structure"]
    #print(structure_WUSS)
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

    sequence_names = set([ record.id for record in align ])
    n_total = len(sequence_names)

    depleted = set()

    print("Calculate sequence identity ...")
    if args.filter is not None:
        depleted = getSimilarSequence(align,threshold = args.filter)
        print("Among {} input sequences, {} with high similarity would be removed ".format(n_total,len(depleted)))
        n_passed = n_total - len(depleted)
        print("{} sequences were retained".format(n_passed))
        sequence_names = sequence_names.difference(depleted)
    else:
        n_passed = n_total


    if args.sampling is not None:
        if args.sampling >= n_passed:
            print("The specified number for subsampling is larger than number of available sequences")
            print("Subsmpling will not be performed")
        else:
            print("Sample {} sequence from {} sequences.".format(args.sampling,n_passed))
            n_to_discard = n_passed - args.sampling
            print("{} sequence will further be discarded".format(n_to_discard))
            depleted = depleted.union(set(np.random.choice(list(sequence_names),n_to_discard,replace = False)))
    print("{} sequences were finally used".format(n_total - len(depleted)))
            
    ## Assign structure to each sequence
    print("Convert stockholm format to dot bracket notations ...")
    for record in align:
        if args.filter is not None:
            if record.id in depleted:
                print(record.annotations["accession"],"was filtered.")
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
        if "start" in record.annotations.keys():
            name  = record.annotations["accession"] + ":" + str(record.annotations["start"]) + ":" + str(record.annotations["end"])
        else:
            name  = record.annotations["accession"]
        
        fout.write(">" + name + "\n")
        fout.write(ungapped_sequence + "\n")
        fout.write(ungapped_structure + "\n")
    fout.close()
    print("Done .")
    
if __name__ == "__main__":
    main()
