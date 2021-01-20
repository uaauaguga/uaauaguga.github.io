#!/usr/bin/env python
import argparse
from Bio import AlignIO

def checkDBN(s):
    stack = []
    for c in s:
        if c == "(":
            stack.append(c)
        elif c == ")":
            _ = stack.pop()
        else:
            continue
    assert len(stack) == 0

def main():
    parser = argparse.ArgumentParser(description='Convert stockholm format to dot bracket format.')
    parser.add_argument('--input', '-i',required=True,help="Input path, read from stdin by default.")
    parser.add_argument('--output','-o',required=True,help="Output path, print to stdout by default.")
    args = parser.parse_args()

    align = AlignIO.read(args.input, "stockholm")
    deWUSS = dict(zip("[{<-._,:>}]","(((.....)))"))
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
    ## Assign structure to each sequence
    for record in align:
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
        name  = "|".join([record.annotations["accession"],str(record.annotations["start"]),str(record.annotations["end"])])
        
        fout.write(">" + name + "\n")
        fout.write(ungapped_sequence + "\n")
        fout.write(ungapped_structure + "\n")
    fout.close()
    
if __name__ == "__main__":
    main()