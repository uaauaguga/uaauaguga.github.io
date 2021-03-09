#!/usr/bin/env python
import argparse
import numpy as np
from utilities import checkDBN,loadRecords
import sys
import pandas as pd

def evaluate(predicted_structure,reference_structure):
    """
    Evaluate the structure prediction performance 
    Currently not consider pseudoknots
    """
    pos_stack_pred = []
    pos_stack_ref = []
    ref_pseudoknots_idx = []
    ref_pairs,pred_pairs = set(),set()
    L = len(predicted_structure)
    for i in range(L):
        if reference_structure[i] == "(":
            pos_stack_ref.append(i)
        elif reference_structure[i] == ")":
            ridx = i
            lidx = pos_stack_ref.pop()
            ref_pairs.add((lidx,ridx)) 
        else:
            if reference_structure[i] != ".":
                ref_pseudoknots_idx.append(i)
        if predicted_structure[i] == "(":
            pos_stack_pred.append(i)
        elif predicted_structure[i] == ")":
            ridx = i
            lidx = pos_stack_pred.pop()
            pred_pairs.add((lidx,ridx))
    if len(ref_pseudoknots_idx) > 0:
        print("{} nt corresponds to pseudoknots among {} nt were excluded from calculation".format(len(ref_pseudoknots_idx),L))
    N_predicted = len(pred_pairs)
    N_reference = len(ref_pairs)
    TP = len(pred_pairs.intersection(ref_pairs))
    FN = N_reference - TP
    FP = N_predicted - TP
    #print("TP","FN","FP",sep="\t")
    #print(TP,FN,FP,sep="\t")
    sensitivity = TP/(TP+FN)
    PPV = TP/(TP+FP)
    F = 2*sensitivity*PPV/(sensitivity+PPV)
    MCC = TP/np.sqrt((TP+FP)*(TP+FN))
    return sensitivity,PPV,F,MCC,TP,FN,FP
    

def main():
    parser = argparse.ArgumentParser(description='Evaluate the accuracy of predicted RNA secondary structure')
    parser.add_argument('--reference', '-r',required=True,help="Reference RNA ground true structure, in dot bracket format")
    parser.add_argument('--predicted', '-p',required=True,help="Predicted RNA structure, in dot bracket format")
    parser.add_argument('--performance','-o',required=True,help="Performance of structure prediction")
    args = parser.parse_args()  
    reference_structures = loadRecords(args.reference,"sequence,structure")
    predicted_structures = loadRecords(args.predicted,"sequence,structure")
    
    ref_seq_ids = set(reference_structures.keys())
    pred_seq_ids = set(predicted_structures.keys())
    print("{} records in reference structure".format(len(ref_seq_ids)))
    print("{} records in predicted structure".format(len(pred_seq_ids)))  
    seq_ids = ref_seq_ids.intersection(pred_seq_ids)
    print("{} records in common".format(len(seq_ids))) 
    records = []
    fields = ["seq_id","sensitivity","PPV","F","MCC","TP","FN","FP"]
    for seq_id in sorted(list(seq_ids)):
        valid = True
        reference_structure = reference_structures[seq_id]["structure"]
        predicted_structure = predicted_structures[seq_id]["structure"]
        try:
            checkDBN(reference_structure)
        except:
            valid = False
            print("The reference structure of {} contains invalid base pairs".format(seq_id))
        try:
            checkDBN(predicted_structure)
        except:
            valid = False
            print("The predicted structure of {} contains invalid base pairs".format(seq_id))
        if len(reference_structure) != len(predicted_structure):
            print("The reference structure and predicted structure of {} has different length".format(seq_id))
            valid = False
        if not valid:
            print("Skip {}".format(seq_id))
            continue
        else:
            print("Processing {} ...".format(seq_id))
        sensitivity,PPV,F,MCC,TP,FN,FP = evaluate(predicted_structure,reference_structure)
        records.append((seq_id,sensitivity,PPV,F,MCC,TP,FN,FP))
    performance = pd.DataFrame.from_records(records)
    performance.columns = fields
    performance = performance.set_index("seq_id")
    performance.to_csv(args.performance,sep="\t")
         
if __name__ == "__main__":
    main()
    


    
