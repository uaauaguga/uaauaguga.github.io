#!/usr/bin/env python
import io
import pandas as pd
import subprocess
import argparse

def get_lengths(path):
    lengths = {}
    with open(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            L = int(fields[2]) - int(fields[1])
            lengths[fields[0]] = L
    return lengths

def main():
    parser = argparse.ArgumentParser(description='Evaluate the accuracy of identified motif locations')
    parser.add_argument('--input','-i',type=str,required=True,help="Input predicted locations in bed format")
    parser.add_argument('--reference','-r',type=str,required=True,help="Actual motif locations in bed format ")
    parser.add_argument('--site-threshold','-t',type=float,default=0.25,help="If the predicted sites overlaps to reference sites equal to or larger than this value , the predicted site is considered as a real hit")
    parser.add_argument('--overlaps','-o',type=str,required=True,help="Performance of the prediction")
    parser.add_argument('--performance','-p',type=str,help="Output path for performance")
    args = parser.parse_args()
    
    identified_lengths = get_lengths(args.input)
    ref_lengths = get_lengths(args.reference)
    overlap_names = set()
    false_names = set()
    cmd = ["bedtools","intersect","-wo","-a",args.reference,"-b",args.input]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    #records = [] 
    fout = open(args.overlaps,"w")
    l_overlap_total = 0
    l_ref_unique_total = 0
    l_pred_unique_total = 0
    
    print("name","overlap","ref-unique","pred-unique",sep="\t",file=fout)
    # Check regions with overlap
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        #name                               start   end     name                            start   end     overlap 
        #AAAB01044800.1=486-531=59-105	59	105	AAAB01044800.1=486-531=59-105	64	80	16
        fields = line.strip().split("\t")
        name = fields[0]
        l_overlap = int(fields[6])
        l_overlap_total += l_overlap
        l_ref_unique = int(fields[2]) - int(fields[1]) - l_overlap
        l_ref_unique_total += l_ref_unique
        l_pred_unique = int(fields[5]) - int(fields[4]) - l_overlap
        l_pred_unique_total += l_pred_unique
        #records.append((name,l_overlap,l_ref_unique,l_pred_unique))
        print(name,l_overlap,l_ref_unique,l_pred_unique,sep="\t",file=fout)
        if l_overlap/(l_overlap + l_ref_unique) >= args.site_threshold:    
            overlap_names.add(name)
        else:
            false_names.add(name)
    false_names = false_names.union(set(identified_lengths.keys()).difference(overlap_names))
    missed_names = set(ref_lengths.keys()).difference(overlap_names)
    
    # Check regions only present in predictions, or false predictions:
    for name in false_names:
        l_pred_unique_total += identified_lengths[name]
        print(name,0,0,identified_lengths[name],sep="\t",file=fout)
    
    # Check regions only present in references, or missed motif:
    for name in missed_names:
        l_ref_unique_total += ref_lengths[name]
        print(name,0,ref_lengths[name],0,sep="\t",file=fout)

    n_hits,n_miss,n_false = len(overlap_names),len(missed_names),len(false_names)

    print("Performance at site level: {} hits, {} were missed, {} were false prediction".format(n_hits,n_miss,n_false))
    print("Performance at position level: {} nt overlap, {} nt uniqe to reference, {} nt uniqe to prediction".format(l_overlap_total,l_ref_unique_total,l_pred_unique_total))

    performance = pd.DataFrame(index=["site","position"],columns=["TP","FN","FP"])
    performance.loc["site",:] = [n_hits,n_miss,n_false]
    performance.loc["position",:] = [l_overlap_total,l_ref_unique_total,l_pred_unique_total]
    print(performance)
    if args.performance is not None:
        performance.to_csv(args.performance,sep="\t")
    fout.close()


if __name__ == "__main__":
    main()


