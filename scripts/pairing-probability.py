#!/usr/bin/env python
import argparse
import subprocess
import os
import tempfile
import sys
import io
import pandas as pd
import shutil
from multiprocessing import Pool
import shutil
import numpy as np
import re
from utilities import *

def RNAstructure(sequence,prefix,shape=None,window=200):
    outdir,basename = os.path.split(prefix)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    partition_file_path = prefix + ".pfs"
    pairing_proba_mat_path = prefix + ".txt"
    if shape is None:
        command = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/partition","-","--maxdistance",str(window),partition_file_path]
    else:
        shape_path = prepareSHAPE(shape,prefix)
        command = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/partition","-","--maxdistance",str(window),"--SHAPE",shape_path,partition_file_path]
        print(" ".join(command))
    p = subprocess.Popen(command, stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    _ = p.communicate(input=sequence.encode())
    _ = subprocess.run(["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/ProbabilityPlot",partition_file_path,pairing_proba_mat_path,"--text"],stdout=subprocess.PIPE)
    with open(pairing_proba_mat_path) as f:
        length = next(f)
        length = int(length.strip())
        pairing_matrix = np.ones((length,length))*(np.inf)
        _ = next(f)
        for line in f:
            line = line.strip()
            i,j,log10p = line.strip().split("\t")
            i,j,log10p = int(i)-1,int(j)-1,float(log10p)
            pairing_matrix[i,j] = pairing_matrix[j,i] = log10p
    pairing_matrix = np.log(10)*pairing_matrix
    pairing_prob_matrix = np.exp(-pairing_matrix)
    pairing_prob_matrix += np.diag(1-pairing_prob_matrix.sum(axis=1))
    #os.remove(partition_file_path)
    #os.remove(pairing_proba_mat_path)
    return pairing_prob_matrix


def RNAstructure(sequence,prefix,shape=None,window=200):
    outdir,basename = os.path.split(prefix)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    partition_file_path = prefix + ".pfs"
    pairing_proba_mat_path = prefix + ".txt"
    if shape is None:
        command = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/partition","-","--maxdistance",str(window),partition_file_path]
    else:
        shape_path = prepareSHAPE(shape,prefix)
        command = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/partition","-","--maxdistance",str(window),"--SHAPE",shape_path,partition_file_path]
        print(" ".join(command))
    p = subprocess.Popen(command, stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    _ = p.communicate(input=sequence.encode())
    _ = subprocess.run(["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/ProbabilityPlot",partition_file_path,pairing_proba_mat_path,"--text"],stdout=subprocess.PIPE)
    with open(pairing_proba_mat_path) as f:
        length = next(f)
        length = int(length.strip())
        pairing_matrix = np.ones((length,length))*(np.inf)
        _ = next(f)
        for line in f:
            line = line.strip()
            i,j,log10p = line.strip().split("\t")
            i,j,log10p = int(i)-1,int(j)-1,float(log10p)
            pairing_matrix[i,j] = pairing_matrix[j,i] = log10p
    pairing_matrix = np.log(10)*pairing_matrix
    pairing_prob_matrix = np.exp(-pairing_matrix)
    pairing_prob_matrix += np.diag(1-pairing_prob_matrix.sum(axis=1))
    #os.remove(partition_file_path)
    #os.remove(pairing_proba_mat_path)
    return pairing_prob_matrix

def ViennaRNA(sequence,prefix,shape=None,window=200):
    outdir,name = os.path.split(prefix)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    name = name + "_" + next(tempfile._get_candidate_names())
    length = len(sequence)
    input = ">{}\n{}".format(name,sequence)
    if shape is None:
        command = ["RNAplfold","--cutoff=0.0001","--print_onthefly","--ulength=1","--winsize="+str(window),"--span="+str(window)]
    else:
        shape_path = prepareSHAPE(shape,prefix)
        command = ["RNAplfold","--cutoff=0.0001","--print_onthefly","--ulength=1","--winsize="+str(window),"--shape",shape_path,"--shapeMethod=D","--span="+str(window)]
    p = subprocess.Popen(command,stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    _ = p.communicate(input=input.encode())
    bp_prob_path = "{}_basepairs".format(name.replace(":","_"))
    un_prob_path = "{}_lunp".format(name.replace(":","_"))
    pairing_prob_matrix = np.zeros((length,length))

    with open(bp_prob_path) as f:
        for line in f:
            line = line.strip()
            i,j,v = re.split("\s+",line)
            i,j,v = int(i)-1,int(j)-1,float(v)
            pairing_prob_matrix[i,j] = pairing_prob_matrix[j,i] = v

    with open(un_prob_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            i,v = re.split("\s+",line)
            i,v = int(i)-1,float(v)
            pairing_prob_matrix[i,i] = v
    x = pairing_prob_matrix.sum(axis=1)
    pairing_prob_matrix = pairing_prob_matrix/pairing_prob_matrix.sum(axis=0).reshape(-1,1) 
    os.rename(bp_prob_path,prefix+".bpp.txt")
    os.rename(un_prob_path,prefix+".upp.txt")
    return pairing_prob_matrix


def getLocalBPPM(bppm,window=200):
    local_bppm = np.zeros((window*2+1,bppm.shape[0]))
    for k in np.arange(-window,window+1):
        flanking = np.zeros(np.abs(k))
        if k < 0:
            # Probability of pairing with upstream bases
            x = np.concatenate([flanking,np.diag(bppm,-k)])
        elif k == 0:
            # Probability of unpaired
            x = np.concatenate([np.diag(bppm)])
        else:
            # Probability of pairing with downstream bases
            x = np.concatenate([np.diag(bppm,-k),flanking])
        local_bppm[k+window,:] = x
    return local_bppm
     

def getEntropy(bppm):
    bppm += 1e-4
    bppm = bppm/bppm.sum(axis=0).reshape((1,-1))
    return -(bppm*np.log(bppm)).sum(axis=0) 


def main():
    parser = argparse.ArgumentParser(description='Fold input RNA with hard constraints or experimental constraints')
    parser.add_argument('--fasta','-f',required=True,help="Input fasta")
    parser.add_argument('--output','-o',required=True,help="Output path for predicted structure in dot bracket notation")
    parser.add_argument('--method','-m', default = "RNAstructure", choices = ["RNAstructure","ViennaRNA"], help = "Which software to use for folding")
    parser.add_argument('--entropy','-e',action="store_true",default=False,help="Whether calculate entropy of pairing probability")
    parser.add_argument('--window','-w',type=int,default=200,help="Max length that a base pair can spanning")
    parser.add_argument('--shape','-sh',help="Path of SHAPE constraint,similar to fasta format, the sequence id should be consistent with fasta input")
    parser.add_argument('--tmp-dir','-t',default="tmp",help = "Directory for store temporary files")
    parser.add_argument('--keep-tmp','-kt',action="store_true",default=False,help="Whether keep the temporary directory")
    parser.add_argument('--jobs','-j',default=1,help="Number of process used")
    args = parser.parse_args()
    fasta = loadFasta(args.fasta)
    if args.shape is not None:
        shapes = loadSHAPE(args.shape)
    else:
        shape = None
    fout = open(args.output,"w")
    for seq_id, sequence in fasta.items():
        print(seq_id)
        if args.shape is not None:
            if seq_id not in shapes:
                shape = None
            else:  
                try:
                    assert len(sequence) == shapes[seq_id].shape[0]
                    shape = shapes[seq_id]
                except:
                    print("The length of {} and its SHAPE reacitvities are different".format(seq_id))
                    continue
        if args.method == "RNAstructure":
            bppm = RNAstructure(sequence,args.tmp_dir+"/"+seq_id,shape=shape,window=args.window)
        elif args.method == "ViennaRNA":
            bppm = ViennaRNA(sequence,args.tmp_dir+"/"+seq_id,shape=shape,window=args.window)        
        pairing_prob = 1 - np.diag(bppm)
        fout.write(">"+seq_id+"\n")
        fout.write(",".join(np.round(pairing_prob,3).astype(str))+"\n")
        local_bppm = getLocalBPPM(bppm,window=args.window)
        if args.entropy:
            entropy = getEntropy(local_bppm)
            fout.write(",".join(np.round(entropy,3).astype(str))+"\n")
        bppm = pd.DataFrame(bppm).round(4)
        bppm.to_csv(os.path.join(args.tmp_dir,seq_id+".bppm.txt"),header=None,sep="\t")
        local_bppm = pd.DataFrame(local_bppm,index=np.arange(-args.window,args.window+1)).round(4)
        local_bppm.to_csv(os.path.join(args.tmp_dir,seq_id+".local.bppm.txt"),header=None,sep="\t")
    fout.close()
    if not args.keep_tmp:
        print("Removing temporary files ...")
        shutil.rmtree(args.tmp_dir)
        print("Done .")
        
    


if __name__ == "__main__":
    main()
