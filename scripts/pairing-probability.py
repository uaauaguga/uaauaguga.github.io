#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
import io
import shutil
from multiprocessing import Pool
import shutil
import numpy as np
from utilities import *

def RNAstructure(sequence,prefix,shape=None):
    outdir,basename = os.path.split(prefix)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    partition_file_path = prefix + ".pfs"
    pairing_proba_mat_path = prefix + ".txt"
    if shape is None:
        command = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/partition","-",partition_file_path]
    else:
        shape_path = prepareSHAPE(shape,prefix)
        command = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/partition","-","--SHAPE",shape_path,partition_file_path]
    p = subprocess.Popen(command, stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    _ = p.communicate(input=sequence.encode())
    _ = subprocess.run(["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/ProbabilityPlot",partition_file_path,pairing_proba_mat_path,"--text"],stdout=subprocess.PIPE)
    epsilon = 1e-14
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
    pairing_prob = pairing_prob_matrix.sum(axis=1)
    os.remove(partition_file_path)
    os.remove(pairing_proba_mat_path)
    return pairing_prob

def ViennaRNA(sequence,prefix,shape=None):
    outdir,name = os.path.split(prefix)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    input = ">{}\n{}".format(name,sequence)
    if shape is None:
        command = ["RNAplfold","--ulength=1","--winsize=200"]
    else:
        shape_path = prepareSHAPE(shape,prefix)
        command = ["RNAplfold","--ulength=1","--winsize=200","--shape",shape_path,"--shapeMethod=D"]
    p = subprocess.Popen(command,stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    _ = p.communicate(input=input.encode())
    ps_path = "{}_dp.ps".format(name.replace(":","_"))
    un_prob_path = "{}_lunp".format(name.replace(":","_"))
    un_prob = []
    with open(un_prob_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            un_prob.append(float(line.split("\t")[1]))
    un_prob = np.array(un_prob)
    os.rename(ps_path,os.path.join(outdir,ps_path))
    os.rename(un_prob_path,os.path.join(outdir,un_prob_path))
    return 1 - un_prob
    

def main():
    parser = argparse.ArgumentParser(description='Fold input RNA with hard constraints or experimental constraints')
    parser.add_argument('--fasta','-f',required=True,help="Input fasta")
    parser.add_argument('--output','-o',required=True,help="Output path for predicted structure in dot bracket notation")
    parser.add_argument('--method','-m', default = "RNAstructure", choices = ["RNAstructure","ViennaRNA"], help = "Which software to use for folding")
    parser.add_argument('--shape','-sh',help="Path of SHAPE constraint,similar to fasta format, the sequence id should be consistent with fasta input")
    parser.add_argument('--tmp-dir','-t',default="tmp",help = "Directory for store temporary files")
    parser.add_argument('--keep-tmp','-kt',action="store_true",default=False,help="Whether keep the temporary directory")
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
            pairing_prob = RNAstructure(sequence,args.tmp_dir+"/"+seq_id,shape=shape)
        elif args.method == "ViennaRNA":
            pairing_prob = ViennaRNA(sequence,args.tmp_dir+"/"+seq_id,shape=shape)
        fout.write(">"+seq_id+"\n")
        fout.write(",".join(np.round(pairing_prob,3).astype(str))+"\n")
    fout.close()
    if not args.keep_tmp:
        print("Removing temporary files ...")
        shutil.rmtree(args.tmp_dir)
        print("Done .")
        
    


if __name__ == "__main__":
    main()
