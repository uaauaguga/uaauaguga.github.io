#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
import io
from multiprocessing import Pool
import shutil
from utilities import *

def concatdbn(indir,output):
    fout = open(output,"w")
    for dbn in os.listdir(indir):
        path = os.path.join(indir,dbn)
        with open(path) as f:
            for line in f:
                fout.write(line)
    fout.close()

def splitFasta(fasta,splitted_fasta_dir):
    sequences = loadFasta(fasta)
    for name,sequence in sequences.items():
        sequence_path = os.path.join(splitted_fasta_dir,name + ".fa")
        with open(sequence_path, "w") as f:
            f.write(">" + name + "\n")
            f.write(sequence + "\n")
    return list(sequences.keys())

def Fold(fasta_path,dbn_path,constraint_path=None,shape_path=None):
    cmd = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/Fold","--MFE","--bracket"]
    f =  open(dbn_path,"w")
    if constraint_path is not None:
        cmd += ["-c",constraint_path]
    if shape_path is not None:
        cmd += ["--SHAPE",shape_path]
    cmd += [fasta_path,"-"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        if line.startswith(">"):
            line = ">" + line.split(" ")[-1]
        f.write(line)
    f.close()
    code = proc.poll()
    if code != 0:
        print(" ".join(cmd) + " failed")
    else:
        print("Folding {} finished".format(fasta_path))
    return code

def runRNAstructureFold(args):
    if args.fasta is None:
        print("RNAstructure requires fasta file in addition to constraint file.")
        sys.exit(2)

    splitted_fasta_dir = os.path.join(args.tmp_dir,"fasta")
    dbn_dir = os.path.join(args.tmp_dir,"dbn")
    os.mkdir(splitted_fasta_dir)
    os.mkdir(dbn_dir)

    print("Split fasta file ...")
    sequence_ids = splitFasta(args.fasta,splitted_fasta_dir)
    print("Done .")

    if args.shape is not None:
        print("Split shape file ...")
        shapes = loadSHAPE(args.shape)
        os.mkdir(os.path.join(args.tmp_dir,"shape"))
        for seq_id in shapes.keys():
            shape_path = prepareSHAPE(shapes[seq_id],os.path.join(args.tmp_dir,"shape",seq_id))
        print("Done .")

    print("Folding sequence using {} threads ...".format(args.jobs))
    pool = Pool(args.jobs)
    workers = []
    for seq_id in sequence_ids:
        sequence_path = os.path.join(splitted_fasta_dir,seq_id + ".fa")
        structure_path = os.path.join(dbn_dir, seq_id + ".dot") 
        constraint_path = os.path.join(args.constraint, seq_id + ".CON") if args.constraint is not None else None
        shape_path =  os.path.join(args.tmp_dir, "shape", seq_id+ ".shape") if args.shape is not None else None
        if (shape_path is not None) and (not os.path.exists(shape_path)):
            shape_path = None
        workers.append(pool.apply_async(func=Fold, args=(sequence_path,structure_path,constraint_path,shape_path)))
    for worker in workers:
        code = worker.get() 
    outdir = os.path.dirname(args.output)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    concatdbn(dbn_dir,args.output)
    #shutil.rmtree(args.tmp_dir)
    print("Done .")
            

def RNAfold(fasta_path,dbn_path,constraint_path=None,shape_path=None,enforce=True):
    f =  open(dbn_path,"w")
    cmd = ["RNAfold", "--noPS"]
    if constraint_path is not None:
        cmd = ["RNAfold", "--noPS","-C",constraint_path]
        if enforce:
            cmd += ["--enforceConstraint"]
    else:
        cmd += [fasta_path]
    if shape_path is not None:
        cmd += ["--shape="+shape_path,"--shapeMethod=D"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    i = 0
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        line = line.strip().split(" ")[0] + "\n"
        if i == 1:
            line = line.replace("U","T")
        f.write(line)
        i += 1
    f.close()
    code =  proc.poll()
    if code != 0:
        print(" ".join(cmd)+" failed .")
    return code


def runViennaRNAPackageRNAfold(args):
    if args.fasta is None and args.constraint is None:
        print("To run ViennaRNA RNAfold, either constraint file or fasta file should be specified")
        sys.exit(3)

    dbn_dir = os.path.join(args.tmp_dir,"dbn")
    os.mkdir(dbn_dir)

    if args.constraint is None:
        print("Split fasta file ...")
        splitted_fasta_dir = os.path.join(args.tmp_dir,"fasta")
        os.mkdir(splitted_fasta_dir)
        sequence_ids = splitFasta(args.fasta,splitted_fasta_dir)
        print("Done .")
    else:
        sequence_ids = [ file[:file.rfind(".")+1] for file in os.listdir(args.constraint) ]

    if args.shape is not None:
        print("Split shape file ...")
        shapes = loadSHAPE(args.shape)
        os.mkdir(os.path.join(args.tmp_dir,"shape"))
        for seq_id in shapes.keys():
            shape_path = prepareSHAPE(shapes[seq_id],os.path.join(args.tmp_dir,"shape",seq_id))
        print("Done .")

    print("Folding sequence using {} threads ...".format(args.jobs))
    pool = Pool(args.jobs)
    workers = []
    for seq_id in sequence_ids:
        constraint_path = os.path.join(args.constraint,seq_id + ".CON") if args.constraint is not None else None
        fasta_path = None if args.constraint is not None else os.path.join(splitted_fasta_dir,seq_id+".fa")
        dbn_path = os.path.join(dbn_dir,seq_id+".dot")
        shape_path = os.path.join(args.tmp_dir,"shape",seq_id+".shape") if args.shape is not None else None
        if (shape_path is not None) and (not os.path.exists(shape_path)):
            shape_path = None
        workers.append(pool.apply_async(func=RNAfold,args=(fasta_path,dbn_path,constraint_path,shape_path)))
    for i,worker in enumerate(workers):
        code = worker.get()
    outdir = os.path.dirname(args.output)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    concatdbn(dbn_dir,args.output)
    #shutil.rmtree(args.tmp_dir)
    print("Done .")

def main():
    parser = argparse.ArgumentParser(description='Fold input RNA with hard constraints or experimental constraints')
    parser.add_argument('--fasta','-f',help="Input fasta")
    parser.add_argument('--output','-o',required=True,help="Output path for predicted structure in dot bracket notation")
    parser.add_argument('--method','-m', default = "RNAstructure", choices = ["RNAstructure","ViennaRNA"], help = "Which software to use for folding")
    parser.add_argument('--constraint', '-c', type=str, help = "Input directory constrains file of RNAstructure or VisnnaRNA package")
    parser.add_argument('--shape','-s',type=str,help = "The shape reactivity file")
    parser.add_argument('--jobs','-j', type=int, default=4, help = "Threads for processing.")
    parser.add_argument('--tmp-dir','-t',required = True,help = "Directory for store temporary files")
    args = parser.parse_args()
    if os.path.exists(args.tmp_dir):
        print("The specified temporary directory already exists")
        sys.exit(1)
    else:
        os.makedirs(args.tmp_dir)
    if args.method == "RNAstructure":
        runRNAstructureFold(args)
    elif args.method == "ViennaRNA":
        runViennaRNAPackageRNAfold(args)
    
    

if __name__ == "__main__":
    main()



