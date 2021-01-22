#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
import io
from multiprocessing import Pool
import shutil


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

def concatdbn(indir,output):
    fout = open(output,"w")
    for dbn in os.listdir(indir):
        path = os.path.join(indir,dbn)
        with open(path) as f:
            for line in f:
                fout.write(line)
    fout.close()


def Fold(fasta_path,dbn_path,constraint_path=None):
    cmd = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py27/bin/Fold","--MFE","--bracket"]
    f =  open(dbn_path,"w")
    if constraint_path is not None:
        cmd += ["-c",constraint_path]
    cmd += [fasta_path,"-"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        if line.startswith(">"):
            line = ">" + line.split(" ")[-1]
        f.write(line)
    f.close()
    code = proc.poll()
    if code != 0:
        print(" ".join(cmd))
    return code



    

def runRNAstructureFold(args):
    splitted_fasta_dir = os.path.join(args.tmp_dir,"fasta")
    dbn_dir = os.path.join(args.tmp_dir,"dbn")
    os.mkdir(splitted_fasta_dir)
    os.mkdir(dbn_dir)

    print("Split fasta file ...")
    sequences = loadFasta(args.fasta)
    for name,sequence in sequences.items():
        sequence_path = os.path.join(splitted_fasta_dir,name + ".fa")
        with open(sequence_path, "w") as f:
            f.write(">" + name + "\n")
            f.write(sequence + "\n")
    print("Done .")

    if args.constraint is not None:
        print("Constraint provided .")

    print("Folding sequence using {} threads ...".format(args.jobs))
    pool = Pool(args.jobs)
    workers = []
    for name in sequences.keys():
        sequence_path = os.path.join(splitted_fasta_dir,name + ".fa")
        structure_path = os.path.join(dbn_dir,name + ".dot") 
        if args.constraint is not None:
            constraint_path = os.path.join(args.constraint,name + ".CON")
        else:
            constraint_path = None
        workers.append(pool.apply_async(func=Fold, args=(sequence_path,structure_path,constraint_path)))
    for worker in workers:
        code = worker.get() 
    concatdbn(dbn_dir,args.output)
    shutil.rmtree(args.tmp_dir)
    print("Done .")
            

def RNAfold(constraint_path,dbn_path,enforce=True):
    f =  open(dbn_path,"w")
    cmd = ["RNAfold", "--noPS","-C",constraint_path]
    if enforce:
        cmd += ["--enforceConstraint"]
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
        print(" ".join(cmd))
    return code




def runViennaRNAPackageRNAfold(args):
    dbn_dir = os.path.join(args.tmp_dir,"dbn")
    os.mkdir(dbn_dir)
    print("Folding sequence using {} threads ...".format(args.jobs))
    if args.constraint is not None:
        print("Constraint provided .")
        pool = Pool(args.jobs)
        workers = []
        names = []
        for constraint_file in os.listdir(args.constraint):
            constraint_path = os.path.join(args.constraint,constraint_file)
            name = ".".join(constraint_file.split(".")[:-1])
            dbn_path = os.path.join(dbn_dir,name+".dot")
            workers.append(pool.apply_async(func=RNAfold,args=(constraint_path,dbn_path,args.enforce)))
            names.append(name)
        for i,worker in enumerate(workers):
            code = worker.get()
        concatdbn(dbn_dir,args.output)
        shutil.rmtree(args.tmp_dir)
    else:
        print("Constraint not provided .")
        cmd = ["RNAfold", "--noPS", "--jobs="+str(args.jobs), "-o", args.output, args.fasta]
        subprocess.run(cmd)
    print("Done .")

def main():
    parser = argparse.ArgumentParser(description='Fold input RNA with hard constraints or experimental constraints')
    parser.add_argument('--fasta','-f',help="Input fasta")
    parser.add_argument('--output','-o',required=True,help="Output path for predicted structure in dot bracket notation")
    parser.add_argument('--method','-m', default = "RNAstructure", choices = ["RNAstructure","ViennaRNA"], help = "Which software to use for folding")
    parser.add_argument('--constraint', '-c', type=str, help = "Input constrains file of RNAstructure or VisnnaRNA package")
    parser.add_argument('--enforce', '-e', action = "store_true", default=False, help = "For ViennaRNA RNAfold, enforceConstraints is a flag indicating whether or not constraints for base pairs should be enforced instead of just doing a removal of base pair that conflict with the constraint")
    parser.add_argument('--reactivity','-r',type=str,help = "The experiment constraint file")
    parser.add_argument('--jobs','-j', type=int, default=4, help = "Threads for processing.")
    parser.add_argument('--tmp-dir','-t',help = "Directory for store temporary files")
    args = parser.parse_args()
    if args.method == "RNAstructure":
        if args.fasta is None:
            print("RNAstructure requires fasta file in addition to constraint file.")
            sys.exit(1)
        if args.tmp_dir is None:
            print("A tmp dir should be specified .")
            sys.exit(2)
        if not os.path.exists(args.tmp_dir):
            os.mkdir(args.tmp_dir)
        else:
            print("The specified dir already exists .")
            sys.exit(3)

        runRNAstructureFold(args)
    elif args.method == "ViennaRNA":
        if args.fasta is None and args.constraint is None:
            print("To run ViennaRNA RNAfold, either constraint file or fasta file should be specified")
            sys.exit(4)
        if args.constraint is not None:
            if args.tmp_dir is None:
                print("A tmp dir should be specified .")
                sys.exit(2)
            if not os.path.exists(args.tmp_dir):
                print(args.tmp_dir)
                os.mkdir(args.tmp_dir)
            else:
                print("The specified dir already exists .")
                sys.exit(3)
        runViennaRNAPackageRNAfold(args)
    
    

if __name__ == "__main__":
    main()



