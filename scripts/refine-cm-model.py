#!/usr/bin/env python
import os
import subprocess
import argparse
import logging
from Bio import AlignIO
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("cm model refinement")

def main():
    parser = argparse.ArgumentParser(description="refine cm model by iterative search within a genome")
    parser.add_argument('--model','-cm',default="data/RNIE.seed-a-14.cm",help="Input cm model")
    parser.add_argument('--genome','-g', help="input bacteria genome", required=True)
    parser.add_argument('--outdir','-o', help="output directory",required=True)
    parser.add_argument('--max-iterations','-n',type=int,default=5,help="number of searching iterations")
    parser.add_argument('--min-improvement','-m',type=int,default=5,help="stop if the last round of lead to mild improvement in hit number")
    parser.add_argument('--e-value','-e',type=float,default=0.1,help="E values cutoff for alignments to build cm model")
    parser.add_argument('--nohmm',action="store_true",help="whether use hmm filter in cmsearch")
    args = parser.parse_args()    
    if not os.path.exists(args.outdir):
        logger.info(f"output directory {args.outdir} does not exists, create it .")
        os.mkdir(args.outdir)
    cmmodel_path = args.model
    last_hits_number = 0
    n_sequence = 0
    with open(args.genome) as f:
        for line in f:
            if line.startswith(">"):
                n_sequence += 1
    logger.info(f"{n_sequence} sequence present in input sequence .")
    devnull = open(os.devnull,"w")
    for n in range(1,args.max_iterations+1):
        logger.info(f"start round {n} of iterative search ...")
        msa_path = f"{args.outdir}/round-{n}.stk"
        tbl_path = f"{args.outdir}/round-{n}.tbl"
        cmd = ["cmsearch","--incE",str(args.e_value),"-A",msa_path,"--tblout",tbl_path,"--noali"]
        if args.nohmm:
            cmd += ["--nohmm"]
        cmd += [cmmodel_path, args.genome]
        subprocess.run(cmd,stdout=devnull)
        gff_path = f"{args.outdir}/round-{n}.gff"
        subprocess.run(["scripts/rfam-tbl2gff.py","-i",tbl_path,"-o",gff_path])
        bed_path = f"{args.outdir}/round-{n}.bed"
        subprocess.run(["scripts/gff2bed.py","--gff",gff_path,"--bed",bed_path,"--feature","ncRNA","--name","ID"]) 
        logger.info(f"updating cm model with MSA {msa_path} ...")
        MSA = AlignIO.read(msa_path,"stockholm") 
        hits_number = len(MSA)
        seq_with_RIT = set()
        for record in MSA:
            seq_id = "/".join(record.id.split("/")[:-1])
            seq_with_RIT.add(seq_id)
        improvement = hits_number-last_hits_number
        logger.info(f"{hits_number} hits detected with E-value < {args.e_value}, increase by {improvement} .")
        logger.info(f"rebuild cm model ...")
        cmmodel_path = f"{args.outdir}/round-{n}.cm"
        cmd = ["cmbuild","-F","--refine",f"{args.outdir}/round-{n}-refined.stk","-n",args.genome.split("/")[-1] + f"-round-{n}",cmmodel_path, msa_path]
        subprocess.run(cmd,stdout=devnull)
        logger.info(f"calibrate cm model {msa_path} ...")
        cmd = ["cmcalibrate",cmmodel_path]
        subprocess.run(cmd,stdout=devnull)
        logger.info(f"Finishing search round {n}")
        if improvement < args.min_improvement and n != args.max_iterations:
            logger.info(f"Seems further searching lead to little improvement, stopping. ")
            break
        logger.info(f"Identify RIT within {len(seq_with_RIT)} of the {n_sequence} input sequences .")
        last_hits_number = hits_number
    logger.info("All done .")
    devnull.close()
    with open(f"{args.outdir}/final-hits.txt","w") as f:
        f.write(f"{len(seq_with_RIT)}/{n_sequence}")

if __name__ == "__main__":
    main()    
