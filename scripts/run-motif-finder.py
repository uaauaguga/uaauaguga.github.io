#!/usr/bin/env python
import subprocess
import argparse
import os
import sys
import logging
import re
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('Run motif finder')



def prepareFasta(args,U2T=False,T2U=False):
    fasta_path = os.path.join(args.outdir,"sequences.fa")
    ffasta = open(fasta_path,"w")
    for dbfile in os.listdir(args.indir):
        with open(os.path.join(args.indir,dbfile)) as f:
            for i,line in enumerate(f):
                if i == 0:
                    line = line.replace(":","=")
                    ffasta.write(line)
                elif i == 1:
                    if U2T:
                        line = line.replace("U","T")
                    elif T2U:
                        line = line.replace("T","U")
                    #line = re.sub(r"[NRKSWY]","A",line) currently used for MEMERIS as a dirty fix
                    ffasta.write(line)
                else:
                    continue
    ffasta.close()
    return os.path.join(os.getcwd(),fasta_path)

def BEAM(args):
    dbn_path = os.path.join(args.outdir,"sequences.dot")
    fdbn = open(dbn_path,"w")    

    if args.use_reference_structure:
        logger.info("Reference structure provided .")
        for dbfile in os.listdir(args.indir):
            with open(os.path.join(args.indir,dbfile)) as f:
                for i,line in enumerate(f):
                    if i == 0:
                        line = line.replace(":","=")
                        fdbn.write(line.split(" ")[0]+"\n")
                    elif i == 2:
                        fdbn.write(line.split(" ")[0]+"\n")
                    else:
                        fdbn.write(line.replace("T","U"))
    else:
        logger.info("Reference structure not provided, fold use RNAfold")
        fastaPath = prepareFasta(args)
        foldedPath = os.path.join(args.outdir,"sequences.folded.dot")
        with open(foldedPath,"w") as f:
            subprocess.run(["RNAfold","--noPS",fastaPath],stdout=f)
        with open(foldedPath) as f:
            for i,line in enumerate(f):
                if i%3 == 2:
                    fdbn.write(line.split(" ")[0]+"\n")
                else:
                    fdbn.write(line)
        os.remove(foldedPath)
    fdbn.close()
       
    bearEncoderPath = "/BioII/lulab_b/jinyunfan/software/beam/BearEncoder.new.jar"
    beamPath = "/BioII/lulab_b/jinyunfan/software/beam/BEAM_release2.5.0.jar"
    
    logger.info("Convert input dot bracket file to BEAR encoding")
    bear_path = os.path.join(args.outdir,"sequences.fb")
    subprocess.run(["java","-jar",bearEncoderPath,dbn_path,bear_path])
    logger.info("Done .")
    
    logger.info("Run motif finder ...")
    cwd = os.getcwd()
    os.chdir(args.outdir)
    bear_path = "sequences.fb"
    cmds = ["java","-jar",beamPath,"-f",bear_path]
    if args.min_length is not None:
        cmds += ["-w",str(args.min_length).split(".")[0]] 
    if args.max_length is not None:
        cmds += ["-W",str(args.max_length).split(".")[0]]
    if args.number is not None:
        cmds += ["-M",str(args.number)]
    logger.info("Running "+" ".join(cmds))
    subprocess.run(cmds)
    
    os.chdir(cwd)
    logger.info("Done .")


def RNApromo(args):
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    logger.info("Prepare input sequences ...")
    fasta_path = prepareFasta(args)
    logger.info("Done .")
    logger.info("Run RNApromo ...")
    RNApromoPath = "/BioII/lulab_b/jinyunfan/software/RNApromo/rnamotifs08_motif_finder.pl"
    cmds = [RNApromoPath,"-positive_seq",fasta_path,"-output_dir",os.path.join(args.outdir,"results")]
    if args.use_reference_structure:
        refStructurePath = os.path.join(args.outdir,"ref_structure.tab")
        refStructuref = open(refStructurePath,"w")
        logger.info("Reference structure provided, use provided structure for motif finding.")
        for dbfile in os.listdir(args.indir):
            with open(os.path.join(args.indir,dbfile)) as f:
                for i,line in enumerate(f):
                    line = line.strip()
                    if i%3 == 0:
                        name = line.replace(">","").replace(":","=")
                    elif i%3 == 2:
                        structure = line.split(" ")[0]
                        print(name,"0",structure,sep="\t",file=refStructuref)
        refStructuref.close()
        cmds += ["-positive_struct",refStructurePath]
    else:
        logger.info("Reference structure not provided, use RNApromo internal structure prediction tool.")
    if args.min_length is not None:
        cmds += ["-min",str(args.min_length)]
    if args.max_length is not None:
        cmds += ["-max",str(args.max_length)]
    if args.number is not None:
        cmds += ["-n",str(args.number)]
    logger.info(" ".join(cmds))
    subprocess.run(cmds)
    logger.info("Done .")

def CMfinder(args):
    logger.info("Prepare input sequences ...")
    fasta_path = prepareFasta(args)
    logger.info("Done .")
    logger.info("Run CMfinder ...")
    #cmfinderPath = "/BioII/lulab_b/jinyunfan/software/cmfinder-0.4.1.18/bin/cmfinder.pl" 
    cmfinderPath = "/BioII/lulab_b/jinyunfan/software/cmfinder-0.4.1.18/bin/cmfinder04.pl"
    cmds = [cmfinderPath]
    #if args.min_length is not None:
    #    cmds += ["-m",str(args.min_length).split(".")[0]]
    #if args.max_length is  not None:
    #    cmds += ["-M",str(args.max_length).split(".")[0]]
    cmds += ["-fragmentary","-combine"]
    cmds += [fasta_path]
    logger.info(" ".join(cmds))
    subprocess.run(cmds)
    logger.info("Done .")


def GraphProt(args):
    logger.info("Prepare positive sequences ...")
    posFastaPath = prepareFasta(args)
    negFastaPath = os.path.join(args.outdir,"sequences-dinuc-shuffled.fa")
    logger.info("Done.")
    logger.info("Prepare dimucleotide shuffled sequence as negative set")
    dinucShuffleScript = "/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/fasta-dinucleotide-shuffle-py3"
    with open(negFastaPath,"w") as negf:
        subprocess.run([dinucShuffleScript,"-f",posFastaPath],stdout=negf)
    logger.info("Done.")
    if not os.path.exists(os.path.join(args.outdir,"results")):
        os.mkdir(os.path.join(args.outdir,"results"))
    prefix = os.path.join(args.outdir,"results","GraphProt")
    logger.info("Training GraphProt model ...")
    #print(" ".join(["GraphProt.pl","-action","train","-fasta",posFastaPath,"-negfasta",negFastaPath,"-prefix",prefix]))
    subprocess.run(["GraphProt.pl","-action","train","-fasta",posFastaPath,"-negfasta",negFastaPath,"-prefix",prefix])
    logger.info("Done .")
    
    modelPath = os.path.join(args.outdir,"results","GraphProt.model")
    logger.info("Predict affinity profile...")
    subprocess.run(["GraphProt.pl","-action","predict_profile","-model",modelPath,"-fasta",posFastaPath,"-prefix",prefix]) 
    logger.info("Done .")
    logger.info("Predict high affinity sites ...")
    subprocess.run(["GraphProt.pl","-action","predict_has","-model",modelPath,"-fasta",posFastaPath,"-prefix",prefix,"-percentile","90"])
    
    maxMergeDistance = 3
    records = []
    record = None
    with open(os.path.join(args.outdir,"results","GraphProt.has")) as f:
        lastSeqIdx,lastPos = -100,-100
        for line in f:
            seqIdx,pos,_ = line.strip().split("\t")
            seqIdx,pos = int(seqIdx),int(pos)
            if seqIdx != lastSeqIdx:
                if record is not None:
                    records.append(record)
                record = [seqIdx,pos,pos+1]
            else:
                flag = False
                for gap in range(1,maxMergeDistance+1):
                    if pos == lastPos+gap:
                        record[2] += gap
                        pos = lastPos+gap
                        break
                else:
                    records.append(record)
                    record = [seqIdx,pos,pos+1]
            lastSeqIdx = seqIdx
            lastPos = pos
    with open(posFastaPath) as f:
        seqNames = [line.strip().replace(">","") for line in f if line.startswith(">")]
    floc = open(os.path.join(args.outdir,"results","locations.txt"),"w")
    for record in records:
        if record[2]-record[1] > args.min_length:
            print(seqNames[record[0]],record[1],record[2],file=floc,sep="\t")
    floc.close()
    logger.info("Done .")
    logger.info("Extract motif from trained model ...")
    cmds = ["GraphProt.pl","-action","motif","-model",modelPath,"-fasta",posFastaPath,"-prefix",prefix]
    logger.info(" ".join(cmds))
    if args.length is not None:
        cmds += ["-motif_len",str(args.length)]
    if args.number is not None:
        cmds += ["-motif_top_n",str(args.number)]
    subprocess.run(cmds)
    logger.info("Done .")
      

def MEMERIS(args):
    logger.info("Prepare input sequences ...")
    fastaPath = prepareFasta(args,U2T=True)
    logger.info("Done .")
    #subprocess.run(["export","MEME_DIRECTORY=/BioII/lulab_b/jinyunfan/software/memeris_1.0/build"])
    memerisPath="/BioII/lulab_b/jinyunfan/software/memeris_1.0/build/bin/memeris"
    output = os.path.join(args.outdir,"results.txt")
    cmds = [memerisPath,'-dna']
    if args.min_length is not None:
        cmds += ["-minw",str(args.min_length)]
    if args.max_length is not None:
        cmds += ["-maxw",str(args.max_length)]
    if args.number is not None:
        cmds += ["-nmotifs",str(args.number)] 
    cmds += [fastaPath]
    #print(cmds)
    with open(output,"w") as f:
        subprocess.run(cmds,stdout=f)

def ssHMM(args):
    from sshmm.structure_prediction import calculate_rna_shapes_from_file, calculate_rna_structures_from_file
    logger.info("Prepare input sequence ...") 
    fastaPath = prepareFasta(args,U2T=True)
    logger.info("Done .")
    #train_seqstructhmm --motif_length 40 --output_directory RF00032-40 RF00032.fa  RF00032.txt
    cmds = ["train_seqstructhmm"]
    if args.length is not None:
        cmds  += ["--motif_length",str(args.length)] 
    cmds += ["--job_name","sequences","--output_directory",args.outdir] 
    structurePath = os.path.join(args.outdir,"structure.txt")
    calculate_rna_shapes_from_file(structurePath,fastaPath, 10)
    cmds += [fastaPath,structurePath]
    subprocess.call(cmds)
   

def RNAcontext(args):
    logger.info("Prepare input sequences ...")
    negFastaPath = os.path.join(args.outdir,"sequences-dinuc-shuffled.fa")
    posFastaPath = prepareFasta(args)
    logger.info("Done.")
    logger.info("Prepare dimucleotide shuffled sequence as negative set")
    dinucShuffleScript = "/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/fasta-dinucleotide-shuffle-py3"
    with open(negFastaPath,"w") as negf:
        subprocess.run([dinucShuffleScript,"-f",posFastaPath],stdout=negf)
    logger.info("Done.")
    trainingSetPath = os.path.join(args.outdir,"training_sequences.txt")
    
    logger.info("Prepare RNAcontext input ...")
    
    def label2seq(fin,fout,label):
        window_size = 200
        for line in fin:
            if line.startswith(">"):
                continue
            line = line.strip().replace("T","U")
            if len(line)<200:
                window_size = len(line) - 1
            fout.write(label+" "+line.strip()+"\n")
        return window_size
    with open(trainingSetPath,"w") as fout:
        with open(posFastaPath) as posfin:
            window_size = label2seq(posfin,fout,"1")
        with open(negFastaPath) as negfin:
            _ = label2seq(negfin,fout,"-1")
    combinedFastaPath = os.path.join(args.outdir,"sequences-combined.fa")
    with open(combinedFastaPath,"w") as fcombined:
        subprocess.run(["cat",posFastaPath,negFastaPath],stdout=fcombined)
    max_dist = min([160,window_size])
    profilePaths = []
    for profileType in ["E","H","I","M"]:
        logger.info("Prepare profile for {}".format(profileType))
        plfoldExePath = "scripts/RNAplfold_scripts/{}_RNAplfold".format(profileType)
        profilePath = os.path.join(args.outdir,"{}_profile.txt".format(profileType))
        profilePaths.append(profilePath)
        outProfilef = open(profilePath,"w")
        inFastaf = open(combinedFastaPath,"r")
        cmds = [plfoldExePath,"-W",str(window_size),"-L",str(max_dist),"-u","1"]
        subprocess.run(cmds,stdin=inFastaf,stdout=outProfilef)
        outProfilef.close()
        inFastaf.close()
    combinedProfilePath = os.path.join(args.outdir,"combined_profile.txt")
    cmds = ["scripts/RNAplfold_scripts/combine_letter_profiles.py"] + profilePaths + ["1"] + [combinedProfilePath]
    
    logger.info("Combine profiles ...")
    subprocess.run(cmds)
    logger.info("Done .")
    
    logger.info("Run RNAcontext motif finder...")
    cwd = os.getcwd()
    os.chdir(args.outdir)
    RNAcontextExe = "/BioII/lulab_b/jinyunfan/software/RNAcontext/bin/rnacontext"
    cmds = [RNAcontextExe,"-a","ACGU","-e","PHIME","-c","training_sequences.txt","-h","combined_profile.txt","-o","results","-s","5"]
    if args.min_length is not None and args.max_length is not None:
        cmds += ["-w","{}-{}".format(args.min_length,args.max_length)]
    else:
        cmds += ["-w","4-10"]
    if not os.path.exists("outputs"):
        os.mkdir("outputs")
    logger.info("Exute: "+" ".join(cmds))
    subprocess.run(cmds)
    os.chdir(cwd)
    logger.info("Done .")


def zagros(args):
    #/BioII/lulab_b/jinyunfan/software/zagros/zagros-1.0.0/bin/zagros
    if args.length is not None:
        if args.length<4 or args.length>12:
            logger.info("zagros only accept motif width between 4 and 12")
            sys.exit(0)
    logger.info("Prepare input sequences ...")
    fastaPath = prepareFasta(args,U2T=True)
    logger.info("Done .")

    zagrosExePath = "/BioII/lulab_b/jinyunfan/software/zagros/zagros-1.0.0/bin/zagros"
    thermalExePath = "/BioII/lulab_b/jinyunfan/software/zagros/zagros-1.0.0/bin/thermo"
    
    logger.info("Predict structure ...")
    structurePath  = os.path.join(args.outdir,"sequences.str")
    subprocess.run([thermalExePath,"-o",structurePath,fastaPath])
    logger.info("Done .")
    
    logger.info("Run zagros motif finder ...")
    output = os.path.join(args.outdir,"results.txt")
    cmds = [zagrosExePath,"-output",output,"-structure",structurePath]

    if args.length is not None:
        cmds += ["-width",str(args.length)]
    if args.number is not  None:
        cmds += ["-number",str(args.number)]
    cmds += [fastaPath]
    
    subprocess.run(cmds)
    logger.info("Done .")

     
    
    
def MEME(args):
    logger.info("Prepare input sequences ...")
    fastaPath = prepareFasta(args,T2U=True)
    logger.info("Done .")
    resultsDir = os.path.join(args.outdir,"results")
    #if not os.path.exists(resultsDir):
    #    os.mkdir(resultsDir)
    cmds = ["meme",fastaPath,"-oc",resultsDir,"-rna"]
    if args.number is not None:
        cmds += ["-nmotifs",str(args.number)]
    #if args.length is not None:
    #    cmds += ["-w",str(args.length)]
    if args.min_length is not None:
        cmds += ["-minw",str(args.min_length)]
    if args.max_length is not None:
        cmds += ["-maxw",str(args.max_length)]
    logger.info("Run meme ...")
    logger.info(" ".join(cmds))
    subprocess.run(cmds)
    logger.info("Done .")


def glam2(args):
    logger.info("Prepare input sequence ...")
    fastaPath = prepareFasta(args,T2U=True)
    logger.info("Done .")
    cmds = ["glam2","-O",os.path.join(args.outdir,"results")]
    if args.min_length is not None:
        cmds += ["-a",str(args.min_length).split(".")[0]]
    if args.max_length is not None:
        cmds += ["-b",str(args.max_length).split(".")[0]]
    #if args.length:
    #    cmds += ["-w",str(args.length)]
    cmds += ["n",fastaPath]
    logger.info("Run glam2 motif finder ...")
    logger.info(" ".join(cmds))
    subprocess.run(cmds)
    logger.info("Done .")   


def main():
    parser = argparse.ArgumentParser(description='Run motif finder')
    parser.add_argument('--algorithm','-a',type=str,required=True,
        help="Which method to use",choices=["CMfinder","BEAM","RNApromo","GraphProt","MEMERIS","ssHMM","RNAcontext","MEME","glam2","zagros"])
    parser.add_argument('--indir','-i',type=str,required=True,help="Input dir contains input sequence.")
    parser.add_argument('--outdir','-o',type=str,required=True,help='Output dir of the result.')
    parser.add_argument('--min-length','-m',type=int,help="Min length for motif to find.")
    parser.add_argument('--max-length','-M',type=int,help="Max length for motif to find.")
    parser.add_argument("--length","-l",type=int,help="Length of the motif.")
    parser.add_argument('--number','-n',type=int,default=1,help="Number of motif to find.")
    parser.add_argument("--use-reference-structure","-r",action="store_true",default=False,help="Whether use reference structure, applicable for RNApromo and BEAM.")
    parser.add_argument("--log",type=str,help="Where to output log information.")
    args = parser.parse_args()
    if args.length is not None:
        if args.min_length is None:
            args.min_length = int(args.length)*0.8
        if args.max_length is None:
            args.max_length = int(args.length)*1.2
    if args.min_length is not None and args.max_length is not None:
        if args.min_length > args.max_length:
            logger.info("The specified min motif length is larger than max motif length.")
            sys.exit(1)
    if args.min_length is not None and args.max_length is not None:
        if args.length is None:
            args.length = int((args.min_length+args.max_length)/2)
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    print(args.min_length,args.max_length)
    globals()[args.algorithm](args)


if __name__ == "__main__":
    main()
