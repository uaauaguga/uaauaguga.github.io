#!/usr/bin/env python
import subprocess
import argparse
import re
import os
import sys
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('Evaluate BEAM motif finder')



def extendSeqName(trucatedSeqNames,fastaPath):
    fullNames = []
    with open(fastaPath) as f:
        for line in f:
            if line.startswith(">"):
                fullNames.append(line.strip().replace(">",""))
    nameMapping = {}
    for trucatedSeqName in trucatedSeqNames:
        for fullName in fullNames:
            if fullName.find(trucatedSeqName) == 0:
                nameMapping[trucatedSeqName] = fullName
                break
    return nameMapping
        


def BEAM(args):
    def getFinalMotif(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#motif"):
                    motif = next(f)
                    motif = motif.strip()
                    break
                else:
                    continue
        return motif

    def getMotifLocations(path):
        ## Return: sequence name, 0 based [start, end)
        records = []
        with open(path) as f: 
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    break
                if line.startswith("BEAR"):
                    continue
                fields = line.split("\t")
                seqName,start,end = fields[1],int(fields[3]),int(fields[4])
                seqName = re.sub(r"\$\d+$","",seqName)
                records.append((seqName,start,end))
        locations = pd.DataFrame.from_records(records)
        locations.columns = ["sequence_id","start","end"]
        locations = locations.set_index("sequence_id")
        return locations

    finalMotifPath = os.path.join(args.indir,"results/sequences/sequences_summary.txt")
    motifLocationPath = os.path.join(args.indir,"results/sequences/benchmark/motifs/sequences_m1_run1.txt")
    #finalMotif = getFinalMotif(finalMotifPath)
    locations = getMotifLocations(motifLocationPath)
    return locations


def CMfinder(args):
    from Bio import AlignIO
    stkPath = os.path.join(args.indir,"sequences.fa.motif.h1_1")
    # stk format is 1 based
    records = []
    with open("sequences.fa.motif.h1_1") as f:
        for record in AlignIO.read(f, 'stockholm'):
            start,end = record.description.split("\t")[0].split("..")
            start,end = int(start)-1,int(end)
            records.append([record.name,start,end]) 
    locations = pd.DataFrame.from_records(records)
    locations.columns = ["sequence_id","start","end"]
    locations = locations.set_index("sequence_id")
    return locations

def RNApromo(args):
    #start,end+1 
    locationsPath = os.path.join(args.indir,"results","locs_1.xls") 
    locations = pd.read_csv(locationsPath,sep="\t")
    fastaPath = os.path.join(args.indir,"sequences.fa")
    nameMapping = pd.Series(extendSeqName(list(locations.index),fastaPath))
    locations.index = nameMapping.loc[locations.index]
    locations = locations.loc[:,["Start","End"]] 
    locations.index.name = "sequence_id"
    locations.columns = ["start","end"]
    locations["end"] = locations["end"]+1
    return locations


def GraphProt(args):
    locationsPath = os.path.join(args.indir,"results","locations.txt")
    locations = pd.read_csv(locationsPath,sep="\t")
    return locations

def MEMERIS(args):
    from Bio import motifs
    records = []
    fastaPath = os.path.join(args.indir,"sequences.fa")
    with open(os.path.join(args.indir,'results.txt')) as f:
        for motif in  motifs.parse(f,'MEME',strict=False):
            for instance in motif.instances:
                start = instance.start - 1
                end = start + len(str(instance))
                records.append((instance.sequence_name,start,end)) 
    locations = pd.DataFrame.from_records(records)
    locations.columns = ["sequence_id","start","end"]
    locations = locations.set_index("sequence_id")
    nameMapping = pd.Series(extendSeqName(list(locations.index),fastaPath))
    locations.index = nameMapping.loc[locations.index]
    return locations


def ssHMM(args):
    
    pass
 

def RNAcontext(args):
    
    pass

def MEME(args):
    xmlPath = os.path.join(args.indir,"results","meme.xml")
    from xml.etree import ElementTree
    root = ElementTree.parse(xmlPath).getroot()
    id2name = {}
    for seq in root.find("training_set").findall("sequence"):
        data = dict(seq.items())
        id2name[data["id"]] = data["name"]
    records = []
    for motifsNode in root.find("motifs").findall("motif"):
        for locationNode in motifsNode.find("contributing_sites").findall("contributing_site"):
            locInfo = dict(locationNode.items())
            seqName = id2name[locInfo['sequence_id']]
            start = int(locInfo['position'])
            length = len(locationNode.find("site").getchildren())
            end = start + length
            records.append((seqName,start,end))
    locations = pd.DataFrame.from_records(records)
    locations.columns = ["sequence_id","start","end"]
    locations = locations.set_index("sequence_id")
    return locations

def glam2(args):
    records = []
    with open(os.path.join(args.indir,"results","glam2.txt")) as f:
        for line in f:
            if line.startswith("Score:"):
                break
        _ = next(f)
        _ = next(f)
        for line in f:
            line = line.strip()
            if len(line) == 0:
                break
            fields = line.split()
            seqName,start,end = fields[0],int(fields[1])-1,int(fields[3])
            records.append((seqName,start,end))
    locations = pd.DataFrame.from_records(records)
    locations = pd.DataFrame.from_records(locations)
    locations.columns = ["sequence_id","start","end"]
    locations = locations.set_index("sequence_id")
    fastaPath = os.path.join(args.indir,"sequences.fa")
    nameMapping = pd.Series(extendSeqName(list(locations.index),fastaPath))
    locations.index = nameMapping.loc[locations.index]
    return locations
    

def main():
    parser = argparse.ArgumentParser(description='Evaluate performance of BEAM motif finder')
    parser.add_argument('--algorithm','-a',type=str,required=True,
        help="Which method to use",choices=["CMfinder","BEAM","RNApromo","GraphProt","MEMERIS","ssHMM","RNAcontext","MEME","glam2"])
    parser.add_argument('--indir','-i',type=str,required=True,help="Input dir contains result of an motif finder")
    parser.add_argument('--bed','-b',type=str,required=True,help="Output bed file")
    args = parser.parse_args()
    locations = globals()[args.algorithm](args)
    locations.to_csv(args.bed,sep="\t",header=None)

if __name__ == "__main__":
    main()
