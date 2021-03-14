#!/usr/bin/env python
from weblogo import *
import numpy as np
import argparse
import re
from collections import defaultdict
from utilities import loadRecords
import pandas as pd
import copy
np.seterr(all='raise')
np.random.seed(0)
sequenceAlphabet = {"A":0,"C":1,"G":2,"T":3}
structureAlphabet = {"U":0,"P":1}

def oneHotEncoding(s,alphabet=sequenceAlphabet):
    s = s.replace("N","A")
    idx = np.array([alphabet[c] for c in s])
    seq = np.zeros((len(alphabet.keys()),idx.shape[0]))
    seq[idx,np.arange(idx.shape[0])] = 1
    return seq

def prob2str(p,cutoff=0.5):
    p = np.array(p)
    s = np.empty_like(p,dtype=str)
    s.fill("U")
    s[p>0.5] = "P"
    return "".join(s)

def getPFM(counts,pseudoCount = 0.1):
    """
    Input:
    Onehot encoded seed sequence of width w
    Output:
    The PFM
    """
    PFM = counts + pseudoCount*np.ones_like(counts)
    return PFM/PFM.sum(axis=0).reshape((1,-1))

def getLogLikelihood(sequence,start,width,PWM,bgModel):
    if start >= 0:
        motifLikelihood = np.log((sequence[:,start:start+width]*PWM).sum(axis=0)).sum()
        bgLikelihood = np.log((bgModel*sequence[:,:start]).sum(axis=0)).sum() + np.log((bgModel*sequence[:,start+width:]).sum(axis=0)).sum() 
        return motifLikelihood + bgLikelihood
    else:
        noMotifLikelihood = np.log((bgModel*sequence).sum(axis=0)).sum()
        return noMotifLikelihood


def E_step(dataset,sequencePFM,structurePFM,sequenceBg,structureBg,gamma,alpha=0.8):
    # alpha: weighting of sequence in the model
    width = sequencePFM.shape[1]
    ZDict = dict()
    # ZDict: store Zs of each sequence
    # Zs[i] is conditional probability of a motif is starting from i nt of a sequence
    expectedJointLogLikelihood = 0
    for seq_id in dataset.keys():
        m = dataset[seq_id]["sequence"].shape[1]-width
        sequenceLoglikelihoods = []
        structureLoglikelihoods = []
        for j in range(m):
            sequenceLoglikelihood = getLogLikelihood(dataset[seq_id]["sequence"],j,width,sequencePFM,sequenceBg)
            structureLoglikelihood = getLogLikelihood(dataset[seq_id]["pairing"],j,width,structurePFM,structureBg)
            sequenceLoglikelihoods.append(sequenceLoglikelihood)
            structureLoglikelihoods.append(structureLoglikelihood)

        sequenceLoglikelihoods = np.array(sequenceLoglikelihoods).reshape(1,-1)
        structureLoglikelihoods = np.array(structureLoglikelihoods).reshape(1,-1)
        combinedLogLikelihoods = sequenceLoglikelihoods*alpha + structureLoglikelihoods*(1-alpha)
        combinedLikelihoods = np.exp(combinedLogLikelihoods)
        
        sequenceNullLoglikelihood = getLogLikelihood(dataset[seq_id]["sequence"],-1,width,sequencePFM,sequenceBg)
        structureNullLoglikelihood = getLogLikelihood(dataset[seq_id]["pairing"],-1,width,structurePFM,structureBg)
        nullCombinedLogLikelihood = alpha*sequenceNullLoglikelihood+(1-alpha)*sequenceNullLoglikelihood
        nullCombinedLikelihood = np.exp(nullCombinedLogLikelihood)

        ZDict[seq_id] = np.concatenate([(gamma/m)*combinedLikelihoods,np.array([[(1-gamma)*nullCombinedLikelihood]])],axis=1)
        ZDict[seq_id] = (ZDict[seq_id]/ZDict[seq_id].sum()).reshape(1,-1)
        Zs = ZDict[seq_id][0,:-1]
        Q = ZDict[seq_id][0,-1]

        # Calculate contribution of each sequence to joint likelihood
        expectedJointLogLikelihood += np.dot(combinedLogLikelihoods,Zs.T)
        expectedJointLogLikelihood += (1-Q)*nullCombinedLogLikelihood
        expectedJointLogLikelihood += (1-Q)*np.log(1-gamma)
        expectedJointLogLikelihood +=  Q*np.log(gamma/m)
    return ZDict,expectedJointLogLikelihood[0]
 
def M_step(dataset,ZDict,sequencePFM,structurePFM,allSequenceCounts,allStructureCounts):
    sequenceMotifCounts = np.zeros_like(sequencePFM)
    structureMotifCounts = np.zeros_like(structurePFM)
    gamma = 0
    width = sequencePFM.shape[1]
    for seq_id in dataset.keys():
        m = dataset[seq_id]["sequence"].shape[1]-width
        for j in range(m):
            sequenceMotifCounts += dataset[seq_id]["sequence"][:,j:j+width]*ZDict[seq_id][0,j]
            structureMotifCounts += dataset[seq_id]["pairing"][:,j:j+width]*ZDict[seq_id][0,j]
        gamma += ZDict[seq_id][0,:m].sum()
    gamma_pseudocount = 1e-4
    gamma = gamma/(len(dataset.keys())+gamma_pseudocount)
    sequencePFM = getPFM(sequenceMotifCounts)
    structurePFM = getPFM(structureMotifCounts)
    bgSequenceCounts = allSequenceCounts - sequenceMotifCounts.sum(axis=1).reshape(-1,1)
    bgStructureCounts = allStructureCounts - structureMotifCounts.sum(axis=1).reshape(-1,1)
    sequenceBg = bgSequenceCounts/bgSequenceCounts.sum()
    structureBg = bgStructureCounts/bgStructureCounts.sum()       
    return sequencePFM,structurePFM,sequenceBg,structureBg,gamma


def randomProjectionInit(dataset,width,size=5,init="structure"):
    global sequenceAlphabet
    global structureAlphabet
    bucket = defaultdict(list)
    sampled = np.random.randint(0,width,size)
    allSequenceCounts = np.zeros((len(sequenceAlphabet.keys()),1))
    allStructureCounts = np.zeros((len(structureAlphabet.keys()),1))
    sequenceCounts = np.zeros((len(sequenceAlphabet.keys()),width))
    structureCounts = np.zeros((len(structureAlphabet.keys()),width))
    for seq_id in dataset.keys():
        allSequenceCounts += oneHotEncoding(dataset[seq_id]["sequence"],sequenceAlphabet).sum(axis=1).reshape((len(sequenceAlphabet.keys()),1))
        allStructureCounts += oneHotEncoding(dataset[seq_id]["pairing"],structureAlphabet).sum(axis=1).reshape((len(structureAlphabet.keys()),1))
        if init == "structure":
            string = dataset[seq_id]["pairing"]
        else:
            string = dataset[seq_id]["sequence"]
        for pos in range(len(string)-width):
            wmer = string[pos:pos+width]
            projected = ""
            for i in sampled:
                projected += wmer[i]
            bucket[projected].append((seq_id,pos))
    L = 0
    for wmer in bucket.keys():
        if len(bucket[wmer]) > L:
            if np.unique(list(wmer)).shape[0] == 1:
                continue
            L = len(bucket[wmer])
            wmer_ = wmer
    print(wmer_,len(bucket[wmer_]))
    seq_ids_used = set()
    for seq_id,pos in bucket[wmer_]:
        if seq_id in seq_ids_used:
            continue
        sequenceCounts += oneHotEncoding(dataset[seq_id]["sequence"][pos:pos+width],sequenceAlphabet)
        structureCounts += oneHotEncoding(dataset[seq_id]["pairing"][pos:pos+width],structureAlphabet)
        seq_ids_used.add(seq_id)
    sequencePFM = getPFM(sequenceCounts)
    structurePFM = getPFM(structureCounts)
    bgSequenceCounts = allSequenceCounts - sequenceCounts.sum(axis=1).reshape((-1,1))
    bgStructureCounts = allStructureCounts - structureCounts.sum(axis=1).reshape((-1,1))
    bgSequenceModel = bgSequenceCounts/bgSequenceCounts.sum()
    bgStructureModel = bgStructureCounts/bgStructureCounts.sum()
    return sequencePFM,structurePFM,bgSequenceModel,bgStructureModel,allSequenceCounts,allStructureCounts


def runMotifFinding(dataset,width=10,gamma=0.9,init ="structure",seed_size=5,alpha=0.9):
    """
    Input: dict contains sequence name and onehot encoded sequences
    Output: PFM
    """
    global sequenceAlphabet
    global structureAlphabet
    sequencePFM,structurePFM,sequenceBg,structureBg,allSequenceCounts,allStructureCounts = randomProjectionInit(dataset,width,size=seed_size,init =init)
    for seq_id in dataset.keys():
        dataset[seq_id]["sequence"] = oneHotEncoding(dataset[seq_id]["sequence"],sequenceAlphabet)
        dataset[seq_id]["pairing"] = oneHotEncoding(dataset[seq_id]["pairing"],structureAlphabet)

    n_epoch = 0

    while True:
        n_epoch += 1
        # E step: get the expectation of each position to become the start of a motif
        ZDict,expectedJointLogLikelihood = E_step(dataset,sequencePFM,structurePFM,sequenceBg,structureBg,gamma,alpha=alpha)
        # M step: update PFM
        print(f"EM round {n_epoch} logl: {expectedJointLogLikelihood} gamma: {gamma}")
        sequencePFM,structurePFM,sequenceBg,structureBg,gamma = M_step(dataset,ZDict,sequencePFM,structurePFM,allSequenceCounts,allStructureCounts)
        if n_epoch > 10 and expectedJointLogLikelihood - _expectedJointLogLikelihood < 1e-3:
            break
        _expectedJointLogLikelihood = expectedJointLogLikelihood

    locations = []
    for seq_id in ZDict:
        start = ZDict[seq_id].argmax()
        if start > dataset[seq_id]["sequence"].shape[1] - width:
             print(f"No motif detected in {seq_id}")
             continue
        #print(seq_id,start,start+width,sep="\t")
        locations.append((seq_id,start,start+width))   
    return sequencePFM,structurePFM,locations
 

def drawLogo(hits,path):
    logoData = LogoData.from_seqs(hits)
    logoOptions = LogoOptions()
    logoFormat = LogoFormat(logoData, logoOptions)
    pngbytes = png_formatter(logoData,logoFormat)
    with open(path,"wb") as f:
        f.write(pngbytes)

def main():
    parser = argparse.ArgumentParser(description='Discovery sequence motif using EM algorithm')
    parser.add_argument('--input','-i',required=True,help="Fasta like file with sequence and pairing probability")
    parser.add_argument('--width','-w',type=int,required=True,help="Width of the motif")
    parser.add_argument('--gamma','-g',type=float,default=0.9,help="Prior of gamma, the fraction of sequences that contains the motif")
    parser.add_argument('--alpha','-a',type=float,default=0.9,help="Weighting of sequence in log likelihood")
    parser.add_argument('--seed-size','-sz',type=int,default=5,help="Seed size for random projection")
    parser.add_argument('--initialize','-init',default="sequence",choices=["sequence","structure"],help="Seed size for random projection")
    parser.add_argument('--locations','-loc',required=True,help="Path of motif location")
    parser.add_argument('--sequence-logo',help="Path of sequence logo")
    parser.add_argument('--structure-logo',help="Path of structure logo")
    args = parser.parse_args()
    dataset = loadRecords(args.input,order="sequence,pairing-probability",dtype="str,float")

    for seq_id in dataset:
        dataset[seq_id]["pairing"] = prob2str(dataset[seq_id]["pairing-probability"])

    dataset0 = copy.deepcopy(dataset)   
    sequencePFM,structurePFM,locations = runMotifFinding(dataset,width=args.width,gamma=args.gamma,seed_size=args.seed_size,init=args.initialize)
            
    if args.sequence_logo:
        sequenceHits = seq.SeqList(alphabet=seq.Alphabet("ACGT"))
    if args.structure_logo:
        structureHits = seq.SeqList(alphabet=seq.Alphabet("UP"))
        

    with open(args.locations,"w") as f:
        for seq_id,start,end in locations:
            f.write(f"{seq_id}\t{start}\t{end}\n")
            if args.sequence_logo:
                sequenceHits.append(dataset0[seq_id]["sequence"][start:end])
            if args.structure_logo:
                structureHits.append(dataset0[seq_id]["pairing"][start:end])
    
    if args.sequence_logo:
        drawLogo(sequenceHits,args.sequence_logo)
    if args.structure_logo:
        drawLogo(structureHits,args.structure_logo)
        

if __name__ == "__main__":
    main()
