#!/usr/bin/env python
import argparse 
import logging
from itertools import product
import numpy as np
import json
from sklearn.metrics import roc_curve,auc
from lightgbm import LGBMClassifier
import pickle
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("evaluate predictive power of k-mer profile")


def load_fasta(path):
    sequences = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                seq_id  = line.strip()[1:].split(" ")[0]
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip().upper().replace("T","U")
    seq_ids = list(sequences.keys())
    for seq_id in seq_ids:
        if "N" in sequences[seq_id]:
            del sequences[seq_id]
    return sequences


def get_kmer_profile(sequences):
    seq_ids = []
    kmer_profiles = []
    token_size = len(tokens_lut)
    omited_tokens = set()
    for seq_id in sequences:
        seq_ids.append(seq_id)
        x = np.zeros(token_size,dtype=int)
        sequence = sequences[seq_id]
        for i in range(len(sequence)-k):
            token = sequence[i:i+k]
            if token not in tokens_lut:
                omited_tokens.add(token)
                continue
            x[tokens_lut[token]] += 1
        kmer_profiles.append(100*x/len(sequence))
    kmer_profiles = np.array(kmer_profiles)
    logger.warning("These tokens were omited:")
    logger.warning(",".join(omited_tokens))
    return seq_ids, kmer_profiles


def main():
    parser = argparse.ArgumentParser(description='sampling negative sequences')
    parser.add_argument('--input', "-i",type=str, required=True, help='positive instance for training')
    parser.add_argument('--output', "-o",type=str, required=True, help='predicted probabilities')
    parser.add_argument('--model','-m',type=str,help="Where to load the model (in pickle format)")
    args = parser.parse_args()
    tokens = []

    logger.info("Load model ...")
    with open(args.model,"rb") as f:
        clf = pickle.load(f)

    global k
    k = int(np.log2(clf.n_features_)/2)
    logger.info(f"The number of input features is {clf.n_features_}, that means k={k}")
    for token in product("ACGU",repeat=k):
        tokens.append("".join(token))
    global tokens_lut
    tokens_lut = {}
    for i, token in enumerate(sorted(tokens)):
        tokens_lut[token] = i
    logger.info("Load sequences ...")
    sequences = load_fasta(args.input)
    logger.info(f"Calculate {k}-mer profile ...")
    seq_ids, kmer_profiles = get_kmer_profile(sequences)

    logger.info("Make predictions ...")
    y_pred_proba = clf.predict_proba(kmer_profiles)[:,1]


    logger.info("Saving results ...")
    
    with open(args.output,"w") as f:
        for seq_id, score in zip(seq_ids,list(y_pred_proba)):
            f.write(f"{seq_id}\t{score}\n")
    logger.info("All done .")

if __name__ == "__main__":
    main()
