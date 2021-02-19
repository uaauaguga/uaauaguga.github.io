import numpy as np


def loadSHAPE(path):
    """
    Load SHAPE reactivities
    """
    SHAPE = {}
    #print("#######These RNAs have NaN lt 0.4######")
    with open(path) as f:
        for line in f:
            assert line.startswith(">")
            name = line.strip().replace(">","")
            reactivity = np.array(next(f).strip().split(",")).astype(float)
            SHAPE[name] = reactivity
            #if np.isnan(reactivity).sum()/reactivity.shape[0]<0.4:
                 #print(name)
    #print("#####################")
    return SHAPE

def prepareSHAPE(shape,prefix):
    shape_path = prefix+".shape"
    with open(shape_path,"w") as f:
        for i,s in enumerate(shape):
            if np.isnan(s):
                s = -999
            f.write("{}\t{}\n".format(i+1,s))
    return shape_path

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

