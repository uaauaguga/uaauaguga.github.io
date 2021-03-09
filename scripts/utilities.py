import numpy as np
from collections import defaultdict
import re

def loadSHAPE(path):
    """
    Load SHAPE reactivities
    """
    SHAPE = {}
    with open(path) as f:
        for line in f:
            assert line.startswith(">")
            name = line.strip().replace(">","")
            reactivity = np.array(next(f).strip().split(",")).astype(float)
            SHAPE[name] = reactivity
    return SHAPE


def loadRecords(path,data_type ="sequence,structure,reactivity"):
    data_types = data_type.split(",")
    n_data_types = len(data_types)
    records = defaultdict(dict)
    with open(path) as f:
        for line in f:
            try:
                assert line.startswith(">")
            except:
                print("The input data is not consistent to specified data_type parameter")
            seq_id = line.strip()[1:].strip()
            for data_type in data_types:
                data = next(f).strip()
                if data_type == "reactivity":
                    data = np.array(data.split(",")).astype(float)
                records[seq_id][data_type] = data
    return records


def writeRecords(dataDict,path,order="sequence,structure,reactivity"):
    order = order.split(",")
    with open(path,"w") as f:
        for seq_id in dataDict.keys():
            f.write(">"+seq_id+"\n")
            for key in order:
                if isinstance(dataDict[seq_id][key],str):
                    line = dataDict[seq_id][key]
                else:
                    data = np.round(np.array(dataDict[seq_id][key]),3)
                    data = data.astype(str)
                    line = ",".join(data)
                f.write(line+"\n")
            


def prepareSHAPE(shape,prefix):
    """
    Prepare SHAPE file for RNAstructure and ViennaRNA package
    """
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
    Each sequence records could span multiple lines
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



def checkDBN(s):
    """
    Check whether a dot bracket notation represents a valid secondary structure
    """
    stack = []
    for c in s:
        if c == "(":
            stack.append(c)
        elif c == ")":
            assert len(stack) > 0
            _ = stack.pop()
        else:
            continue
    assert len(stack) == 0


def getKmer(structure,reactivity,flankingL = 2):
    """
    Get k-mer, where k = flankingL*2 + 1
    structure is a string
    reactivity is a one dimensional numpy array with same length
    return:
    reactivity at positions where reactivity is not None
    structure k-mer flanking these positions
    """
    indices = np.where(~np.isnan(reactivity))[0]
    indices = indices[(indices >= flankingL)&(indices<(reactivity.shape[0]- flankingL))]

    structures_ = []
    for idx in indices:
        s = structure[idx-flankingL:idx+flankingL+1]
        structures_.append(s)
    reactivity_ = reactivity[indices]

    return structures_,reactivity_


def annotatePairing(dataDict):
    """
    Input:
    a dict of dict , sequence id as keys
    Items in the inner dict should contain a 'structure' key, 
    and the corresponding value corresponds to secondary structure in  dot bracket notation 
    Return:
    the input dict, with an additional 'pairing' field, 
    'P' indicate this position is paired
    'U' indicate this position is unpaired
    """
    for seq_id in dataDict.keys():
        annotation = dataDict[seq_id]["structure"].replace(".","U")
        annotation = re.sub("[^U]","P",annotation)
        dataDict[seq_id]["pairing"] = annotation
    return dataDict

def pairingDigitToString(X):
    """
    Convert pairing state of a set of k-mer into string
    Input:
    numpy array of shape (n,k) with binary value (0 for unpaired, 1 for paired)
    where n is number of k mer
    k is length of k-mer
    Return:
    numpy array of length n
    With item are length k string [PU]+ 
    """
    s = np.empty_like(X).astype(str)
    s[X==1],s[X==0] = "P","U"
    return np.apply_along_axis(lambda x:"".join(x),arr=s,axis=1)
