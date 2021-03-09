#!/usr/bin/env python
import sys
sys.path.append('scripts')
from utilities import writeRecords,loadRecords 
import os
def main():
    structures = loadRecords("test/motif-finder-demo/data/RF00037-100-hardconstraint-ViennaRNA.dot","sequence,structure")
    bpp = loadRecords("test/motif-finder-demo/data/RF00037-100-ViennaRNA-with-SHAPE.txt","pairing-probability")
    for seq_id in structures:
        structures[seq_id]["pairing-probability"] = bpp[seq_id]["pairing-probability"]
    writeRecords(structures,"test/motif-finder-demo/data/RF00037-100-ViennaRNA-with-SHAPE-demo.txt",order="sequence,structure,pairing-probability")


if __name__ == "__main__":
    main()
