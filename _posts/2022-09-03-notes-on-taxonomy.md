---
layout: post
title:  "Notes on Taxonomy"
date:   2022-09-03 10:39:16 +0800
usemathjax: true
categories: jekyll update
---


- Download mapping between genbank id and taxonomy <ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid>

- Extract taxnomy information from blast database

```bash
blastdbcmd -entry_batch seq_ids.txt -db nt -dbtype nucl -out taxo.info -outfmt "%a %i %T %L"
# AC189400.2 gb|AC189400.2| 51351 Chinese cabbage
# %a means accession
# %g means gi
# %l means sequence length
# %T means taxid
# %X means leaf-node taxids (what does this mean?)
# %i means sequence id
# %L means common taxonomic name
# %C means common taxonomic names for leaf-node taxids
```


- A very useful module in [ete](http://etetoolkit.org/) toolkits <http://etetoolkit.org/docs/2.3/reference/reference_ncbi.html>

```python
#!/usr/bin/env python
import pickle 
import os
def main():
    from ete3 import NCBITaxa,Tree
    # taxo id of species you want to extract
    species_ids = open("taxo-ids-unique.txt").read().strip().split("\n")
    species_ids = set(species_ids)
    if not os.path.exists("plant-tree.pkl"):
        print("Load ncbi taxo tree ...")
        ncbi = NCBITaxa()
        print("Get plant subtree ...")
        # get descendants of a nodes
        plant_taxo_ids = ncbi.get_descendant_taxa("58024") 
        tree = ncbi.get_topology(plant_taxo_ids,intermediate_nodes=True)
        print("Saving subtree to pickle file .")
        f = open("plant-tree.pkl","wb")
        pickle.dump(tree,f)
        f.close()
    else:
        with open("plant-tree.pkl","rb") as f:
            tree = pickle.load(f)
    used_nodes = []
    genera_counter = {}
    # traverse the tree to extract desired nodes
    for node in tree.traverse():
        if node.rank == "species":
            if node.name in species_ids:
                node.name = node.sci_name.replace(" ",".")
                genus = node.sci_name.split(" ")[0]
                if genus not in genera_counter:
                    genera_counter[genus] = 1
                else:
                    genera_counter[genus] += 1
                    if genera_counter[genus] > 2:
                        continue
                used_nodes.append(node)
    # only keep a subset of nodes and corresponding topology in extracted tree
    tree.prune(used_nodes)
    print(tree.write())


if __name__ == "__main__":
    main()
```

- <https://bioinf.shenwei.me/taxonkit/> seems worth a try


- Visualization of taxonomy tree
  - [itol](https://itol.embl.de/): interactive tree of life, a web tool
  - [ggtree](https://yulab-smu.top/treedata-book/)