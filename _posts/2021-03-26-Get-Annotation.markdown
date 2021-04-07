---
layout: post
title:  "Get Annotations"
date:   2021-03-25 22:47:20 +0800
usemathjax: true
categories: jekyll update
---

## Play with Genome annotation

### Where to download genome annotation in gtf/gff

- [Gencode](https://www.gencodegenes.org/)
  - mice and human
- [Ensembl](https://uswest.ensembl.org/info/about/species.html)

### Annotation resources in bioconductor

- Some material for reading
- <https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf> 
- <https://www.bioconductor.org/help/course-materials/2019/BSS2019/05_Annotations.html>
- Also see discussion here <https://www.biostars.org/p/287871/> for difference between biomart and bioconductor annotation resources

### Biomart
- The biomaRt package is simply a way to programmatically access the Biomart server and get the results back into R


### Mapping between homolog
- <https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/>


## Functional annotation

- So call functional annotation is actually mappings from genes to functions
- If we have a gene set, to see whether they are related to a function, we use genes that annotated to such function, and use the overlap between known gene set and given gene set as a proxy for functional relevance
- If we have a ranked list, we can do similar things

### Get KEGG annotations with REST API
- Read <https://www.kegg.jp/kegg/rest/keggapi.html> for the native version of KEGG's REST API
- Several packages provides wrappers for these API

- In Biopython

```python
from Bio.KEGG import REST
from tqdm import tqdm
import re

# Get KEGG pathway id and description of human
pathWayIds = REST.kegg_list("pathway","hsa").read().strip().split("\n")

# Retrieve data of a KEGG pathway
pathWayId = pathWayIds[0].split("\t")[0].split(":")[1]
# hsa00010
pathWayData = REST.kegg_get(pathWayId).read()
# A structured string,quite complicate
# Here we get entrez id from this string as an example
def getEntrezId(result):  
    entrezs = []
    flag = 0
    for line in result.split("\n"):
        data = re.split("\s+",line)
        if len(data)<2:
            continue
        if data[0]=="GENE":
            flag = 1
            entrezs.append(data[1])
        elif flag == 1:
            if data[0]=="":
                entrezs.append(data[1])
            else:
                break
    return entrezs

entrezIds = getEntrezId(pathWayData)
```

- In bioconductor package [KEGGREST](http://bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html)

```R
library(KEGGREST)
homo.pathways <- keggList("pathway","hsa")
pathway.id <- names(homo.pathways)[1]
pathway.data <- keggGet(c(pathway.id))
# keggGet parsed the returned string to a list
entrez.ids <- pathway.data[[1]]$GENE[c(T,F)]
# Here the [c(T,F)] trick take element with odds index
```

### Get GO annotation in R

- See discussion here <https://www.biostars.org/p/52101/#68158>

```R
library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
go.ids <- c(names(goterms)[1:2])
go.gene.list <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=go.ids, column="ENTREZID",
                                  keytype="GOALL", multiVals='list')
```

- For enrichment of non-model organism, see [Annotationhub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html)



```R
library(AnnotationDbi)
library(AnnotationHub)
hub <- AnnotationHub()
# Query species you want 
dm <- query(ah, c("cgriseus"))
# Retrieve database from dm object with its id
org.cgriseus.db <- hub[["AH70586"]]
# THen simply replace org.Hs.eg.db with your database for retrieve go term
```



- This [this post](https://guangchuangyu.github.io/cn/2017/07/clusterprofiler-maize/) for further reading

### Softwares and tools for enrichment analysis

- Under most situation, its not necessary to download gene set and perform hypergeometric test your self for enrichment analysis
- Here list two useful tools
- [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler)
  - More bioinformatician oriented
- [metascape](https://metascape.org/gp/index.html)
  - Web tool, more bench scientist oriented