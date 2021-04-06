---
layout: post
title:  "Dive into RNA-seq Quantification"
date:   2021-03-31 07:05:55 +0800
usemathjax: true
categories: jekyll update
---

- How to quantify RNA-seq data and perform differential analysis
- Some strange phenomenon we may meet

## Expression quantification

- Reads mapping is not required for expression level quantification

### Duplication handling
- Seems it has become a consensus that we've better not remove duplicated reads in RNA-seq when there is no UMI
- <https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/>
- <https://www.biostars.org/p/55648/>


### Exon level

- Count directly
- Some tools like [DEXSeq](https://www.bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html) use counts at exon level to infer differential exon usage

### Transcript level / isoform level

- Counts at isoform level cannot be resolved precisely from RNA-seq data as the read is short, and there are large overlap between different isoforms of same gene
- Multiple tools is able to give a estimated TPM or FPKM at transcript level
  - [Rsem](https://deweylab.github.io/RSEM/)
  - [Salmon](https://combine-lab.github.io/salmon/)
  - [Kallisto](https://pachterlab.github.io/kallisto/)

### Gene level

- Count directly
- Summarize estimated expression level at transcripts level to gene
  - [tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html) is very useful for such task



### Normalization

- We shall apply different normalization methods for different propose
- For differential expression, here shows case in edgeR, limma and DEseq2
- Differential expression compare difference between samples, hence we shall pay no attention to things that believed to be the same across different samples 
  - GC content
  - Gene length

### Variance stabilizing  transformation



### Some background in statistics

#### **Over-dispersed Poisson / Poisson gamma mixture / Negative binomial model**

- Poisson GLMs

  $$P(y;\mu) = \frac{e^{-\mu}\mu^{y}}{y!}$$

- The most common link function used for Poisson glms is the logarithmic
link function, that means, fit $$log(\mu)$$ with a linear model
- In Poisson GLMs, we have $$Var(y) = \mu$$, in practice the apparent
variance of the data often exceeds $$\mu$$, that is called "overdispersion"
- See <http://llc.stat.purdue.edu/2014/41600/notes/prob1804.pdf> for deduction of Poisson mean and variance
- The presence of overdispersion will lead to under estimation of standard deviation
- How to detect overdispersion?
- Hierarchical model: instead of assuming $$y_{i}{\sim}Pois(\mu_{i})$$, we can add a second layer of variability by allowing
$$μ_{i}$$ itself to be a random variable. We suppose, where $$G(\mu_{i},\psi)$$ denotes a distribution with mean $$\mu_{i}$$ and **[coefficient of variation](https://en.wikipedia.org/wiki/Coefficient_of_variation)** $$\psi$$

$$y_{i}|\lambda_{i}{\sim}Pois(\lambda_{i})$$

$$\lambda_{i}{\sim}G(\mu_{i},\psi)$$

- According to [law of total expectation](https://en.wikipedia.org/wiki/Law_of_total_expectation)

$$E(y) = \mathop{E}_{\lambda}[\mathop{E}_{y|\lambda}(y|\lambda)] = \mu$$

- According to [law of total variance](https://en.wikipedia.org/wiki/Law_of_total_variance)


$$Var(Y) = E(Var(Y|X)) + Var(E(Y|X))$$

$$Var(y) = E(Var(y|\lambda)) + Var(E(y|\lambda)) = \mu + \mu^2\psi$$

- So the variance contains an overdisperion term $$\psi\mu^2$$ . The larger $$\psi$$, the greater
the overdispersion.
- When $$G$$ is a gamma distribution, we call the overall distribution negative binomial
- Negative bionomial distribution can be explicitly written as: (Here $$k=\frac{1}{\psi}$$)

$$P(y_i;\mu_i,k) = \frac{\Gamma(y_i+k)}{\Gamma(y_i+1)\Gamma(k)}(\frac{\mu_i}{\mu_i+k})^{y_i}(1-\frac{\mu_i}{\mu_i+k})^k  $$ 


#### **Empirical bayes and moderated variance**

#### **Design matrix and contrast in R**
- Factors and levels in R

```R
#If level is not specified, the level is determined by R as alphabet order
factor(c("c","b","a","c","c","b"))
#[1] c b a c c b
#Levels: a b c
factor(c("c","b","a","c","c","b"),levels = c("c","b","a"))
#[1] c b a c c b
#Levels: c b a
```
- Reset reference level

```R
groups <- factor(c("c","b","a","c","c","b"),levels = c("a","c","b"))
# Levels: a c b
relevel(groups,ref="c")
# Moce c to first level, or reference level
# Levels: c a b
```

- Design matrix in R

```R
# The first level is the reference level in design matrix
groups <- factor(c("c","b","a","c","c","b"),levels = c("a","c","b"))
design.matrix <- model.matrix(~groups,levels = c("c","b","a"))
#  (Intercept) groupsc groupsb
#1           1       1       0
#2           1       0       1
#3           1       0       0
#4           1       1       0
#5           1       1       0
#6           1       0       1
# l - 1 binary value is needed for coding l levels of categorical data
# Three columns in the output dataframe, first columns corresponds to intercept term
# Next 2 columns corresponds to factors other then reference level
```


```R
gender <- c("M","M","M","M","F","F","F","U","U")
country <- c("B","C","A","B","C","A","A","B","C")
metainfo <- data.frame(gender=gender,country=country)
# According to R's alphabet order behavior, F and A is the reference level
design.matrix <- model.matrix(~gender+country,metainfo)
# 5 columns: 1 (intercept) + (3-1) (gender) + (3-1) (country)
```

- Read more at <http://genomicsclass.github.io/book/pages/interactions_and_contrasts.html>


### Differential analysis (a paired comparison example)

#### Load counts and design matrix
- Note only base R fucntion are used

```R
count.path <- "TCGA-CRC-counts.txt"
count.matrix <- read.table(count.path, header = T, sep ="\t",row.names=1,stringsAsFactors = F,check.names = F)
count.matrix <- as.matrix(count.matrix)
sample.ids <- colnames(count.matrix)
tissue.types <- rep("normal",length(sample.ids))
# level=c("normal","tumor") means normal is the reference level
# 注意别搞反了！
# If level is not specified, the level is determined by R as alphabet order
tissue.types[grep("01A$",sample.ids)] <- "tumor"
tissue.types <- factor(tissue.types,level=c("normal","tumor"))
patient.ids <- factor(unlist(lapply(sample.ids,function(x){strsplit(x,"-")[[1]][3]})))
# model.matrix is an R function
# If you'd like to consider some batch effect, just add these factors to desig matrix
design <- model.matrix(~tissue.types+patient.ids)
```

#### **edgeR**
#### Prepare DGEList

```R
library(edgeR)
y <- edgeR::DGEList(counts=count.matrix)
keep <- edgeR::filterByExpr(y,group = tissue.types)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- edgeR::calcNormFactors(y,method="TMM")
```

- The `DGEList` S4 object contains `counts` and `samples`
- `counts` is a matrix with integer value
- `samples` is a data.frame with three fields:
  - `group`: introduce in `filterByExpr`
  - `lib.size`
  - `norm.factors`: calculated in `edgeR::calcNormFactors`, all ones if not perform this step 


#### Dispersion estimation and model fitting
- We have shown (gene `g` in sample `i`), CV means coefficient of variation, or
standard deviation divided by mean

$$\frac{Var(y_{gi})}{\mu_{gi}^2} = \frac{1}{\mu_{gi}} + \phi_{g}$$

$$CV^2(y_{gi})= \frac{1}{\mu_{gi}} + \phi_{g}$$

- $$\phi_{g}$$ is called dispersion, suppose technical varistion follows Poisson law, then $$\phi_{g}$$ equals to square to Biological CV

$$CV_{Total}^2 = CV_{Technical}^2 + CV_{Biological}^2$$

- Estimate biological CV is equivalent to estimate dispersion parameter in negative bionomial distribution

```R
y <- estimateDisp(y, design)
# This is actually short for the following three steps
# Estimate trended dispersion, add a common.dispersion term to the DEGList
#y <- estimateGLMCommonDisp(y, design)
# Estimate trended dispersion, add a trended.dispersion term to the DEGList
#y <- estimateGLMTrendedDisp(y, design)
# tagwise NB dispersion, add a tagwise.dispersion term to the DEGList
#y <- estimateGLMTagwiseDisp(y, design)

# Attempt to use trended.dispersion, use common.dispersion if not present
fit.ql <- glmQLFit(y, design)
# The first column is intercept, the second corresponds to coefficient to be tested (tumor vs. normal), the remaining columns are different patients
tumor.vs.normal.ql <- glmQLFTest(fit, coef=2)
tumor.vs.normal.ql.de <- topTags(tumor.vs.normal.ql,n=Inf)

# order of precedence: genewise dispersion, trended dispersions, common dispersion
fit.lr <- glmFit(y, design)
tumor.vs.normal.lr <- glmLRT(fit, coef=2)
tumor.vs.normal.lr.de <- topTags(tumor.vs.normal.lr,n=Inf)

# Summarize number of up and down regulated genes
summary(limma::decideTests(tumor.vs.normal.lr,lfc=1))
```

- Dispersions
  - Estimated by quantile-adjusted conditional maximum likelihood (qCML), not a ML estimation
  - `common.dispersion`: a single numeric value
  - `trended.dispersion`: numeric value vector, with length same as number of gene, what its does is stratifying dispersion estimation by gene abundance
  - `tagwise.dispersion`: numeric value vector, with length same as number of gene
  - `glmQLFit` first tries to use `trended.dispersion`, `glmFit` first tries to use `tagwise.dispersion`
  - `glmQLFTest` is more conservative than `glmLRT`


#### **limma**
- limma was originally designed for microarray data analysis, and was then extended to RNA-seq differential expression and differention splicing analysis <https://academic.oup.com/nar/article/43/7/e47/2414268>
- The data input to limma should be counts, rather than popular expression summaries such as reads-per-kilobase-per-million (RPKM), so that limma can estimate the appropriate mean-variance relationship. 

- Count data always show marked mean-variance relationships. Raw counts show increasing variance with increasing count size, while log-counts typically show a decreasing mean-variance trend. 
- `limma` provides two choices, work on logCPM and consider mean - variance trend, or work on `voom` transformed data without considering mean - variance trend
- If counts is a DGEList object from the edgeR package, then voom will use the normalization factors found in the object when computing the logCPM values. In other words, the logCPM values are computed from the effective library sizes rather than the raw library sizes. That means `voom(y, design)` is not equivalent to `voom(y$counts, design)` if not all `DGEList$samples$norm.factors` are ones

- Prepare DGEList
  - Same as ``edgeR


- Perform differential analysis with limma - trend

```R
# logCPM is a numeric matrix
logCPM <- edgeR::cpm(y, log=TRUE, prior.count=3)
fit.limma <- limma::lmFit(logCPM, design)
# trend: should an intensity-trend be allowed for the prior variance?
fit.limma.trend <- limma::eBayes(fit.limma, trend=TRUE)
de.table.trend <- limma::topTable(fit.limma.trend, coef=2,number=Inf)
```

- Perform differential analysis with limma - voom

```R
v <- limma::voom(y, design, plot=FALSE)
fit.limma.voom <- lmFit(v, design)
# trend is FALSE by default
fit.limma.voom <- limma::eBayes(fit.limma)
de.table.voom <- limma::topTable(fit, coef=2,number=Inf)
```


#### DESeq2
- <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>

```R
library(DESeq2)
coldata <- data.frame(tissue.types=tissue.types,patient.ids=patient.ids)
dds <- DESeqDataSetFromMatrix(countData = count.matrix,
                              colData = coldata,
                              design= ~ tissue.types+patient.ids)
dds <- DESeq(dds)
#DESeq is short for the following three function

# 1) dds <- estimateSizeFactors(dds)
# Estimation of size factor 

# A size factor is added to dds as dds@colData$sizeFactor

# 2) dds <- estimateDispersions(dds)
# estimateDispersions add dds@dispersionFunction() for expression level dispersion trend
# estimateDispersions is short for
#   1. estimateDispersionsGeneEst
#   2. estimateDispersionsFit
#   3. estimateDispersionsMAP

# 3) dds <- nbinomWaldTest(dds)
# Perform wald test

# Show comparisons that can be extracted
resultsNames(dds) 
# "Intercept","tissue.types_tumor_vs_normal","patient.ids_2679_vs_2671",...
# Extract differential expresion result
deseq.res <- results(dds,name="tissue.types_tumor_vs_normal")
# An alternative way is
#deseq.res <-results(dds,contrast=c("tissue.types", "tumor", "normal"))
```


### Data Transformation for downstream analysis
- Downstream analysis
  - PCA
  - heatmap
  - clustering
  - machine learning
- Data used for down stream analysis can be quiet different from differential expression
- Simplist solution might be `edgeR::cpm(y, log=TRUE, prior.count=1)`
- Things get more complicated when considering batch effect for task other than differential expression
- Further reading
  - Caveats for correcting for batch effect <https://academic.oup.com/biostatistics/article/17/1/29/1744261>

### Reference

- limma package related
  - [limma user guide](https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)

- edgeR package related
  - [edgeR manual](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

