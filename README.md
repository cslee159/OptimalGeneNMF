# NMFOptiGene

NMFOptiGene is a R script for finding .... (추가예정)

논문 참고 : description->논문에 언급

## Gene Optimization using NMF based on Missing Value Imputation

### Script Command
```       
Rscript  script.optimization.R \
            target.gene.list \
            related.gene.list \
            expression.tsv \

```
      
### Input Description
1. *target.gene.list*   
Gene optimization are conducted to impute these target genes well. Transcriptional signatures related to target genes are going to be extracted. In example, anti-apoptotic BCL2 family genes (BCL2, MCL1, BFL1,  BCLXL, and BCLW) are used. 

2. *related.gene.list*   
Genes related with target genes. Transcriptionally-, regulationally-, or functionally-related with target genes can be included. As more genes are included, the optimization process takes exponentially longer, so we recommend less than 300 genes. In example, NanoString panel genes which selected based on our study about Acute Myeloid Leukemia (AML) data sets.

3. *expression.tsv*   
Expression matrix with row and column as genes and samples respectively. First colume name must be "Gene". Negative values are not allowed. We recommend RNA-seq or NanoString values. In example, normalized expression of NanoString of 40 AML samples are used.

### Output Description
1. *optimization.result.tsv*   

2. *optimization.MAPE.png*   

3. *H.matrix.tsv*, *W.matrix.tsv*, and *W.matrix.r1.tsv*   

