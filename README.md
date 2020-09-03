# NMFOptiGene

NMFOptiGene is a R script for finding .... (추가예정)

논문 참고 : description->논문에 언급

## Gene Optimization Example

### Script Command
```       
Rscript  script.optimization.R \
            target.gene.list \
            related.gene.list \
            expression.tsv \

```
      
### Input Description
1. target.gene.list
Gene optimization are conducted to impute these target genes well. Transcriptional signatures related to target genes are going to be extracted. In our study,for example, we use anti-apoptotic BCL2 family genes (BCL2, MCL1, BFL1,  BCLXL, and BCLW). 

2. related.gene.list
Genes related with target genes. Transcriptionally-, regulationally-, or functionally-related genes can be included. As more genes are included, the optimization process takes exponentially longer, so we recommend less than 300 genes.

3. expression.tsv
