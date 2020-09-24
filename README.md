# NMFOptiGene

NMFOptiGene is an R script for finding optimal genes to represent transciprional patterns of genes of interest (target genes) and finally extracting transcriptional signatures of the optimal genes by utilizing Non-negative Matrix Factorization (NMF), Missing Value Imputation (Xihui Lin (2018). [Fast Nonnegative Matrix Factorization and Applications to
Pattern Extraction, Deconvolution and Imputation](https://www.biorxiv.org/content/10.1101/321802v1.full). *bioRXiv*), and backward selection. Details are described  in our paper in preparation for publication (to be updated later).

## Gene Optimization using NMF based on Missing Value Imputation

### Script Command
```       
Rscript  script.optimization.R \
            target.gene.list \
            related.gene.list \
            expression.tsv 

```
      
### Input Description
1. *target.gene.list*   
Gene optimization is conducted to impute these target genes well. Transcriptional signatures related to target genes are going to be extracted. In example data, anti-apoptotic BCL2 family genes (BCL2, MCL1, BFL1,  BCLXL, and BCLW) are used. 

2. *related.gene.list*   
Genes related with target genes. Transcriptionally-, regulationally-, or functionally-related with target genes can be included. As more genes are included, the optimization process takes exponentially longer, so we recommend less than 300 genes. In example data, NanoString panel genes which selected based on our study about Acute Myeloid Leukemia (AML) data sets.

3. *expression.tsv*   
Expression matrix with row and column as genes and samples respectively. First column name must be "Gene". Negative values are not allowed. We recommend RNA-seq or NanoString data. In example data, the normalized expressions of NanoString of 40 AML samples are used.

### Output Description
1. *optimization.result.tsv*   
Record of the optimization process. In each step, a removed gene and Mean Absolute Percentage Error (MAPE) between original and imputed entries of target genes are recorded.

2. *optimization.MAPE.png*   
Visualization of optimization with decreasing MAPE for each step.

3. *H.matrix.tsv*, *W.matrix.tsv*, and *W.matrix.r1.tsv*   
*H.matrix.tsv* is profile of signatures. *W.matrix.tsv* presents definition of signatures. *W.matrix.r1.tsv* is a scaled (row sum to 1) matrix of *W.matrix.tsv*.


## Citing NMFOptiGene
~~If you use NMFOptiGene for published work, please cite our paper~~   
(Publishing of NMFOptiGene is in progress.)

## Contact
Chansub Lee (cslee159@snu.ac.kr)





