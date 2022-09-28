# OptimalGeneNMF
[![DOI](https://zenodo.org/badge/292222008.svg)](https://zenodo.org/badge/latestdoi/292222008)

OptimalGeneNMF is an R script for optimizing genes, resulting in represent a complicated regulation network of the given "target genes". Finally, from the "optimal genes", it extracts transcriptional signatures using Non-negative Matrix Factorization (NMF). Details are described  in our paper in preparation for publication (to be updated later).

An overall scheme of our gene optimization is depicted in the below figure in which "target genes" indicate BCL2 family genes (BCL2, MCL1, BFL1, BCLXL, and BCLW) as an example. In order to select a dataset-specific "optimal genes" that represents transcriptional signatures of the "target genes", a backward selection was performed with the input gene set. As a result, a subset of input genes was selected so as to optimally impute the expression profiles of target genes. 

Here, to evaluate imputation power, we constructed 10 sets of simulated missing datasets in which "target genes" of each fold in 10-fold split samples were set to missing value (set to NA in R). After applying NMF to the each simulated missing dataset, their results were multiplied again to produce a “imputed matrix” whose missing values were imputed by NMF. In each imputed matrix, mean absolute percentage error (MAPE) of "target genes" was calculated by comparing it with those in the original matrix. We used an average of mean absolute percentage errors (aMAPE) from 10 imputed matrices to measure imputation accuracy.

Using aMAPE as a selection criterion, we performed a backward selection with the input genes. For each round, genes except for the "target genes" were individually eliminated to calculate aMAPE after removal of each gene, and the genes with minimal aMAPE were removed. This process was repeated until the number of genes reached 20, and we called the genes with minimal aMAPE “optimal genes”. Subsequently, we extracted the signatures using NMF from the expression matrix consisting of optimized genes. In outcomes of NMF, the W matrix represents the coefficient of the optimized genes in the signatures, and the H matrix represents the signatures profile of the samples.


![image](https://user-images.githubusercontent.com/70630535/177025236-b5a930f6-6832-4f08-bb00-ebf579c48063.png)




## Gene Optimization using NMF

### Script Command
```       
Rscript  script.optimization.R \
            target.gene.list \
            related.gene.list \
            expression.tsv 

```
      
### Input Description
1. *target.gene.list*   
Gene optimization is conducted to filter out the genes that do not contribute to the recovery of these target genes. As a result, transcriptional signatures that reflect the regulation factors of these target genes are extracted. In the example data, anti-apoptotic BCL2 family genes (BCL2, MCL1, BFL1,  BCLXL, and BCLW) are used. 

2. *related.gene.list*   
Genes related to regulation of target genes. Transcriptionally-, regulationally-, or functionally-related to target genes can be included. As more genes are included, the optimization process takes exponentially longer, so we recommend less than 300 genes. In example data, NanoString panel genes which selected based on our study about Acute Myeloid Leukemia (AML) data sets.

3. *expression.tsv*   
Expression matrix with row and column as genes and samples, respectively. First column name must be "Gene". Negative values are not allowed. We recommend RNA-seq or NanoString data. In the example data, the normalized expressions of NanoString of 40 AML samples are used.

### Output Description
1. *optimization.result.tsv*   
Record of the optimization process. In each step, a removed gene and Mean Absolute Percentage Error (MAPE) between original and imputed entries of target genes are recorded.

2. *optimization.MAPE.png*   
Visualization of optimization with decreasing MAPE for each step.

3. *optimized.gene.list*
Final optimized genes representing regulation network of the target genes.

4. *H.matrix.tsv*, *W.matrix.tsv*, and *W.matrix.r1.tsv*   
*H.matrix.tsv* is profile of signatures. *W.matrix.tsv* presents definition of signatures. *W.matrix.r1.tsv* is a scaled (row sum to 1) matrix of *W.matrix.tsv*.



## Citing NMFOptiGene
~~If you use NMFOptiGene for published work, please cite our paper~~   
Lee, C., Lee, S., Park, E. et al. Transcriptional signatures of the BCL2 family for individualized acute myeloid leukaemia treatment. Genome Med 14, 111 (2022). https://doi.org/10.1186/s13073-022-01115-w

## Contact
Chansub Lee (cslee159@snu.ac.kr or cslee159@naver.com)





