# OptimalGeneNMF

NMFOptiGene is an R script for optimizing genes, resulting in represent a complicated regulation network of the given "target genes". Finally, from the "optimized genes", it extracts transcriptional signatures using Non-negative Matrix Factorization (NMF). Details are described  in our paper in preparation for publication (to be updated later).

An overall scheme of our gene optimization is depicted in the below figure. In order to select an optimal set of genes that represent transcriptional signatures of the target genes, a backward selection was performed with the input gene set. In this procedure, we assumed better signatures show better performance of recovering the expression profiles of the target genes.  In this respect, we sought to minimize the recovery error of the target genes by reconstructing their profiles using NMF with simulated missing datasets. In this study, the recovery error was defined as aMAPE (an average of mean absolute percentage errors) of each fold in 10-fold split samples. For MAPE, this error was calculated from comparing values of the original matrix and that of the imputed matrix. For each fold of the dataset, all expressions of the target genes were nullified to construct the simulated missing dataset.  After applying NMF to the simulated dataset, their results were multiplied again to produce a “recovered matrix” whose missing values were imputed by NMF. Finally, the MAPE is calculated by comparing the recovered matrix to the original matrix. Using aMAPE as a selection criterion, we performed a backward selection with the input genes, as follows. For each round, genes except for the target genes were individually eliminated to calculate aMAPE after removal of each gene, and then the gene with minimal aMAPE was removed. This process was repeated until the number of genes reaches 20, and we call the genes with minimal aMAPE as “optimized genes” in this study. Eventually, we calculate H and W matrices from the expression matrix consisting of optimized genes.

![EMF적용_FigS14](https://user-images.githubusercontent.com/70630535/106703046-c71b9580-662c-11eb-8226-db2a75e20414.jpg)




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
Gene optimization is conducted to filter out the genes that do not contribute to the recovery of these target genes. As a result, transcriptional signatures that reflect the regulattion factors of these target genes are extracted. In example data, anti-apoptotic BCL2 family genes (BCL2, MCL1, BFL1,  BCLXL, and BCLW) are used. 

2. *related.gene.list*   
Genes related to regulation of target genes. Transcriptionally-, regulationally-, or functionally-related with target genes can be included. As more genes are included, the optimization process takes exponentially longer, so we recommend less than 300 genes. In example data, NanoString panel genes which selected based on our study about Acute Myeloid Leukemia (AML) data sets.

3. *expression.tsv*   
Expression matrix with row and column as genes and samples respectively. First column name must be "Gene". Negative values are not allowed. We recommend RNA-seq or NanoString data. In example data, the normalized expressions of NanoString of 40 AML samples are used.

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
(Publishing of NMFOptiGene is in progress.)

## Contact
Chansub Lee (cslee159@snu.ac.kr)





