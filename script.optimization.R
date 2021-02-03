library(NNLM)
library(data.table)
library(ggplot2)
library(parallel)

seed <- 777
rnk <- 3
fold <- 10
cost <- 'mkl'
n.cores <- 30
output.path <- 'output'
output.file <- 'output/optimization.result.tsv'
output.file2 <- 'output/optimization.MAPE.png'
output.file3 <- 'output/H.matrix.tsv'
output.file4 <- 'output/W.matrix.tsv'
output.file5 <- 'output/W.matrix.r1.tsv'
output.file6 <- 'output/optimized.gene.list'

# Inputs
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
	target.list.file <- args[1]
	gene.list.file <- args[2]
	expr.file <- args[3]
} else {stop('invalid inputs')}
#=============================================================================================#
## STEP1 :  Optimizing "Relateed Genes" (Backward Selection)
cat('\n## Gene Optimization...\n')

target.list <- read.csv(target.list.file, header=F, stringsAsFactors=F)$V1
gene.list <- read.csv(gene.list.file, header=F, stringsAsFactors=F)$V1
expr <- fread(expr.file, data.table=F)
try(dir.create(output.path))

gene.list <- gene.list[gene.list %in% expr$Gene]
if(any(!target.list %in% gene.list)) {stop('some of "target genes" are not exist in "related genes"...')}
row.names(expr) <- expr$Gene
expr <- as.matrix(expr[gene.list, 2:ncol(expr)])

# Start Gene Set
total.gene.list <- sort(setdiff(gene.list, target.list))
now.exist.gene.list <- total.gene.list
n.now.exist.gene <- length(now.exist.gene.list)

# Split into "Target Genes" and "Rest Genes"
expr.part1 <- expr[target.list,]
expr.part2 <- expr[now.exist.gene.list,]

# Missing for Target Genes
n.sample <- ncol(expr)
fold.n.sampe <- floor(n.sample / fold)
set.seed(seed)
random.sample <- sample(1:(fold.n.sampe * fold), fold.n.sampe * fold)
na.idx.list <- lapply(1:fold,
	function(i){
		idx <- random.sample[(fold.n.sampe*(i-1)+1) : (fold.n.sampe*(i))]
		return(idx)
})
expr.part1.mask.list <- lapply(na.idx.list,
	function(idx){
		missing.expr.part1 <- expr.part1
		missing.expr.part1[,idx] <- NA
		return(missing.expr.part1)
})

# Function : get MAPEs of target genes about given masked expression matrix
get.mape <- function(expr.part2.sub.1gn){
	cand.gene.target.mape <- lapply(target.list, function(x){c()})
	names(cand.gene.target.mape) <- target.list

        for (i in 1:fold){
		# masked index (sample index)
		na.idx <- na.idx.list[[i]]
		# masked Target Genes expression (sample index)
		expr.part1.mask <- expr.part1.mask.list[[i]]
		expr.sub.1gn.mask <- rbind(expr.part1.mask, expr.part2.sub.1gn)

		# sort by genes
		expr.sub.1gn.mask <-  expr.sub.1gn.mask[order(rownames(expr.sub.1gn.mask)),]

		# NMF
		set.seed(seed)
		z <- try(nnmf(expr.sub.1gn.mask, rnk, loss=cost, max.iter=5000, verbose = FALSE, check.k=F), T)
		if (class(z)=='try-error') {cat('NMF error!!!\n');next}
		A2 <- (with(z, W%*%H))
		A2.part1 <- A2[target.list,]
		A1.part1 <- expr.part1

		# MAPE
		for (g in target.list){
			g.mape <- mean(
				abs( (A1.part1[g,na.idx] - A2.part1[g,na.idx]) / (A1.part1[g,na.idx]) )
				) * 100
			cand.gene.target.mape[[g]] <- c(cand.gene.target.mape[[g]], g.mape)
		}
	}
	# sumary (merge fold's MAPE)
	cand.gene.target.mape.mean <- lapply(cand.gene.target.mape, mean)
	cand.gene.mape.mean <- mean(unlist(cand.gene.target.mape.mean))

	return(list(cand.gene.mape.mean, cand.gene.target.mape.mean))
}

# Function : find a gene which least helpful to impute Target Genes
backward <- function(now.exist.gene.list){
	expr.part2.sub = expr.part2[now.exist.gene.list, ]

	# get MAPEs about each remain genes
	mape.res.list <- mclapply(now.exist.gene.list, do.gene, expr.part2.sub, mc.cores=n.cores)
	mape.mean.list <- lapply(mape.res.list, function(x){x[[1]]})
	target.mape.mean.list <- lapply(mape.res.list, function(x){x[[2]]})

	# summary
	trash.gene.idx <- which.min(unlist(mape.mean.list))
	trash.gene <- now.exist.gene.list[trash.gene.idx]
	trash.gene.mape <- mape.mean.list[[trash.gene.idx]]
	trash.gene.bcl2Fam.mape <- target.mape.mean.list[[trash.gene.idx]]
	return(list(trash.gene, trash.gene.mape, trash.gene.bcl2Fam.mape))
}

# Function : get MAPE about given gene
do.gene <- function(cand.gene, expr.part2.sub){
	n.now.gene <- nrow(expr.part2.sub)
	g <- sort(setdiff(rownames(expr.part2.sub), cand.gene))
	expr.part2.sub.1gn <- expr.part2[g,]
	mape.res <- get.mape(expr.part2.sub.1gn)
	return(mape.res)
}

# Main Code
steps <- c()
genes <- c()
mapes <- c()
now.exist.gene.list.pastes <- c()
n.now.exist.genes <- c()
target.mapes <- lapply(target.list, function(x){c()})
names(target.mapes) <- target.list

step <- 0
while(n.now.exist.gene > 20){
	res <- backward(now.exist.gene.list)
	gene <- res[[1]]
	mape <- res[[2]]
	text = paste0("  step ", step, ": ",gene," removed (MAPE=", round(mape,2) ,  ")\n")
	cat(text)
	target.mape <- res[[3]]
	now.exist.gene.list <- sort(setdiff(now.exist.gene.list,gene))
	
	# update
	now.exist.gene.list.paste <- paste0(now.exist.gene.list, collapse=',')
	n.now.exist.gene <- length(now.exist.gene.list)

	genes <- c(genes, gene)
	mapes <- c(mapes, mape)
	now.exist.gene.list.pastes <- c(now.exist.gene.list.pastes, now.exist.gene.list.paste)
	n.now.exist.genes <- c(n.now.exist.genes, n.now.exist.gene)
	steps <- c(steps, step)
	for (g in target.list){
		target.mapes[[g]] <- c(target.mapes[[g]], target.mape[[g]])
	}
        step <- step + 1

        # save
	res1 <- data.frame(step=steps, removed_gene=genes,MAPE=mapes)
	res3 <- data.frame(n_gene=n.now.exist.genes, gene_list=now.exist.gene.list.pastes)
	res2 <- as.data.frame(matrix(unlist(target.mapes), nrow=length(steps)))
	colnames(res2) <- target.list
	res <- cbind(res1, res2, res3)

	fwrite(res, output.file, sep='\t')
}

#=============================================================================================#
## STEP2 :  Plotting MAPE
cat('\n## Plotting MAPE...\n')

best.n_gene <- res[which.min(res$MAPE),]$n_gene
#result2 <- melt(result, id.vars=c('step','gene','n_gene'), )
gg <- ggplot(res)
#gg <- gg + geom_line(aes(x=n_gene, y=mape), size=0.7)
gg <- gg + geom_line(aes(x=n_gene, y=MAPE), size=1)
#gg <- gg + geom_vline(xintercept=best.n_gene, linetype='dashed', color='red2', size=1)
gg <- gg + geom_vline(xintercept=best.n_gene, linetype=2, color='red', size=1)
#gg <- gg + theme_minimal()
gg <- gg + theme_classic()
gg <- gg + xlab('GeneSet Size') + ylab('Average MAPE')
gg <- gg + ggtitle('Gene Optimization')
gg <- gg + theme(text=element_text(size=20), axis.text.y=element_text(size=18,colour='black'), axis.text.x=element_text(size=18,colour='black'))
ggsave(output.file2)

#=============================================================================================#
## STEP3 :  Calculating Signatures
cat('\n## Calculating Signatures...\n')
minMAPE.gene.list <- as.character(res[which.min(res$mape),]$gene_list)
minMAPE.gene.list <- unlist(strsplit(minMAPE.gene.list, ','))

expr.sub <- expr[sort(c(target.list, minMAPE.gene.list)),]

set.seed(seed)
z <- try(nnmf(expr.sub, rnk, loss=cost, max.iter=5000, verbose = FALSE, check.k=F), T)
w <- z$W
h <- z$H

rownames(h) <- paste0('MetaGene_',1:rnk)
colnames(w) <- paste0('MetaGene_',1:rnk)

#scaling (row sum to 1)
rsum <- apply(w, 1, sum)
w2 <- lapply(1: nrow(w),
		function(i){
			res <- w[i,]/rsum[i]
			return(res)
		})
w2 <- as.matrix(do.call(rbind, w2))
rownames(w2) <- rownames(w)

fwrite(as.data.frame(h), output.file3, sep='\t', row.names=T)
fwrite(as.data.frame(w), output.file4, sep='\t', row.names=T)
fwrite(as.data.frame(w2), output.file5, sep='\t', row.names=T)
write.table(c(target.list, minMAPE.gene.list), output.file6, col.names=F, row.names=F, quote=F)

cat('\nDone...\n')


