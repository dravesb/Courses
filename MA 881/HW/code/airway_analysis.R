#-----------------------------
# Load Data + Libraries
#-----------------------------
library("airway")
data("airway")
library("DESeq2")
library("edgeR")
library("limma")

#-----------------------------
# Format Plots
#-----------------------------
par(mfrow = c(3,2))

#-----------------------------
# DESeq2
#-----------------------------

#preprocessings in DESeq lang.
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- estimateSizeFactors(dds)
dds$dex <- relevel(dds$dex, "untrt")

#run pipeline
dds <- DESeq(dds)
res <- results(dds)

#get pvalue plot
hist(res$pvalue[res$baseMean > 1], breaks=20, col="lightgray", border="blue", main = "DESeq2 P-values")

#map to normal 
hist(qnorm(1 - (res$pvalue[res$baseMean>1]/2)), breaks=20, col="lightgray", border="blue", main = "DESeq2 Z-values")
#-----------------------------
# edgeR
#-----------------------------
countdata <- assay(airway)
coldata <- colData(airway)
genetable <- data.frame(gene.id=rownames(airway))
dge <- DGEList(counts=countdata, samples=coldata, genes=genetable)

dge <- dge[rowSums(dge$counts)>1, ]

design <- model.matrix(~cell+dex, dge$samples)

dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)

#negative binomial model
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef=ncol(design))
tt <- topTags(lrt, n=nrow(dge), p.value=0.1)

#plot pvalues
hist(tt$table$PValue, breaks=20, col="lightgray", border="blue", main= "edgeR P-values")

#map to normal 
hist(qnorm(1 - (tt$table$PValue/2)), breaks=20, col="lightgray", border="blue", main = "edgeR Z-values")

#-----------------------------
# Limma + Voom
#-----------------------------
design <- model.matrix(~cell+dex, data=dge$samples)
colnames(design)

colnames(design) <- gsub("cell", "", colnames(design))
colnames(design)[1] <- "Intercept"

v <- voom(dge, design, plot=FALSE)
vfit <- lmFit(v, design)
efit <- eBayes(vfit)

#plot pvalues
hist(efit$p.value[,5], breaks=20, col="lightgray", border="blue", main= "Limma+Voom P-Values")

#map to normal 
hist(qnorm(1 - (efit$p.value[,5]/2)), breaks=20, col="lightgray", border="blue", main = "Limma+Voom Z-values")





