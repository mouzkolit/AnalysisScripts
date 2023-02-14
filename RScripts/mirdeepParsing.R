rm(list=ls())
library(DEGreport)
library(dplyr)
library(tibble)
library(xml2)
library(clusterProfiler)
library(DOSE)
library(data.table)
library(ggplot2)
library("RColorBrewer")
library("gplots")
library(org.Mm.eg.db)
library("biomaRt")
dev.off()
outputPrefix <- "wt_flox_gp130_mirna_"

countData = read.csv("/home/data-science/Desktop/mirna_gp130/miRNA_expressed_counts.csv")

sampleNames <- c("flox_1"
                 ,"flox_2"
                 ,"flox_3"
                 ,"flox_4"
                 ,"gp130_1"
                 ,"gp130_2"
                 ,"gp130_3"
                 ,"gp130_4"
)

sampleCondition <-c(rep("flox",4), rep("gp130",4))


sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments = c("flox","gp130")

library("DESeq2")
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

#guts
dds <- DESeq(ddscol)
res <-results(dds)

# here we contrast res_flox_sni

res <- lfcShrink(dds=dds, res=res, type = "ashr")
rld <- rlog(dds, blind=TRUE)
mat_rld  = assay(dds)
des = colData(ddsHTSeq)

# Subset the LRT results to return genes with padj < 0.05
ressig = subset(res, padj < 0.05)

# save data results and normalized reads to csv
# insert gene_name_chromosome to the table 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata1 <- merge(as.data.frame(ressig), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(dds, ylim=c(-8,8),main = "mRNAseq experiment", alpha = 0.05)
dev.copy(svg, paste0(outputPrefix, "-MAplot_initial_analysis.svg"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=F)
vsd <- varianceStabilizingTransformation(dds, blind=F)

rld_counts = as.data.frame(assay(rld))
vsd_counts = as.data.frame(assay(vsd))

write.csv(rld_counts, paste0(outputPrefix, "-rlog.csv"))
write.csv(vsd_counts, paste0(outputPrefix, "-vsd.csv"))


# save normalized values
write.table(as.data.frame(assay(rld),file = paste0(outputPrefix, "-rlog-transformed-counts.csv")))
write.table(as.data.frame(assay(vsd),file = paste0(outputPrefix, "-vst-transformed-counts.csv")))



# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(svg,paste0(outputPrefix, "-variance_stabilizing.svg"))
dev.off()

# clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
dev.copy(svg, paste0(outputPrefix, "-clustering.svg"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(10,10))
dev.off()

dev.copy(svg, paste0(outputPrefix, "-pca.svg"))
degPCA((vsd_counts), colData(dds),
       condition="condition", name="condition", shape="condition")
dev.off()



# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
# plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,6:5],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

# heatmap of data
library("RColorBrewer")
library("gplots")
# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.3, labRow=F,
          main="1000 Top Expressed mRNAs")
dev.copy(svg, paste0(outputPrefix, "-HEATMAP.svg"))
dev.off()


library(EnhancedVolcano)

# make a volcano plot for the genes that changed between the conditions
EnhancedVolcano(resdata,
                lab = resdata$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1.5,
                pointSize = 3,
                title = "SNI vs. gp130-/- SNI",
                xlim = c(-2, 4),
                legendPosition = "right",
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = "grey30"
)

write.csv(resdata1, file = paste0(outputPrefix, "-results-with-normalized_significant.csv"))



write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))