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
library(stringr)
library(tidyr)
library(BioCircos)
library(EnhancedVolcano)
library("DESeq2")
library(hrbrthemes)
library(viridis)
library("pvca")
library("Biobase")
dev.off()

outputPrefix <- "CCH_stimulated_ncRNA_"

cts = read.csv("CCH_mirna_counts.csv", header = T)
annotation = cts[,c("Transcript.Name","biotype")]


# give sample names 
colnames_bar <- c("Transcript.Name",
                 "T00_1",
                 "T00_2",
                 "T00_3",
                 "T00_4",
                 "T00_5",
                 "T00_6",
                 "T15_1",
                 "T15_2",
                 "T15_3",
                 "T15_4",
                 "T15_5",
                 "T15_6",
                 "T30_1",
                 "T30_2",
                 "T30_3",
                 "T30_4",
                 "T30_5",
                 "T30_6",
                 "T60_1",
                 "T60_2",
                 "T60_3",
                 "T60_4",
                 "T60_5",
                 "T60_6",
                 "biotype"
                 
)

# here a copy of a dataframe will be provided
barplot_ncRNA = cts
colnames(barplot_ncRNA) = colnames_bar


# here we summarize the counts to the biotype features
library("dplyr")      
barplot_all = gather(barplot_ncRNA, Sample, Counts, 2:25) 
barplot_all <-barplot_all %>%
  group_by(Sample, biotype) %>%
  summarise(percentage = sum(Counts)) %>%
  mutate(freq = percentage/ sum(percentage))


ggplot(barplot_all, aes(x = Sample, y = freq, fill = biotype)) + 
  geom_bar(stat='identity') +
  scale_fill_viridis(discrete = T) +
  xlab("Samples")+ 
  ylab("Frequency")+
  facet_grid(. ~ "Frequency per ncRNA type")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# retrieve the annotation space for each subtype of ncRNA
#annotation_sno <- subset(annotation, biotype == "snoRNA")
#annotation_sn <- subset(annotation, biotype =="snRNA")
#annotation_trna <- subset(annotation, biotype == "tRNA")
annotation_miRNA <- subset(annotation, biotype =="miRNA")


# set the rownames and the colnames between the sample table and the count table equal
rownames(cts) = cts$Transcript.Name
cts = subset(cts, select = -c(Transcript.Name,biotype))

# give sample names 
sampleNames <- c(
                  "T00_1",
                  "T00_2",
                  "T00_3",
                  "T00_4",
                  "T00_5",
                  "T00_6",
                  "T15_1",
                  "T15_2",
                  "T15_3",
                  "T15_4",
                  "T15_5",
                  "T15_6",
                  "T30_1",
                  "T30_2",
                  "T30_3",
                  "T30_4",
                  "T30_5",
                  "T30_6",
                  "T60_1",
                  "T60_2",
                  "T60_3",
                  "T60_4",
                  "T60_5",
                  "T60_6"
                  
)

colnames(cts) <- sampleNames

time = rep(c("T00", "T15", "T30","T60"), each = 6)
sampleTable <- data.frame(sampleName = sampleNames, time = time)
rownames(sampleTable) <- colnames(cts)

# run the differential expression analysis 
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sampleTable,
                              design = ~ time)
dds <- DESeq(dds, test = "LRT", reduced=~1)
res <- results(dds)

# shrink the foldchange
res <- lfcShrink(dds=dds, coef=2, res=res, type = "ashr") 

# extract only significant genes
ressig = subset(res, padj < 0.05)


# save data results and normalized reads to csv
# insert gene_name_chromosome to the table 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata1 <- merge(as.data.frame(ressig), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)




write.csv(resdata, file = paste0(outputPrefix, "normalized_expression_counts.csv"))
write.csv(resdata1, file = paste0(outputPrefix, "significant_normalized_expression_counts.csv"))
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

#############################################################################################
library(dittoSeq)
data_pca <- plotPCA(vsd,intgroup=c("time","sizeFactor"),returnData = T)
head(data_pca)

bulk_ditto <- importDittoBulk(x = dds)
bulk_ditto <- addDimReduction(object = bulk_ditto, embeddings = data_pca, name = "pca")
dittoDimPlot(bulk_ditto,"time", size = 3, do.ellipse = T, reduction.use = data_pca)

#####miRNA only PCA
pca_mirna = plotPCA(vsd[annotation_miRNA$Transcript.Name,],intgroup=c("time","sizeFactor"),returnData = T)
head(pca_mirna)

bulk_ditto <- importDittoBulk(x = dds)
bulk_ditto <- addDimReduction(object = bulk_ditto, embeddings = pca_mirna, name = "pca")
dittoDimPlot(bulk_ditto,"time", size = 3, do.ellipse = T, reduction.use = pca_mirna)

#####miRNA only PCA
pca_trna = plotPCA(vsd[annotation_trna$Transcript.Name,],intgroup=c("time","sizeFactor"),returnData = T)
head(pca_trna)

bulk_ditto <- importDittoBulk(x = dds)
bulk_ditto <- addDimReduction(object = bulk_ditto, embeddings = pca_trna, name = "pca")
dittoDimPlot(bulk_ditto,"time", size = 3, do.ellipse = T, reduction.use = pca_trna)
###################################################################################################
resOrdered <- resdata
lncRNA_res <- resOrdered[order(resOrdered$padj),]
top_lnc <- head(lncRNA_res$Row.names,50)



#topgenes <- head(rownames(resOrdered),50)
mat <- assay(vsd)[top_lnc,] 
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("time")])

colnames(df)[1] <- "time"

rownames(df) <- colnames(mat)

df<-arrange(df,time)



pheatmap::pheatmap(mat,annotation = df,cluster_cols=T, color=colorRampPalette(c("purple", "black", "yellow"))(75),
                   scale="row",border_color=NA, cluster_rows = T)

####################################################################################################
#counts and report

counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))
degQC(counts, design[["time"]], pvalue = res[["pvalue"]])
resCov <- degCovariates(log2(counts(dds)+0.5),
                        colData(dds))



