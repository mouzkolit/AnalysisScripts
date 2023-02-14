library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
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
library(org.Hs.eg.db)
library("biomaRt")
library(stringr)
outputPrefix = "CCH_stimulated_mRNA_"

cts <- read.csv("summary_unique.tsv.csv", header = TRUE, sep = "\t", row.names = 1)

# give sample names 
sampleNames <- c("T00_1",
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
         
                 
) # please add here the name of the samples

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



# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')


names(resdata1)[1] <- 'ensembl_gene_id'
mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",         
                            mart    = useMart("ENSEMBL_MART_ENSEMBL",       
                                              host    = "www.ensembl.org"))      
a_sig <- resdata1
b_sig <- getBM(attributes=c("ensembl_gene_id",
                            "external_gene_name",
                            "gene_biotype"
                            
),
filters = "ensembl_gene_id",
values=a_sig$ensembl_gene_id,
mart=mart)

resdata1 <- merge(a_sig, b_sig, by="ensembl_gene_id")

write.csv(resdata, file = paste0(outputPrefix, "normalized_expression_counts.csv"))
write.csv(resdata1, file = paste0(outputPrefix, "significant_normalized_expression_counts.csv"))

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

# use the mart to switch to gene symbols 

symbols <- vsd_counts
symbols$ensembl_gene_id <- rownames(symbols)

symbols_genes <- getBM(attributes=c("ensembl_gene_id",
                                    "external_gene_name",
                                    "gene_biotype"
                                    
                                    
),
filters = "ensembl_gene_id",
values=symbols$ensembl_gene_id,
mart=mart)

vsd_counts_symbols <- merge(symbols, symbols_genes, by="ensembl_gene_id")


write.csv(rld_counts, paste0(outputPrefix, "-rlog.csv"))
write.csv(vsd_counts_symbols, paste0(outputPrefix, "-vsd_with_symbols.csv"))


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

######################################################################################

#######################################################################################
#use dittoseq for analysis 
######################################################################################
library(dittoSeq)

data_pca <- plotPCA(vsd,intgroup=c("time","sizeFactor"), returnData = T)

head(data_pca)

bulk_ditto <- importDittoBulk(x = dds)
bulk_ditto <- addDimReduction(object = bulk_ditto, embeddings = data_pca, name = "pca")
dittoDimPlot(bulk_ditto,"time", size = 3, do.ellipse = T, reduction.use = data_pca)


#can thing about this

lnc_annotation <- resdata[vsd_counts_symbols$gene_biotype == "lncRNA",]$Row.names

#scRNA
pca_lncRNA = plotPCA(vsd[lnc_annotation,],returnData = T,intgroup=c("time","sizeFactor"))

proc_ditto <- importDittoBulk(x = dds)
proc_ditto <- addDimReduction(object = proc_ditto, embeddings = pca_lncRNA, name = "pca")
dittoDimPlot(proc_ditto,"time", size = 3, do.ellipse = T, reduction.use = pca_lncRNA)


######################################################################################
# make a heatmap of the top 50 genes
#######################################################################################


library("genefilter")

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



##############################################################################################

library("WGCNA")

#create the object for analysis of the expression with WGCNA
wgcna_input <- vsd_counts_symbols

high_input_genes = resdata[resdata$baseMean > 10,]$Row.names
wgcna_input_final = wgcna_input[wgcna_input$ensembl_gene_id %in% high_input_genes,]
wgcna_input_final = wgcna_input_final[-c(26,27)]
rownames(wgcna_input_final) <- wgcna_input_final$ensembl_gene_id
wgcna_input_final <- wgcna_input_final[-c(1)]
wgcna_input_final <- as.data.frame(t(wgcna_input_final))

#check all genes if good quality
gsg = goodSamplesGenes(wgcna_input_final, verbose = 3);
gsg$allOK


#get the sample clustering 
sampleTree = hclust(dist(wgcna_input_final), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


## get the soft threshold power 

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(wgcna_input_final,
                        powerVector = powers,
                        verbose = 5, 
                        networkType = "signed",
                        corFnc = "bicor",
                        blockSize = 20000,
                        corOptions = list(use = 'p', maxPOutliers = 0.1))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


enableWGCNAThreads(8)

net = blockwiseModules(wgcna_input_final, power = 7,
                       maxBlockSize = 20000,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.20,
                       numericLabels = FALSE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "CCH_neuron",
                       verbose = 3, 
                       randomSeed = 1234,
                       networkType = "signed",
                       deepSplit = 2,
                       corType = "bicor",
                       corOptions = list(use = 'p', maxPOutliers = 0.1))


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
df.module.labels <- as.data.frame(moduleLabels)



MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "CCH-01-networkConstruction-auto.RData")


#calculate the module relatinship using kME

hub.genes <- signedKME(wgcna_input_final, MEs, corFnc = "bicor")

rownames(vsd_counts_symbols) <- vsd_counts_symbols$ensembl_gene_id

final_gene_dataframe = merge(x = df.module.labels, y = vsd_counts_symbols, by=0, ,
                             all.x = TRUE)

final_hub = final_gene_dataframe[,-c(1)]

write.csv(final_hub, file = paste0(outputPrefix, "gene_expression_with_module.csv"))
write.csv(hub.genes, file = paste0(outputPrefix, "kME_modules.csv"))

hub_genes <- chooseTopHubInEachModule(wgcna_input_final, moduleColors)



#########################################################################
#load the trait data as binary data representing the time to evaluate positiv relationships
data_traits <- read.csv("trait_time.csv", header = TRUE, row.names = 1)
nGenes = ncol(wgcna_input_final);
nSamples = nrow(wgcna_input_final);


MEs0 = moduleEigengenes(wgcna_input_final, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, data_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(data_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#########################################################################
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(wgcna_input_final, power = 7);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


#########################################################################
#heatmap of genes 

eigenvector_modules = as.data.frame(t(net$MEs))
pheatmap::pheatmap(eigenvector_modules, cluster_cols = F,color=colorRampPalette(c("purple", "black", "yellow"))(75),
                   scale="row",border_color=NA)


