rm(list=ls())
library(dplyr)
library(tibble)
library(xml2)
library(clusterProfiler)
library(data.table)
library(ggplot2)
library("RColorBrewer")
library("gplots")
library(stringr)
library(tidyr)
library(BioCircos)
library(EnhancedVolcano)
library("DESeq2")
library(hrbrthemes)
library(viridis)
library("pvca")
library("Biobase")
library("rtracklayer")
library("dplyr")
library("comprehenr")
library("genefilter")
library("WGCNA")
dev.off()

##################################################################
#Load the Data
##################################################################
outputPrefix <- "excerpt_RNA_merged_trna"

## read matrix of the ncRNAs
# read the count fle for ncRNAs
object_excerpt <- load("/home/data-science/Manatee/iPSC_excerpt/FinalResult/exceRpt_smallRNAQuants_ReadCounts.RData")

## make the barplot for the analysis
barplot_ncRNA = read.table("/home/data-science/Manatee/iPSC_excerpt/FinalResult/exceRpt_biotypeCounts.txt")
barplot_ncRNA$biotype = rownames(barplot_ncRNA)


##################################################################
#Preprocess Data
##################################################################
vector_biotypes = c("miRNA","piRNA","tRNA","rRNA","snoRNA","snRNA", "lincRNA","protein_coding","misc_RNA")

# make the vector for the new biotypes
empty_vector = c()
for(i in barplot_ncRNA$biotype){
  if(i %in% vector_biotypes){
    empty_vector <- append(empty_vector, i)
  }
  else{
    empty_vector <- append(empty_vector, "other")
  }
}

colnames_bar = c( "ReferenceID",
                  "D00_1",
                  "D00_2",
                  "D00_3",
                  "D05_1",
                  "D05_2",
                  "D05_3",
                  "D09_1",
                  "D09_2",
                  "D09_3",
                  "D16_1",
                  "D16_2",
                  "D16_3",
                  "D26_1",
                  "D26_2",
                  "D26_3",
                  "D36_1",
                  "D36_2",
                  "D36_3",
                  "D00_4",
                  "D00_5",
                  "D00_6",
                  "D05_4",
                  "D05_5",
                  "D05_6",
                  "D09_4",
                  "D09_5",
                  "D09_6",
                  "D16_4",
                  "D16_5",
                  "D16_6",
                  "D26_4",
                  "D26_5",
                  "D26_6",
                  "D36_4",
                  "D36_5",
                  "D36_6",
                  "D00_7",
                  "D00_8",
                  "D00_9",
                  "D05_7",
                  "D05_8",
                  "D05_9",
                  "D09_7",
                  "D09_8",
                  "D09_9",
                  "D16_7",
                  "D16_8",
                  "D16_9",
                  "D26_7",
                  "D26_8",
                  "D26_9",
                  "D36_7",
                  "D36_8",
                  "D36_9",
                  "biotype"
)


###############################################################
# Basic Exploration for ncRNAs all
###############################################################

barplot_ncRNA$biotype = empty_vector
colnames(barplot_ncRNA) = colnames_bar[2:56]
# make an aggregation of the plots to get the time course
barplot_agg = gather(barplot_ncRNA, Sample, Counts, 1:54) 


#make a barplot to show the percentages of distribution 
barplot_bar = barplot_agg %>%
  group_by(Sample, biotype) %>%
  summarise(percentage = sum(Counts)) %>%
  mutate(freq = percentage/ sum(percentage))


# creating the ggplot object
group.colors <- c(miRNA = "#262424", 
                  lincRNA = "#eb1b0c",
                  rRNA ="#4a3d3d",
                  snoRNA = "#dcde57",
                  snRNA = "#3dd12a",
                  misc_RNA = "blue",
                  tRNA = "#aa34cf",
                  protein_coding = "#878787",
                  other = "lightgreen",
                  piRNA = "skyblue")

ggplot(barplot_bar, aes(x = Sample, y = freq, fill = biotype)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values=group.colors) +
  xlab("Samples")+ 
  ylab("Frequency")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(. ~ "Frequency per ncRNA type")

#make a pie chart for the different ncRNA types aggregated on the Day
pieplot_chart <- barplot_bar
pieplot_chart$time <- str_split_fixed(pieplot_chart$Sample,"_",2)[,1]

# group by the biotype and per timepoint for the creation of the pie-chart
pieplot_chart_trial <- pieplot_chart %>% 
  dplyr::group_by(time, biotype) %>% 
  dplyr::summarise(pie = sum(percentage)) %>%
  dplyr::mutate(freq = round(100* pie/ sum(pie),2))

# define the positions of the label within the pie chart
df2 <- pieplot_chart_trial %>% 
  mutate(csum = rev(cumsum(rev(freq))), 
         pos = freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), freq/2, pos))


#drawing the pie charts with custom 
ggplot(data=as.data.frame(pieplot_chart_trial), aes(x="", y=freq, color = biotype, fill = biotype)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  coord_polar("y") + 
  geom_label_repel(data = df2,
                  aes(y = pos, label = paste0(freq, "%")),
                  size = 3, nudge_x = 1, show.legend = FALSE,
                  color = "white"
                  ) +
  guides(fill = guide_legend(title = "biotype")) +
  facet_wrap(.~ time, ncol = 3)+
  theme(legend.position="bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        strip.background = element_rect(colour="black", fill = "grey"))

#############################################################
#Preprocess tRNA Expresssion
#############################################################

trna_expression <- read.csv("tRNA_files_expression.csv")

# reduce heterogenity by merging and aggregating counts belonging to the
# same isoacceptor
rownames(trna_expression) = trna_expression$ReferenceID
trna_expression$trna_acceptor = str_split_fixed(rownames(trna_expression), "-(?=[^-]+$)",2)[,1] 
trna_expression = trna_expression %>% group_by(trna_acceptor) %>% summarise_if(is.numeric,sum)
trna_expression = as.data.frame(trna_expression)

# make a pie chart for the percentage wise expression of 
trna_pie <- trna_expression
trna_pie$trna_acceptor <- str_split_fixed(trna_pie$trna_acceptor, "-(?=[^-]+$)",3)[,1] 
trna_pie$trna_acceptor <- str_split_fixed(trna_pie$trna_acceptor, "-(?=[^-]+$)",3)[,1] 
colnames(trna_pie) <- colnames_bar[1:55]
pie_trna <- gather(trna_pie, Sample, Counts, 2:55) 
pie_trna$time <- str_split_fixed(pie_trna$Sample, "_",2)[,1]

###########################################################
# pie chart table tRNA
###########################################################

# merge and groupby data to obtain frequency per timepoint
pieplot_trna <- pie_trna %>% 
  dplyr::group_by(time, ReferenceID) %>% 
  arrange(desc(Counts)) %>%
  dplyr::summarise(pie = sum(Counts)) %>%
  filter(across(everything(), ~ !grepl("nm", .))) %>%
  dplyr::mutate(freq = round(100* pie/ sum(pie),2)) 
 
# define the positions of the label within the pie chart
position_trna <- pieplot_trna %>% 
  mutate(csum = rev(cumsum(rev(freq))), 
         pos = freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), freq/2, pos)) 

# only select the ones that represent a high number of counts
position_trna <- position_trna[position_trna$freq > 3,]

# make a personal color map
colourCount = 24
getPalette = colorRampPalette(brewer.pal(9, "Paired"))

#drawing the pie charts with custom 
ggplot(data=as.data.frame(pieplot_trna), aes(x="", y=freq, color = ReferenceID, fill = ReferenceID)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") + 
  geom_label_repel(data = position_trna,
                   aes(y = pos, label = paste0(ReferenceID, " ",freq, "%")),
                   size = 3, nudge_x = 1, show.legend = FALSE,
                   color = "white"
  ) +
  guides(fill = guide_legend(title = "ReferenceID")) +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_color_manual(values = getPalette(colourCount)) +
  facet_wrap(.~ time, ncol = 3)+
  theme(legend.position="bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        strip.background = element_rect(colour="black", fill = "grey"))

#group by the acceptor
trna_expression$biotype = "tRNA"
rownames(trna_expression) = trna_expression$trna_acceptor
colnames(trna_expression) = colnames_bar
trna_expression = trna_expression[,2:56]

###############################################################
# merge and concat different ncRNA group tables for DGE
###############################################################

biotype_gen <- as.data.frame(as.matrix(str_split_fixed(rownames(exprs.gencode),":",2)))
exprs.gencode <- as.data.frame(exprs.gencode)
exprs.miRNA <- as.data.frame(exprs.miRNA)
exprs.piRNA <- as.data.frame(exprs.piRNA)
exprs.piRNA$biotype = "piRNA"
exprs.miRNA$biotype = "miRNA"
exprs.gencode$biotype = biotype_gen$V2
exprs.gencode = exprs.gencode[exprs.gencode$biotype != "miRNA",]
colnames(exprs.piRNA) = colnames_bar[2:56]
colnames(exprs.miRNA) = colnames_bar[2:56]
colnames(exprs.gencode)= colnames_bar[2:56]

# get the annotation data
cts <- rbind(exprs.gencode,exprs.miRNA, exprs.piRNA, trna_expression)

# make an annotation data table
annotation = cts
annotation$mean = rowMeans(annotation[1:54], na.rm = FALSE)
annotation =  annotation[annotation$mean > 1,]
annotation$ReferenceId = rownames(annotation)
annotation = annotation[,c("ReferenceId","biotype")] 

# retrieve the annotation space for each subtype of ncRNA
annotation_sno <- subset(annotation, biotype == "snoRNA")
annotation_sn <- subset(annotation, biotype =="snRNA")
annotation_trna <- subset(annotation, biotype == "tRNA")
annotation_rrna <- subset(annotation, biotype == "rRNA")
annotation_lincRNA <- subset(annotation, biotype == "lincRNA")
annotation_proc <- subset(annotation, biotype == "processed_transcript")
annotation_vault <- subset(annotation, biotype == "vaultRNA")
annotation_scRNA <- subset(annotation, biotype =="scaRNA")
annotation_miRNA <- subset(annotation, biotype =="miRNA")
annotation_piRNA <- subset(annotation, biotype =="piRNA")
annotation_mttrna <- subset(annotation, biotype =="Mt_tRNA")
annotation_protein_coding <- subset(annotation, biotype =="protein_coding")


annotation_ncrnas = subset(annotation, biotype != "protein_coding")
annotation_ncrnas = subset(annotation_ncrnas, biotype %in% c(
                                                             "snoRNA",
                                                             "tRNA",
                                                             "rRNA",
                                                             "vaultRNA",
                                                             "lincRNA",
                                                             "piRNA",
                                                             "Mt_tRNA",
                                                             "scaRNA",
                                                             "snRNA"))
# remove the biotype again for differential gene expression analysis
cts = subset(cts, select = -c(biotype))
cts = cts[rownames(cts) %in% annotation$ReferenceId,]
cts <- mutate_all(cts, function(x) as.integer(as.character(x)))

###################################################################
# construct a metadata table for the analysis
###################################################################

cell_line = rep(c("AD2", "840", "AD3"), each = 18)
time = rep(c("DAY00", "DAY05", "DAY09", "DAY16", "DAY26","DAY36"), each = 3, times = 3)
sampleTable <- data.frame(sampleName = colnames(cts), condition = time, cellline = cell_line)

#build the deseq2 object iwth metadata and design
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sampleTable,
                              design = ~ cellline + condition + condition:cellline)

#run the differential expression using the LRT test
dds <- DESeq(dds, test = "LRT", reduced=~cellline)
res <-results(dds)
res <- lfcShrink(dds=dds, res=res, type = "ashr")
# Subset the LRT results to return genes with padj < 0.05
ressig = subset(res, padj < 0.01)


# save data results and normalized reads to csv
# insert gene_name_chromosome to the table 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata1 <- merge(as.data.frame(ressig), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
basemean <- resdata[resdata1$baseMean > 20,]

# save intermediate output
write.csv(resdata, file = paste0(outputPrefix, "normalized_expression_counts.csv"))
write.csv(resdata1, file = paste0(outputPrefix, "significant_normalized_expression_counts.csv"))

#########################################################
# visualization Analysis 
#########################################################

plotMA(dds, ylim=c(-8,8),main = "ncRNA expression", alpha = 0.05)
dev.copy(svg, paste0(outputPrefix, "-MAplot_initial_analysis.svg"))
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
vsd <- varianceStabilizingTransformation(dds, blind=F)
vsd_counts = as.data.frame(assay(vsd))
#vsd = vsd[rownames(vsd) %in% basemean$Row.names]
write.csv(vsd_counts, paste0(outputPrefix, "-vsd.csv"))
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

#use dittoseq for analysis 
######################################################################################
# Principle Component Analysis splitted on the ncRNA biotypes
######################################################################################


data_pca <- plotPCA(vsd, returnData = T)
data_pca$celline <- cell_line
percentVar <- round(100 * attr(data_pca, "percentVar"))
head(data_pca)

# snoRNAs only
pca_sno = plotPCA(vsd[annotation_sno$ReferenceId,],returnData = T)
pca_sno$celline <- cell_line
percentVar_sno <- round(100 * attr(pca_sno, "percentVar"))
pca_sno

# snRNA PCA
pca_sn = plotPCA(vsd[annotation_sn$ReferenceId,],returnData = T)
pca_sn$celline <- cell_line
percentVar_sn <- round(100 * attr(pca_sn, "percentVar"))
head(pca_sn)

# rRNA PCA
pca_rrna = plotPCA(vsd[annotation_rrna$ReferenceId,],returnData = T)
pca_rrna$celline <- cell_line
percentVar_rrna <- round(100 * attr(pca_rrna, "percentVar"))
head(pca_rrna)

#lincRNA PCA
pca_linc = plotPCA(vsd[annotation_lincRNA$ReferenceId,],returnData = T)
pca_linc$celline <- cell_line
percentVar_linc <- round(100 * attr(pca_linc, "percentVar"))
head(pca_linc)

#processed_trnascripts
pca_prot = plotPCA(vsd[annotation_protein_coding$ReferenceId,],returnData = T)
pca_prot$celline <- cell_line
percentVar_prot <- round(100 * attr(pca_prot, "percentVar"))
head(pca_prot)

#scRNA
pca_scRNA = plotPCA(vsd[annotation_scRNA$ReferenceId,],returnData = T)
pca_scRNA$celline <- cell_line
percentVar_scRNA <- round(100 * attr(pca_scRNA, "percentVar"))
head(pca_scRNA)

#piRNA
pca_piRNA = plotPCA(vsd[annotation_piRNA$ReferenceId,],returnData = T)
pca_piRNA$celline <- cell_line
percentVar_pirna <- round(100 * attr(pca_piRNA, "percentVar"))
head(pca_piRNA)

#tRNA
pca_tRNA = plotPCA(vsd[annotation_trna$ReferenceId,],returnData = T)
pca_tRNA$celline <- cell_line
percentVar_trna <- round(100 * attr(pca_tRNA, "percentVar"))
head(pca_tRNA)

ggplot_pca <- function(data, percentVar){
  p1 <- ggplot(data, aes(PC1, PC2, color = condition))+
    geom_point(aes(shape = celline), size=3)+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))+
    stat_ellipse()+
    facet_grid(. ~ "DESeq2 ncRNA distribution")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p1)
}

# pca plot for the analysis
all_nc <- ggplot_pca(data_pca, percentVar)
sno_nc <- ggplot_pca(pca_sno, percentVar_sno)
sn_nc <- ggplot_pca(pca_sn, percentVar_sn)
rr_nc <- ggplot_pca(pca_rrna, percentVar_rrna)
trna_nc <- ggplot_pca(pca_tRNA, percentVar_trna)
pirna_nc <- ggplot_pca(pca_piRNA, percentVar_pirna)
prot_nc <- ggplot_pca(pca_prot, percentVar_prot)


######################################################################################
# make a heatmap of the top 50 genes
######################################################################################
resOrdered <- res[(rownames(res) %in% annotation_ncrnas$ReferenceId),]
snoRNA_res <- resOrdered[order(resOrdered$padj),]
top_sno <- head(rownames(snoRNA_res),50)
#topgenes <- head(rownames(resOrdered),50)
vsd <- vsd[, order(colnames(vsd))]
mat <- assay(vsd)[top_sno,] 
mat <- mat - rowMeans(mat)
dds_trial = dds[, order(colnames(dds))]
df <- as.data.frame(colData(dds))
df <- df[order(rownames(df)),]
df <- df[,c("sampleName","condition")]
colnames(df)[2] <- "time"
rownames(df) <- colnames(mat)
df<-arrange(df,time)
df<- subset(df, select = time)


pheatmap::pheatmap(mat,
                   annotation = df,
                   cluster_cols=F,
                   clusterin_distance_rows = "correlation",
                   annotation_legend = F,
                   color = hcl.colors(50, "Lisbon"),
                   scale="row",
                   border_color="grey",
                   show_colnames = F,
                   fontsize_number = 7,
                   gaps_col = c(9,18,27,36,45),
                   cellwidth =10,
                   cellheight = 10,
                   filename = "all_nc_rna_expression.pdf"
                   )


###############################################################################
# Select columns sort on Columns for further downstream analysis and retrieve the vsd files and padj values
###############################################################################

vsd_significant <- vsd_counts[resdata1$Row.names,] # get only the significant ones
column_vector <- c("D00_1",
                   "D00_2",
                   "D00_3",
                   "D00_4",
                   "D00_5",
                   "D00_6",
                   "D00_7",
                   "D00_8",
                   "D00_9",
                   "D05_1",
                   "D05_2",
                   "D05_3",
                   "D05_4",
                   "D05_5",
                   "D05_6",
                   "D05_7",
                   "D05_8",
                   "D05_9",
                   "D09_1",
                   "D09_2",
                   "D09_3",
                   "D09_4",
                   "D09_5",
                   "D09_6",
                   "D09_7",
                   "D09_8",
                   "D09_9",
                   "D16_1",
                   "D16_2",
                   "D16_3",
                   "D16_4",
                   "D16_5",
                   "D16_6",
                   "D16_7",
                   "D16_8",
                   "D16_9",
                   "D26_1",
                   "D26_2",
                   "D26_3",
                   "D26_4",
                   "D26_5",
                   "D26_6",
                   "D26_7",
                   "D26_8",
                   "D26_9",
                   "D36_1",
                   "D36_2",
                   "D36_3",
                   "D36_4",
                   "D36_5",
                   "D36_6",
                   "D36_7",
                   "D36_8",
                   "D36_9")
vsd_significant <- vsd_significant[,column_vector]

resdata1 <- resdata1[ , c("Row.names",
                          "baseMean",
                          "log2FoldChange",
                          "lfcSE",
                          "pvalue",
                          "padj",
                          "D00_1",
                          "D00_2",
                          "D00_3",
                          "D00_4",
                          "D00_5",
                          "D00_6",
                          "D00_7",
                          "D00_8",
                          "D00_9",
                          "D05_1",
                          "D05_2",
                          "D05_3",
                          "D05_4",
                          "D05_5",
                          "D05_6",
                          "D05_7",
                          "D05_8",
                          "D05_9",
                          "D09_1",
                          "D09_2",
                          "D09_3",
                          "D09_4",
                          "D09_5",
                          "D09_6",
                          "D09_7",
                          "D09_8",
                          "D09_9",
                          "D16_1",
                          "D16_2",
                          "D16_3",
                          "D16_4",
                          "D16_5",
                          "D16_6",
                          "D16_7",
                          "D16_8",
                          "D16_9",
                          "D26_1",
                          "D26_2",
                          "D26_3",
                          "D26_4",
                          "D26_5",
                          "D26_6",
                          "D26_7",
                          "D26_8",
                          "D26_9",
                          "D36_1",
                          "D36_2",
                          "D36_3",
                          "D36_4",
                          "D36_5",
                          "D36_6",
                          "D36_7",
                          "D36_8",
                          "D36_9"
)]

# retrieve the biotype for the significant normalized counts
annotation_sno_significant <- resdata1[resdata1$Row.names %in% annotation_sno$ReferenceId,]
annotation_sn_significant <- resdata1[resdata1$Row.names %in% annotation_sn$ReferenceId,]
annotation_trna_significant <- resdata1[resdata1$Row.names %in% annotation_trna$ReferenceId,]
annotation_rrna_significant <- resdata1[resdata1$Row.names %in% annotation_rrna$ReferenceId,]
annotation_proc_significant <- resdata1[resdata1$Row.names %in% annotation_proc$ReferenceId,]
annotation_linc_significant <- resdata1[resdata1$Row.names %in% annotation_lincRNA$ReferenceId,]
annotation_mirna_significant <- resdata1[resdata1$Row.names %in% annotation_miRNA$ReferenceId,]
annotation_pirna_significant <- resdata1[resdata1$Row.names %in% annotation_piRNA$ReferenceId,]
annotation_prot_significant <- resdata1[resdata1$Row.names %in% annotation_protein_coding$ReferenceId,]

#save the dataframe
write.csv(annotation_sno_significant,"snoRNA_normalized_counts.csv")
write.csv(annotation_sn_significant,"snRNA_normalized_counts.csv")
write.csv(annotation_trna_significant,"tRNA_normalized_counts.csv")
write.csv(annotation_rrna_significant,"rRNA_normalized_counts.csv")
write.csv(annotation_proc_significant,"processed_normalized_counts.csv")
write.csv(annotation_linc_significant,"linc_normalized_counts.csv")
write.csv(annotation_pirna_significant,"pirna_normalized_counts.csv")
write.csv(annotation_prot_significant,"protein_normalized_counts.csv")
#get the variance stabilized counts for the biotype
vsd_sno <- vsd_significant[annotation_sno_significant$Row.names,]
vsd_snRNA <- vsd_significant[annotation_sn_significant$Row.names,]
vsd_tRNA <- vsd_significant[annotation_trna_significant$Row.names,]
vsd_rRNA <- vsd_significant[annotation_rrna_significant$Row.names,]
vsd_processed <- vsd_significant[annotation_proc_significant$Row.names,]
vsd_linc <- vsd_significant[annotation_linc_significant$Row.names,]
vsd_pirna <- vsd_significant[annotation_pirna_significant$Row.names,]

# save the plots
write.csv(vsd_sno,"snoRNA_vsd_counts.csv")
write.csv(vsd_snRNA,"snRNA_vsd_counts.csv")
write.csv(vsd_tRNA,"tRNA_vsd_counts.csv")
write.csv(vsd_rRNA,"rRNA_vsd_counts.csv")
write.csv(vsd_processed,"processed_vsd_counts.csv")
write.csv(vsd_linc,"linc_vsd_counts.csv")
write.csv(vsd_pirna,"pirna_vsd_counts.csv")

#save the datafile
save.image(file= "excerpt_ncRNA_workspace.RData")

###########################################################
# WGCNA Analysis to obtain networks of correlated ncRNA structures
###########################################################


#create the object for analysis of the expression with WGCNA
wgcna_input <- vsd_counts
wgcna_input <- vsd_counts[!rownames(vsd_counts) %in% annotation_miRNA$ReferenceId,]
wgcna_input <- wgcna_input[!rownames(wgcna_input) %in% annotation_protein_coding$ReferenceId,]
#wgcna_input <- wgcna_input[rownames(wgcna_input) %in% annotation_nc_only$ReferenceId,]
high_input_genes = resdata[resdata$baseMean > 20,]$Row.names
wgcna_input_final = wgcna_input[rownames(wgcna_input) %in% high_input_genes,]
#wgcna_input_final = wgcna_input_final[-c(26,27)]
#rownames(wgcna_input_final) <- wgcna_input_final$ensembl_gene_id
#wgcna_input_final <- wgcna_input_final[-c(1)]

wgcna_input_final <- as.data.frame(t(wgcna_input_final))
# 
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
plot(sampleTree, main = "Sample clustering (w/o miRNA) to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


## get the soft threshold power 

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(wgcna_input_final,
                        powerVector = powers,
                        verbose = 5, 
                        networkType = "signed",
                        blockSize = 20000
)
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

# Run the WGCNA blockwise module fucntion to detect the modules
enableWGCNAThreads(8)

net = blockwiseModules(wgcna_input_final, power = 16,
                       maxBlockSize = 20000,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = FALSE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ncRNA_without_mirna_module",
                       verbose = 3, 
                       randomSeed = 1234,
                       networkType = "signed",
                       deepSplit = 2)


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# retrieve the module labels 
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
df.module.labels <- as.data.frame(moduleLabels)



MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "ncRNA_wo_mirna_manatee-01-networkConstruction-auto.RData")


#calculate the module relatinship using kME

hub.genes <- signedKME(wgcna_input_final, MEs)

#rownames(vsd_counts_symbols) <- vsd_counts_symbols$ensembl_gene_id

# get the dataframe expression with labels
final_gene_dataframe = merge(x = df.module.labels, y = vsd_counts, by=0,
                             all.x = TRUE)


write.csv(final_gene_dataframe, file = paste0(outputPrefix, "gene_expression_with_module_wo_mirna.csv"))
write.csv(hub.genes, file = paste0(outputPrefix, "kME_wo_mirna_modules.csv"))

hub_genes <- chooseTopHubInEachModule(wgcna_input_final, moduleColors)
#########################################################################
#load the trait data as binary data representing the time to evaluate positiv relationships
# here trait data are the transcription factors
data_traits <- read.csv("Transcripton_Factor_expression.csv", header = TRUE, row.names = 1)
nGenes = ncol(wgcna_input_final);
nSamples = nrow(wgcna_input_final);

write.csv(MEs, "Eigenvector_modalities_wo_mirna.csv")

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
dissTOM = 1-TOMsimilarityFromExpr(wgcna_input_final, power = 14);
adj <- dissTOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
adjacency_graph = igraph::graph.adjacency(dissTOM)
network <- graph.adjacency(adj)
network <- simplify(network)  # removes self-loops
V(network)$color <- moduleLabels
network <- delete.vertices(network, degree(network)==0)
plot(network, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)

#########################################################################
#trajectories of genes
#########################################################################

# give sample names 
columns <- c("D_01",
             "D_01",
             "D_01",
             "D_05",
             "D_05",
             "D_05",
             "D_09",
             "D_09",
             "D_09",
             "D_16",
             "D_16",
             "D_16",
             "D_26",
             "D_26",
             "D_26",
             "D_36",
             "D_36",
             "D_36",
             "D_01",
             "D_01",
             "D_01",
             "D_05",
             "D_05",
             "D_05",
             "D_09",
             "D_09",
             "D_09",
             "D_16",
             "D_16",
             "D_16",
             "D_26",
             "D_26",
             "D_26",
             "D_36",
             "D_36",
             "D_36",
             "D_01",
             "D_01",
             "D_01",
             "D_05",
             "D_05",
             "D_05",
             "D_09",
             "D_09",
             "D_09",
             "D_16",
             "D_16",
             "D_16",
             "D_26",
             "D_26",
             "D_26",
             "D_36",
             "D_36",
             "D_36"
)

# retrieve the Eigenvector modules
eigenvector_modules = as.data.frame(t(net$MEs))
colnames(eigenvector_modules) <- columns
eigenvector_modules = eigenvector_modules[,order(names(eigenvector_modules))]

# Draw the Heatmap for the Eigenvector Modules to obtain module trajectories
pheatmap::pheatmap(eigenvector_modules, 
                   cluster_cols = F,
                   color=colorRampPalette(c("purple", "black", "orange"))(75),
                   scale="row",
                   border_color="grey",
                   show_colnames = F,
                   fontsize_number = 7,
                   gaps_col = c(9,18,27,36,45),
                   cellwidth =10,
                   cellheight = 10)









