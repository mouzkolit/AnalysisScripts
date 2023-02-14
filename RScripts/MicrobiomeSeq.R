library(dada2)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(tidyr)
library("ape")

# set the path
path <- "./21172_RawData"
list.files(path) # list all the file names

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# plot the quality of the foward and reverse strand
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])


#adjust the naming for the filtering procedure
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filtering was performed with suggested parameters form a dada2 tutorial
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,240), trimLeft = c(17,21),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread = T) # On Windows set multithread=FALSE

# dereplicated
derepF1 <- derepFastq(filtFs, verbose=TRUE)
derepR1 <- derepFastq(filtRs, verbose=TRUE)

# start learning the error rates
errF <- learnErrors(derepF1, multithread=TRUE)
errR <- learnErrors(derepR1, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# get the filtered and error corrected
dadaFs <- dada(derepF1, err=errF, multithread=TRUE)
dadaRs <- dada(derepR1, err=errR, multithread=TRUE)

# check it 
dadaFs[[1]]


# merge the pairs
mergers <- mergePairs(dadaFs, derepF1, dadaRs, derepR1, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# make sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


#check reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# retrieve the taxonomy of the species
tax <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa", multithread=TRUE)
tax <- addSpecies(tax, "silva_species_assignment_v138.1.fa")

saveRDS(seqtab.nochim, "seqtab_gp130_nochim.rds") 
saveRDS(tax, "tax_final.rds")

###########################################################################
#seqtab.nochim <- readRDS("seqtab_gp130_nochim.rds")
#tax   <- readRDS("tax_final.rds")
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(ID = samples.out)
row.names(samdf) = samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

count_table_tax  = t(rbind(ps@otu_table, t(ps@tax_table)))

count_table_one_word = data.frame(count_table_tax[,1:(nrow(ps@otu_table))])
View(head(count_table_one_word))
count_table_one_word$taxonomy = paste0(count_table_tax[,"Kingdom"], ";" , 
                                       count_table_tax[,"Phylum"], ";" , 
                                       count_table_tax[,"Class"], ";" , 
                                       count_table_tax[,"Order"], ";" , 
                                       count_table_tax[,"Family"], ";" , 
                                       count_table_tax[,"Genus"], ";")


# get the aggregated table
agg <- aggregate_taxa(ps, level = "Family", verbose = TRUE)
agg_kingdom <- aggregate_taxa(ps, level = "Kingdom", verbose = TRUE)


## save the aggregated table
colap_kingdom <- unite(as.data.frame(agg_kingdom@tax_table@.Data), newCol, -unique) 
kingdom_table <- agg_kingdom@otu_table@.Data

colap <- unite(as.data.frame(agg@tax_table@.Data), newCol, -unique) 

genus_table <- agg@otu_table@.Data
row.names(genus_table) <- colap$newCol
write.csv(x = genus_table, quote = F, file = "genus_table_from_dada2_pain_family.csv")

######################################################################
seqtab.nochim <- readRDS("seqtab_gp130_nochim.rds")
tax   <- readRDS("tax_final.rds")
metadata <- read.csv("sample-metadata.csv", header = TRUE)
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(ID = samples.out, condition = metadata$condition, gender = metadata$gender)
row.names(samdf) = samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

count_table_tax  = t(rbind(ps@otu_table, t(ps@tax_table)))

count_table_one_word = data.frame(count_table_tax[,1:(nrow(ps@otu_table))])
count_table_one_word$taxonomy = paste0(count_table_tax[,"Kingdom"], ";" , 
                                       count_table_tax[,"Phylum"], ";" , 
                                       count_table_tax[,"Class"], ";" , 
                                       count_table_tax[,"Order"], ";" , 
                                       count_table_tax[,"Family"], ";" , 
                                       count_table_tax[,"Genus"], ";")


#Use phyloseq for first overview
ps <- prune_samples(sample_names(ps) != "Mock", ps)

#get alpha diversity plot using shannon
plot_richness(ps, x="condition", measures=c("Shannon", "Simpson","Chao1"), color="gender")

#transform data into proportions for distance measurement
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

# get the bray curtis plot for the analysis
plot_ordination(ps.prop, ord.nmds.bray, color="gender", title="Bray NMDS", shape = "condition")

# get the top20 families 
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="condition",fill="Family") + facet_wrap(~gender, scales="free_x")

####
# retrieve the phylogenetic trees
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
physeq1 = merge_phyloseq(ps,random_tree)
plot_bar(physeq1, fill="Family")
# plot the taxa phylotree
myTaxa = names(sort(taxa_sums(physeq1), decreasing = TRUE)[1:10])
ex1 = prune_taxa(myTaxa, physeq1)
plot(phy_tree(ex1), show.node.label = TRUE)
plot_heatmap(physeq1, taxa.label="Phylum")

# make the differential expression analysis using DESeq2
library("DESeq2")

diagdds = phyloseq_to_deseq2(physeq1, ~ condition*gender)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

# have a look at significances
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
head(sigtab)
