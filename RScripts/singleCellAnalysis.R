library(data.table)
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(gprofiler2)


# preprocess the count matrix 
trial = readRDS("GSE154659_C57_Raw_counts.RDS")
columns = colnames(trial)
# get the indeces
index <- 0
vector <- c()

#only retrieve the ones associated with crush injury
for(i in columns){
  index <- index + 1
  if (grepl("^male_C57_Crush_168", i) == TRUE) {
    print("hello")
    print(index)
    vector <- c(vector, i)
  }
  
  if (grepl("^male_C57_Naive_0", i) == TRUE) {
    print("hello")
    print(index)
    vector <- c(vector, i)
  }
}
matrix <- trial[,vector]

# make a metadata table
metadata <- data.frame(do.call(rbind, strsplit(colnames(matrix), "_", fixed=TRUE)))
colnames(metadata) <- c("sex","mouse","condition","timepoint","batch","celltype","barcode","gender")
metadata$unique_barcodes <- make.unique(metadata$barcode)
rownames(metadata) <- metadata$unique_barcodes
colnames(matrix) <- metadata$unique_barcodes

# load seurat object 
crush <- CreateSeuratObject(
  matrix,
  project = "CreateSeuratObject",
  assay = "RNA",
  meta.data = metadata
)
Idents(crush) <- crush@meta.data$celltype
crush_neuron <- subset(crush, idents = c("NF1", "NF2", "NF3", "NP","p", "PEP1", "PEP2", "SST", "cLTMR1"))

experiment <- PercentageFeatureSet(crush_neuron, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

experiment <- SCTransform(experiment, vars.to.regress = "percent.mt",  verbose = T, vst.flavor = "v2")
experiment <- RunPCA(experiment, verbose = FALSE)
experiment <- RunUMAP(experiment, dims = 1:30, verbose = FALSE)


# Find Markers fucntion
experiment <- FindNeighbors(experiment, dims = 1:30, verbose = FALSE)
experiment <- FindClusters(experiment, verbose = FALSE, resolution = 0.3)


DefaultAssay(experiment) <- "RNA"
experiment = NormalizeData(experiment)
experiment = ScaleData(experiment)
saveRDS(experiment, "crush_planned.RDS")

experiment <- readRDS("crush_planned.RDS")
experiment@meta.data$injure_condition <- ifelse(experiment@meta.data$seurat_clusters %in% c(1,5),"injured", "uninjured")
# Check the Results
p2 = DimPlot(experiment, group.by = "celltype")
p3 = DimPlot(experiment, group.by = "condition")
p5 = DimPlot(experiment, group.by = "injure_condition")
p2 + p3 + p5
##################################################
# Downstream Analysis
# check the cell type frequency per condition to retrieve injured neuron clusters
condition_counts = experiment@meta.data
counts = as.data.frame(table(condition_counts$seurat_cluster,condition_counts$condition))
final_counts = counts %>% group_by(Var1, Var2) %>% summarise(percentage = sum(Freq)) %>%
  mutate(freq = percentage/ sum(percentage))

ggplot(final_counts, aes(fill = Var2, y=freq, x=Var1, label = round(freq, 2))) +  
  geom_bar(stat="identity")  +  
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  theme_light()+
  theme(text = element_text(size = 17))     
  

################################################
# add an injured and uninjured tag

# load 6540 and sni markers
table_6540 = read.csv("mir_6540_-results-with-normalized_significant.csv", header = T)
sni_markers = read.csv("_mRNA_with_cluster.csv")
sni_markers = sni_markers[sni_markers$cluster %in% c(3,5),]
table_6540 = table_6540[table_6540$log2FoldChange > 0,]


# markers_based on injured uninjured state
Idents(experiment) <- experiment@meta.data$injure_condition
markers_injured = FindAllMarkers(experiment, assay = "RNA", slot = "data")
markers_injured = markers_injured[markers_injured$cluster =="injured",]
DotPlot(experiment, features = c("Atf3","Sprr1a","Gal","Gap43","Cacna2d1","Tac1","Grik1"), group.by = "injure_condition", cols = c("yellow", "purple"))
DotPlot(experiment, features = c("Atf3","Sprr1a","Gal","Gap43","Cacna2d1","Tac1","Grik1"), group.by = "seurat_clusters", cols = c("yellow", "purple"))
FeaturePlot(experiment, features = c("Atf3","Sprr1a","Gal","Gap43","Cacna2d1","Tac1","Grik1"))
markers_6540 = markers_injured[markers_injured$gene %in% table_6540$external_gene_name,]
markers_6540 = markers_injured[markers_injured$gene %in% table_6540$external_gene_name,]

#based on clustering
# Test used is Mast
Idents(experiment) <- experiment@meta.data$seurat_clusters
markers = FindAllMarkers(experiment, assay = "RNA", slot = "data", test.use = "MAST")
markers = markers[markers$cluster %in% c(1,5),]
markers = markers[markers$p_val_adj < 0.1,]
markers_cluster_6540 = markers[markers$gene %in% table_6540$external_gene_name,]

#get intersections
markers_crush = unique(markers$gene)
sni_6540 = unique(intersect(sni_markers$ensembl_gene_id, table_6540$ensembl_gene_id))
sni_crush_overlap = unique(intersect(markers$gene, sni_markers$external_gene_name))
crush_6540 = unique(intersect(markers$gene, table_6540$external_gene_name))
all_overlap = Reduce(intersect, list(markers$gene, table_6540$external_gene_name, sni_markers$external_gene_name))

# draw a venn_diagram over the distribution
library("VennDiagram")    
grid.newpage()
draw.triple.venn(3460, 
                 2985, 
                 2103, 
                 1079, 
                 365,
                 645,
                 160,
                 c("SNI",
                   "Crush Injury",
                   "6540-/-"),
                 fill = c("lightyellow", "lightgreen", "skyblue"))


# retrieve the gprofiling results for the overlap
save_table <- function(genes, table_name){
  enrichments <- gost(query = genes,
                                 organism = "mmusculus") 
  subsetted <- subset(enrichments$result, select = -c(parents))
  write.csv(subsetted, table_name)
}

save_table(sni_6540, "sni_6540_overlap.csv")
save_table(crush_6540, "crush_6540_overlap.csv")
save_table(all_overlap, "crush_sni_6540_overlap.csv")

markers_6540_ordered = markers[markers$gene %in% crush_6540,]
markers_6540_ordered = markers_6540_ordered[order(markers_6540_ordered$avg_log2FC),]
negative_fold_change = head(markers_6540_ordered$gene, 100)
positive_fold_change = tail(markers_6540_ordered$gene, 100)
all_negative_positive = c(negative_fold_change, positive_fold_change)


DoHeatmap(experiment, features = all_negative_positive, draw.lines = T) + NoLegend()+ 
  theme(text = element_text(size = 0))


# here I should 
top_down_features = c("Meg3",
                    "Syt7",
                    "Scn10a",
                    "Kcna2",
                    "Prkca",
                    "Phf24",
                    "Als2",
                    "Scn11a",
                    "Tshz2",
                    "Ahnak",
                    "Pde11a",
                    "Deptor",
                    "Tshz2",
                    "Syt1",
                    "Cpne3",
                    "Itch",
                    "Prkca",
                    "Osbpl3")

top_up_features <- c("Cacna2d1",
                       "Sez6l",
                       "Plxna4",
                       "Myo10",
                       "Mmp16",
                       "Plcxd2",
                       "Uhmk1",
                       "Nrip1",
                       "Kcnq3",
                       "Irs2",
                       "Slc24a2",
                       "Phactr2",
                       "Irs2",
                       "Ahnak2",
                       "Map1b",
                       "Arrb1",
                       "Ndst1")
p1 = VlnPlot(experiment, features = top_up_features , stack=T, sort = TRUE, flip = TRUE, group.by = "injure_condition")
p2 = VlnPlot(experiment, features = top_down_features, stack=T, sort = TRUE, flip = TRUE, group.by = "injure_condition")
p1+ p2

FeaturePlot(experiment, features = top_up_features)
FeaturePlot(experiment, features = top_down_features)

fwrite(list(all_overlap), file = "overlapping_genes.csv")
fwrite(list(crush_6540), file = "gene_crush_6540.csv")

library(gprofiler2)

gostres <- gost(query = sni_6540, organism = "mmusculus")
enrichment <- as.data.frame(gostres$result)
enrichment <- subset(enrichment, select = -c(parents))
write.csv(enrichment, "SNI_6540_overlap_enrichments.csv")
