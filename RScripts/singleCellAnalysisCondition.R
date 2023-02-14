library("Seurat")
library(dplyr)
library(purrr)
library(metap)
library(tibble)
library(DESeq2)
library(ggplot2)
library(patchwork)
library(limma)
library(EnhancedVolcano)
library(gprofiler2)
library(ggplot2)
library(celldex)
library(SingleR)
library(clusterProfiler)
library(org.Mm.eg.db)
library("monocle3")
library(AnnotationHub)
library("glmGamPoi")
hubCache(AnnotationHub())
#load the dataset object
ctrl.data <- Read10X(data.dir = "./20098-0001/filtered_feature_bc_matrix")
fabry.data <- Read10X(data.dir = "./20098-0002/filtered_feature_bc_matrix")


fabry_ctrl = CreateSeuratObject(counts = ctrl.data, project = "ctrl", min.cells = 3, min.features = 200)
fabry_cog <- CreateSeuratObject(counts = fabry.data, project = "ctrl", min.cells = 3, min.features = 200)


# check the quality of the dataset and remove high mitochondrial counts
fabry_ctrl[["percent.mt"]] <- PercentageFeatureSet(fabry_ctrl, pattern = "^mt-")
fabry_cog[["percent.mt"]]  <- PercentageFeatureSet(fabry_cog, pattern = "^mt-")

#make a control
FeatureScatter(fabry_ctrl, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
FeatureScatter(fabry_cog, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)


v1 = VlnPlot(fabry_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
v2 = VlnPlot(fabry_cog, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
v1
v2

ctrl<- subset(fabry_ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 6)
fabry<- subset(fabry_cog, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 6)

#experiment list

fabry_ctrl$condition <- "fabry"
fabry_cog$condition <- "cognitive"


experiment.list <- c(fabry_ctrl, fabry_cog)

# run the SCTransform

for (i in 1:length(experiment.list)) {
  experiment.list[[i]] <- SCTransform(experiment.list[[i]], method = "glmGamPoi", verbose = T, vars.to.regress = c("percent.mt"), return.only.var.genes = T,variable.features.n = 3000)
}

# select the feeatures for downstream integration 
experiment.features <- SelectIntegrationFeatures(object.list = experiment.list, nfeatures = 3500)
experiment <- PrepSCTIntegration(object.list = experiment.list, anchor.features = experiment.features, 
                                 verbose = TRUE)

experiment.anchors <- FindIntegrationAnchors(object.list = experiment, normalization.method = "SCT", 
                                             anchor.features = experiment.features, verbose = TRUE, dims = 1:30)
experiment.integrated <- IntegrateData(anchorset = experiment.anchors, normalization.method = "SCT", 
                                       verbose = TRUE, dims = 1:35)


# set the assay to an integrated analysis 

# Run the standard workflow for visualization and clustering
experiment <- RunPCA(experiment.integrated, npcs =25,verbose = FALSE)
# t-SNE and Clustering
experiment <- RunUMAP(experiment, reduction = "pca", dims = 1:25)
experiment <- FindNeighbors(experiment, reduction = "pca", dims = 1:25)
experiment <- FindClusters(experiment, resolution = 0.5,random.seed = 42)


#Draw the plots to identify if integration worked properly
p1 <- DimPlot(experiment, reduction = "umap", group.by = "condition")
p2 <- DimPlot(experiment, reduction = "umap", label = TRUE)
p1 +p2

DimPlot(experiment, reduction = "umap", split.by = "condition")

###########################################################################

DefaultAssay(experiment) <- "RNA"

experiment <- NormalizeData(experiment, verbose = TRUE)
experiment <- ScaleData(experiment, verbose = TRUE)


#get the cluster annotation based on singleR
counts <- GetAssayData(experiment, assay = "RNA", slot = "data")
ref <- celldex::MouseRNAseqData()
singler <- SingleR(counts,
                   ref = ref,
                   labels = ref$label.fine)

plotScoreHeatmap(singler)
plotDeltaDistribution(singler)

##### get the diagnostics
all.markers <- metadata(singler)$de.genes
olig.markers <- unique(unlist(all.markers$Oligodendrocytes))
experiment$labels <- singler$labels

# plot the cell types from the experiment
plot_experiment_singler <- experiment
plot_experiment_singler <- SetIdent(plot_experiment_singler, value ="labels")
DimPlot(plot_experiment_singler) # plot the umap after cell type identificaiton


saveRDS(experiment, "Seurat_analysis_integration_cognition.rds")

###############################################################################

experiment <- readRDS("Seurat_analysis_integration.rds") # read in the file if you prefer

experiment <- RenameIdents(experiment, 
                           `0` = "T-cells",
                           `1` = "T-cells",
                           `2` = "T-cells", 
                           `3` = "Immune_cells_2",
                           `4` = "Oligo_2",
                           `5` = "Ependymal",
                           `6` = "qNSC_1",
                           `7` = "Ependymal",
                           `8` = "Microglia_2",
                           `9` = "Oligo_1", 
                           `10` = "qNSC_2",
                           `11` = "Ependymal",
                           `12` = "aNSC",
                           `13` = "T_cells",
                           `14` = "Unknown",
                           `15` = "Endothelial_cells",
                           `16` = "Microglia_1",
                           `17` = "Monocytes",
                           `18` = "Neurons",
                           `19` = "T_cells",
                           `20` = "aNSC",
                           `21` = "Endothelial_cells",
                           `22` = "Immune_cells_1",
                           `23` = "Monocytes",
                           `24` = "T-cells"
)


# get conserved markers throughout both conditions
get_conserved <- function(cluster){
  try(FindConservedMarkers(experiment,
                           ident.1 = cluster,
                           grouping.var = "condition",
                           only.pos = TRUE) %>%
        rownames_to_column(var = "gene") %>%
        cbind(cluster_id = cluster, .))
}
conserved_markers <- map_dfr(unique(Idents(experiment)), get_conserved)

#reduce top conserved markers
top_conserved_marker <- conserved_markers %>% group_by(cluster_id) %>% top_n(n = -10, wt = max_pval) %>% top_n(n = 10, wt = fabry_avg_log2FC)

# heatmap of top conserved marker
DoHeatmap(experiment, features = top_conserved_marker$gene)

#stacked violinplot of conserved marker for expression
StackedVlnPlot(obj = experiment, features = c("Mbp","Sept4","Cldn5","Ly6c1","Cx3cr1",
                                              "Rgs10","C1qb","Camk2a","Camk2b","Ntsr2","Htra1",
                                              "Ctsw","H2-Q7"
))


#Draw a distributino table for all cell clusters
pt <- table(Idents(experiment), experiment$condition)
pt <- as.data.frame(pt)
pt$Cluster <- as.character(pt$Var1)
pt <- pt[order(pt$Freq),]
colors <- colorRampPalette(c("red", "green","blue","yellow","skyblue","orange","grey","violet","ivory","pink"))(32)

ggplot(pt, aes(x = Var2, y = Freq, fill = Cluster, cols = colors)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill",size=0.1, alpha = 0.7,  width = 0.5, color = "black") +
  xlab("Sample") +
  ylab("Proportion")


#get differentially expressed genes for everything
all_degs<- experiment
Idents(all_degs) <- "condition"
avg.all_degs<- as.data.frame(log1p(AverageExpression(all_degs, verbose = FALSE)$RNA))
avg.all_degs$gene <- rownames(avg.all_degs)

all_degs.tested <- FindMarkers(all_degs, ident.1 = "fabry", ident.2 = "cognitive")

all_degs.up <- all_degs.tested[all_degs.tested$avg_log2FC > 0,] # upregulated genes
all_degs.down <- all_degs.tested[all_degs.tested$avg_log2FC < 0,] # downregulated genes

#Draw Volcano for up and downregulated genes 
EnhancedVolcano(all_degs.tested,
                lab = rownames(all_degs.tested),
                x = "avg_log2FC",
                y = "p_val_adj",
                FCcutoff = 1,
                drawConnectors = T)


# get enrichment for all degs
enrichments <- enrichGO(gene = rownames(all_degs.down),
                        OrgDb =  org.Mm.eg.db,
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.01,
                        keyType = "SYMBOL")

dotplot(enrichments, showCategory = 15)

FeaturePlot(experiment, reduction = "umap", c("Ttr"), min.cutoff = "q10", cols = c("skyblue","red"), split.by = "condition") # olig
