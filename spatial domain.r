library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(RColorBrewer)
library(pheatmap)

sp_harmony <- readRDS('/data/work/MP/scRNA/spatial_domain/shiyan-you_duizhao-you.Seurat.rds')

my_colors <- c(
  "0" = "#8DD3C7",  # domain 0
  "1" = "#FFFFB3",  # domain 1
  "2" = "#BEBADA",  # domain 2
  "3" = "#FB8072",  # domain 3
  "4" = "#80B1D3"   # domain 4
)
SpatialDimPlot(object = sp_harmony, group.by = "seurat_clusters",cols=my_colors,pt.size.factor = 1.5)
ggsave('/data/work/MP/scRNA/spatial_domain/shiyan-you_duizhao-you_domain.pdf', width = 12, height = 12)

#spatial DEG
DefaultAssay(sp_harmony) <- "Spatial"
Idents(sp_harmony) <- sp_harmony$`seurat_clusters`
sp_harmony.markers <- FindAllMarkers(sp_harmony, only.pos = TRUE)
markers <- sp_harmony.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC >0.5)  %>%
    dplyr::filter(p_val_adj <0.05)
write.csv(markers, file = "/data/work/MP/scRNA/spatial_domain/domain_markers.csv")

# split
sp_list <- SplitObject(sp_harmony, split.by = "stim")
names(sp_list)  
shiyan_obj <- sp_list[["shiyan-you"]]
duizhao_obj <- sp_list[["duizhao-you"]]

### control group
rctd_weights <- read.csv("/data/work/MP/scRNA/RCTD/control/control_norm_weights.csv", row.names = 1)
rownames(rctd_weights) <- gsub("\\.", "-", rownames(rctd_weights))
head(rctd_weights)

composition <- rctd_weights
rownames(composition) <- paste0(rownames(composition), "_2")

common_cells <- intersect(rownames(composition), colnames(duizhao_obj))

control_ST_matched <- AddMetaData(
  duizhao_obj,
  metadata = composition[common_cells, ]  # ensure consistent ordering
)

# Step 1: Extract the cell type columns of interest
celltypes <- c(
  "T.cells", "B.cells", "Monocytes", "Endothelial.cells",
  "Erythrocytes", "Fibroblasts", "Macrophages",
  "Pericytes", "Platelets", "Epithelial.cells"
)  # replace with your actual column names

meta <- control_ST_matched@meta.data

# Step 2: Convert seurat_clusters to factor (to ensure controllable ordering)
meta$seurat_clusters <- as.factor(meta$seurat_clusters)

# Step 3: Compute the mean value of each cell type column within each cluster
mean_scores <- meta %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(celltypes), mean, na.rm = TRUE))

# Step 4: Convert to matrix and set row names as cell types
score_matrix <- as.matrix(t(mean_scores[, -1]))  # remove cluster column before transposing
colnames(score_matrix) <- paste0("Cluster_", mean_scores$seurat_clusters)
rownames(score_matrix) <- celltypes

# Step 5: Plot heatmap
pheatmap(
  score_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "red"))(100),
  main = "Avg RCTD Cell Type Score per Seurat Cluster"
)

write.csv(
  score_matrix,
  file = "/data/work/MP/scRNA/spatial_domain/control_score_martix.csv"
)

### model group
rctd_weights <- read.csv("/data/work/MP/scRNA/RCTD/model/model_norm_weights.csv", row.names = 1)

# Replace '.' with '-' in row names
rownames(rctd_weights) <- gsub("\\.", "-", rownames(rctd_weights))

# Inspect the result
head(rctd_weights)

composition <- rctd_weights
rownames(composition) <- paste0(rownames(composition), "_1")

common_cells <- intersect(rownames(composition), colnames(shiyan_obj))

shiyan_ST_matched <- AddMetaData(
  shiyan_obj,
  metadata = composition[common_cells, ]  # ensure consistent ordering
)

# Step 1: Extract the cell type columns of interest
celltypes <- c(
  "T.cells", "B.cells", "Monocytes", "Endothelial.cells",
  "Erythrocytes", "Fibroblasts", "Macrophages",
  "Pericytes", "Platelets", "Epithelial.cells", "Neutrophils"
)  # replace with your actual column names

meta <- shiyan_ST_matched@meta.data

# Step 2: Convert seurat_clusters to factor (to ensure controllable ordering)
meta$seurat_clusters <- as.factor(meta$seurat_clusters)

# Step 3: Compute the mean value of each cell type column within each cluster
mean_scores <- meta %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(celltypes), mean, na.rm = TRUE))

# Step 4: Convert to matrix and set row names as cell types
score_matrix <- as.matrix(t(mean_scores[, -1]))  # remove cluster column before transposing
colnames(score_matrix) <- paste0("Cluster_", mean_scores$seurat_clusters)
rownames(score_matrix) <- celltypes

# Step 5: Plot heatmap
pheatmap(
  score_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "red"))(100),
  main = "Avg RCTD Cell Type Score per Seurat Cluster"
)

write.csv(
  score_matrix,
  file = "/data/work/MP/scRNA/spatial_domain/model_score_martix.csv"
)
