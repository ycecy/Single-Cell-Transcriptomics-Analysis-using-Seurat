# Install the necessary packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("enrichplot")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", force = TRUE)

BiocManager::install("DESeq2")

BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("openxlsx", quietly = TRUE))
  install.packages("openxlsx")

if (!requireNamespace("biomaRt", quietly = TRUE))
  
  install.packages("Seurat")
install.packages("SeuratObject")

install.packages("VennDiagram")
install.packages("ggVennDiagram")

BiocManager::install("goseq")

install.packages("cowplot")

BiocManager::install("speckle")


# Set the directory
setwd("~/ycecy/file")

library(dplyr)
library(Seurat)
library(patchwork)

# Upload the COPD - Influenza A Infection dataset
data_copd <- readRDS("/Users/ycecy/file/data_copd.rds") 

# Annotation of data_copd before mapping
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

pl <- DimPlot(data_copd, reduction = "umap") +
  theme_minimal() +
  ggtitle("Annotation Before Mapping")
ggsave(filename = "Annotation Before Mapping.png",
       plot=pl, width = 10, height = 8, dpi = 300)

# Annotate the cells from data_copd with reference dataset below (Re-do the mapping)
data_ref1 <- readRDS("/Users/ycecy/file/Madissoon et al.rds")
data_ref2 <- readRDS("/Users/ycecy/file/Madissoon et al_fibroblast.rds")
data_ref <- merge(data_ref1, y = data_ref2, project = "data_ref_all")
bronch_data_ref = subset(data_ref, subset = tissue == "bronchus")

SaveSeuratRds(bronch_data_ref,"/Users/ycecy/file/data_ref_merged.rds")

# Upload the merged bronchus reference dataset
bronchdata_ref <- readRDS("/Users/ycecy/file/data_ref_merged.rds")

head(rownames(bronchdata_ref))
head(rownames(data_copd))
nrow(data_copd)

# Extract the gene symbols for getBM function
gene_symbols <- rownames(GetAssayData(data_copd, assay = "RNA", layer = "data"))

library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)

# Converting ENSG IDs to Gene Symbols for getBM function
sytoens <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

biomaRt_results <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"), 
  filters = "hgnc_symbol",                         
  values = gene_symbols,                        
  mart = sytoens)

head(biomaRt_results)

# Extract the count matrix from data_copd Seurat object to create the new Seurat data
# to obtain non-normalized data and prevent re-normalizing
count.data = data_copd@assays[["RNA"]]@counts
count.data <- as.matrix(count.data)
uniKeys_data_copd <- rownames(count.data)

# Map the HGCN symbols to Ensembl ID's 
gene_mapping <- biomaRt_results[match(uniKeys_data_copd, biomaRt_results$hgnc_symbol), "ensembl_gene_id"]

# Remove NA values
gene_mapping[is.na(gene_mapping)] <- uniKeys_data_copd[is.na(gene_mapping)]

# Remove the duplicates
rownames(count.data) <- make.unique(gene_mapping)

# Create the new Seurat object
ref <- CreateSeuratObject(
  count.data,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_")

# Normalize the data
ref <- NormalizeData(ref)

# Feature selection identification
ref <- FindVariableFeatures(ref)

# Scale the data
ref <- ScaleData(ref)

# Find the anchors between the 'bronchdata_ref' reference data and 'ref' query object
bronch.anchors <- FindTransferAnchors(
  reference = bronchdata_ref,
  query = ref,
  reference.assay = "RNA",
  query.assay = "RNA",
  project.query = TRUE)

# Transfer cell type labels
predictions <- TransferData(
  anchorset = bronch.anchors,
  refdata = bronchdata_ref$Celltypes,
  dims = 2:30)

# Add the metadata of interest for new Seurat object
ref <- AddMetaData(ref, metadata = predictions)
ref <- AddMetaData(ref, metadata = data_copd$sample,col.name = "sample")
ref <- AddMetaData(ref, metadata = data_copd$viralLoad,col.name = "viralLoad")
ref <- AddMetaData(ref, metadata = data_copd$donor,col.name = "donor")


# Principal component analysis for linear dimensional reduction on the scaled data
ref <- RunPCA(ref, features = VariableFeatures(object = ref))

# Visualization of the dimensions after PCA
DimPlot(ref, reduction = "pca") + NoLegend()
VizDimLoadings(ref, dims = 1:2, reduction = "pca") 
DimHeatmap(ref, dims = 1, cells = 300, balanced = TRUE)

# Check the dimensionality of the dataset
ElbowPlot(ref)

# Cluster the cells
ref <- FindNeighbors(ref, dims = 1:10)
ref <- FindClusters(ref, resolution = 0.5)

# Run UMAP as non-linear dimensional reduction
ref <- RunUMAP(ref, dims = 1:10)

SaveSeuratRds(ref,"/Users/ycecy/file/copd_annotated_new.rds")

# Load the newly annotated reference dataset 
ref <- readRDS("/Users/ycecy/file/copd_annotated_new.rds")


# Newly annotated cells before merging into bigger groups
pl <- DimPlot(ref, reduction = "umap", group.by = "predicted.id")+
  theme_minimal() +
  ggtitle("Annotated Cells")
ggsave(filename = "Annotated Cells.png",
       plot=pl, width = 10, height = 8, dpi = 300)

# Merge subgroups to bigger ones for Differential Expression Analysis 
ref$celltype.new <- ref$predicted.id

ref$celltype.new[ref$predicted.id %in% c("Basal", "Dividing_Basal")] <- "Basal"
ref$celltype.new[ref$predicted.id %in% c("Secretory_Club", "Secretory_Goblet")] <- "Secretory"
ref$celltype.new[ref$predicted.id %in% c("Fibro_alveolar", "Fibro_adventitial", "Fibro_perichondrial",
                                         "Mesothelia", "Fibro_myofibroblast", "Fibro_peribronchial")] <- "Fibroblasts"
#ref$celltype.new[ref$predicted.id %in% c("SMG_Duct", "SMG_Mucous", "SMG_Serous")] <- "SMG"

table(ref$celltype.new)

# Newly annotated cells after merging
pl <- DimPlot(ref, reduction = "umap", group.by = "predicted.id")+
  theme_minimal() +
  ggtitle("Annotated Cells")
ggsave(filename = "Annotated Cells.png",
       plot=pl, width = 10, height = 8, dpi = 300)

# Influenza A infection across samples
pl <- DimPlot(ref, group.by = "sample",reduction = "umap",
              cols = c("lightgray", "darkred"), pt.size = 0.5) +
  theme_minimal() +
  ggtitle("Influenza A Infection Across Samples")
ggsave(filename = "Influenza A Infection Sample.png",
       plot=pl, width = 10, height = 8, dpi = 300)


# Create metadata for COPD condition of the donors
ref$copd_condition <- ifelse(ref$donor %in% c("donor4", "donor5", "donor6"),
                             "copd+", 
                             "copd-")

# Create metadata for infection status according to the set threshold
# Log transformed viral load
ref$virallog10 <- log10(ref$viralLoad)

# Omit NA values
ref$virallog10[is.infinite(ref$virallog10)] <- NA

# Influenza A infection across samples over viral load
pl <- FeaturePlot(ref, features = "virallog10",reduction = "umap",
                  cols = c("lightgray", "darkred"), pt.size = 0.5) +
  theme_minimal() +
  ggtitle("Log Transformed Viral Load Across Samples")
ggsave(filename = "Log Transformed Viral Load Across Samples.png",
       plot=pl, width = 10, height = 8, dpi = 300)

# Set the threshold
ref$inf_status <- ifelse(ref$viralLoad == 0, "uninfected",
                         ifelse(ref$viralLoad > 0 & ref$viralLoad <= 0.0001, "low",
                                ifelse(ref$viralLoad > 0.0001 & ref$viralLoad <= 0.001, "medium", "high")))
table(ref$inf_status)

# Influenza A infection degree across samples
pl <- DimPlot(ref, group.by = "inf_status",reduction = "umap",
              cols = c("black", "red", "lightpink", "lightblue"), pt.size = 0.1) +
  theme_minimal() +
  ggtitle("Influenza A Infection Degree Across Samples")
ggsave(filename = "Influenza A Infection Degree Across Samples.png",
       plot=pl, width = 10, height = 8, dpi = 300)


# Influenza A infection degree shaped across cell types
# Create the DimPlot
umap_plot <- DimPlot(ref, reduction = "umap", group.by = "celltype.new", pt.size = 1)

# Extract the UMAP coordinates
umap_data <- FetchData(ref, vars = c("umap_1", "umap_2", "celltype.new", "inf_status"))

pl <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = celltype.new, shape = inf_status)) +
  geom_point(size = 0.8, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Infection degree-shaped across cell types") +
  scale_shape_manual(values = c(16, 3, 22, 6))
ggsave(filename = "Infection degree-shaped across cell types.png",
       plot=pl, width = 10, height = 8, dpi = 300)


# Create a table of counts of infection status of the cell types
inf_by_celltype <- table(umap_data$celltype.new, umap_data$inf_status)
print(inf_by_celltype)

# Save the result table as excel file
library(openxlsx)
write.xlsx(inf_by_celltype, "/Users/ycecy/file/infbycelltype.xlsx")


# Visualization of donors
pl <- DimPlot(ref, group.by = "donor",reduction = "umap",
        cols = c("lightgray", "darkred", "lightpink", "lightblue", "darkgreen"), pt.size = 0.5) +
  theme_minimal() +
  ggtitle("Donor Information")
ggsave(filename = "Donor Information.png",
       plot=pl, width = 10, height = 8, dpi = 300)

# Visualization of donors and infection status
umap_plot <- DimPlot(ref, reduction = "umap", group.by = "donor", pt.size = 1)
umap_data <- FetchData(ref, vars = c("umap_1", "umap_2", "donor", "virallog10"))

# Create ggplot
pl <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = virallog10, shape = donor)) +
  geom_point(size = 0.8, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Donor and viral load across cell types") +
  scale_color_viridis_c() +
  scale_shape_manual(values = 1:10)
ggsave(filename = "Donor and viral load across across cell types.png",
       plot=pl, width = 10, height = 8, dpi = 300)


# COPD across cell types
umap_plot <- DimPlot(ref, reduction = "umap", group.by = "celltype.new", pt.size = 1)
umap_data <- FetchData(ref, vars = c("umap_1", "umap_2", "celltype.new", "copd_condition"))

pl <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = celltype.new, shape = copd_condition)) +
  geom_point(size = 0.8, alpha = 0.8) +
  theme_minimal() +
  labs(title = "COPD-shaped across cell types") +
  scale_shape_manual(values = c(16, 3))
ggsave(filename = "COPD-shaped across cell types.png",
       plot=pl, width = 10, height = 8, dpi = 300)


# FURTHER DATA ANALYSIS

# 1. Differential Expression Analysis (DEG)
# 1.1 Find marker genes for each cell type
# 1.2.1 BASAL
# Create a new column for celltype - sample
ref$celltype.infstat <- paste(ref$celltype.new, ref$sample, sep = "_")
Idents(ref) <- "celltype.infstat"

markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "Basal_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "Basal_infected"),
                       test.use = "wilcox")

# Add 'gene' column to the dataframe 
markers$gene <- rownames(markers)
# Extract the results to an Excel file
write.xlsx(markers, "/Users/ycecy/file/basal.markers.xlsx")

# 1.2.2 CILIATED
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "Ciliated_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "Ciliated_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/ciliated.markers.xlsx")

# 1.2.3 DEUTEROSOMAL
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "Deuterosomal_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "Deuterosomal_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/deuterosomal.markers.xlsx")

# 1.2.4 FIBROBLASTS
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "Fibroblasts_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "Fibroblasts_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/fibroblasts.markers.xlsx")

# 1.2.5 IONOCYTE_n_BRUSH
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "Ionocyte_n_Brush_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "Ionocyte_n_Brush_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/ionocyte_n_brush.markers.xlsx")

# 1.2.6 MYOEPITHELIAL
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "Myoepithelial_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "Myoepithelial_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/myoepithelial.markers.xlsx")

# 1.2.7 SECRETORY
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "Secretory_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "Secretory_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/secretory.markers.xlsx")

# 1.2.8 SMG_DUCT
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "SMG_Duct_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "SMG_Duct_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/SMG_duct.markers.xlsx")

# 1.2.9 SMG_MUCOUS
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "SMG_Mucous_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "SMG_Mucous_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/SMG_mucous.markers.xlsx")

# 1.2.10 SMG_SEROUS
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "SMG_Serous_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "SMG_Serous_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/SMG_serous.markers.xlsx")

# 1.2.11 SUPRABASAL
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.infstat == "Suprabasal_control"),
                       ident.2 = WhichCells(ref, expression = celltype.infstat == "Suprabasal_infected"),
                       test.use = "wilcox")
markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/suprabasal.markers.xlsx")

# 1.2 Create a new column for donor - infection 
ref$donor.treat <- paste(ref$donor, ref$sample, sep = "_")
Idents(ref) <- "donor.treat"

# Find differentially expressed genes for COPD condition of donors
ref.de <- FindMarkers(object = ref, 
                      ident.1 = c("donor1_control", "donor2_control", "donor3_control"), 
                      ident.2 = c("donor1_infected", "donor2_infected", "donor3_infected"),
                      test.use = "wilcox")

ref.de$gene <- rownames(ref.de)
head(ref.de, n = 10)

write.xlsx(ref.de, file = "DEG_results_ref.xlsx", rowames = TRUE)

# 1.3 Find marker genes for the infection status
ref$vload <- ref$inf_status
ref$vload[ref$inf_status %in% c("high", "medium")] <- "high+medium"
table(ref$vload)

ref$celltype.vload <- paste(ref$celltype.new, ref$vload, sep = "_")
Idents(ref) <- "celltype.vload"
table(ref$celltype.vload)

# 1.3.1 BASAL
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "Basal_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "Basal_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_Basal_markers.xlsx")

# 1.3.2 CILIATED
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "Ciliated_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "Ciliated_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_Ciliated_markers.xlsx")

# 1.3.3 DEUTEROSOMAL
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "Deuterosomal_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "Deuterosomal_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_Deuterosomal_markers.xlsx")

# 1.3.4 FIBROBLASTS
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "Fibroblasts_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "Fibroblasts_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_Fibroblasts_markers.xlsx")


# 1.3.5 IONOCYTE_n_BRUSH
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "Ionocyte_n_Brush_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "Ionocyte_n_Brush_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_Ionocyte_n_Brush_markers.xlsx")

# 1.3.6 SECRETORY
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "Secretory_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "Secretory_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_Secretory_markers.xlsx")


# 1.3.7 SMG_DUCT
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "SMG_Duct_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "SMG_Duct_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_SMG_Duct_markers.xlsx")


# 1.3.8 SMG_MUCOUS
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "SMG_Mucous_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "SMG_Mucous_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_SMG_Mucous_markers.xlsx")


# 1.3.9 SMG_SEROUS
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "SMG_Serous_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "SMG_Serous_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_SMG_Serous_markers.xlsx")


# 1.3.10 SUPRABASAL
markers <- FindMarkers(object = ref,
                       ident.1 = WhichCells(ref, expression = celltype.vload == "Suprabasal_high+medium"),
                       ident.2 = WhichCells(ref, expression = celltype.vload == "Suprabasal_uninfected"),
                       test.use = "wilcox")

markers$gene <- rownames(markers)
write.xlsx(markers, "/Users/ycecy/file/inf.degree_Suprabasal_markers.xlsx")


# 2. Gene Ontology Analysis (GO)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


# 2.2 Gene Ontology Analysis over infection degree
options(enrichplot.colours = c("darkorange","yellow"))


# 2.2.1 BASAL
markers <- read.xlsx("/Users/ycecy/file/inf.degree_Basal_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_basal_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_Basal_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.2 CILIATED
markers <- read.xlsx("/Users/ycecy/file/inf.degree_Ciliated_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_ciliated_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_ciliated_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.3 DEUTEROSOMAL
markers <- read.xlsx("/Users/ycecy/file/inf.degree_Deuterosomal_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_deuterosomal_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_deuterosomal_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.4 FIBROBLASTS
markers <- read.xlsx("/Users/ycecy/file/inf.degree_Fibroblasts_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_fibroblasts_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_fibroblasts_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.5 IONOCYTE_n_BRUSH
markers <- read.xlsx("/Users/ycecy/file/inf.degree_Ionocyte_n_Brush_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_ionocyte_n_brush_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_ionocyte_n_brush_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.6 SECRETORY
markers <- read.xlsx("/Users/ycecy/file/inf.degree_Secretory_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_secretory_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_secretory_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.7 SMG_DUCT
markers <- read.xlsx("/Users/ycecy/file/inf.degree_SMG_Duct_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_SMG_duct_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_SMG_duct_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.8 SMG_MUCOUS
markers <- read.xlsx("/Users/ycecy/file/inf.degree_SMG_Mucous_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_SMG_mucous_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_SMG_mucous_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.9 SMG_SEROUS
markers <- read.xlsx("/Users/ycecy/file/inf.degree_SMG_Serous_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_SMG_serous_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_SMG_serous_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# 2.2.10 SUPRABASAL
markers <- read.xlsx("/Users/ycecy/file/inf.degree_Suprabasal_markers.xlsx")

GO_results <- enrichGO(gene = markers$gene,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "inf.degree_suprabasal_GO.xlsx", rowames = TRUE)

pl <- dotplot(GO_results, showCategory = 20)
ggsave(filename = "inf.degree_suprabasal_GO.png",
       plot=pl, width = 9, height = 10, dpi = 300)


# Proportion plot 
library(speckle)

Idents(ref) <- "celltype.new"

ref@meta.data$group <- ifelse(ref@meta.data$sample %in% c("cond_1", "cond_2"), "control", "infected")

# Rename columns
ref@meta.data$sample <- ref@meta.data$donor  
ref@meta.data$cluster <- ref@meta.data$celltype.new 

# Plot the cell type proportions by donor and condition
pl <- plotCellTypeProps(ref)
ggsave(filename = "ProportionPlotCellTypes.png",
       plot=pl, width = 9, height = 7, dpi = 300) 
