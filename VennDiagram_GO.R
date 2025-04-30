# Set the directory
setwd("~/ycecy/file")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

# Load the dataset
ref <- readRDS("/Users/ycecy/file/copd_annotated_new.rds")

# Merge subgroups to bigger ones for Differential Expression Analysis 
ref$celltype.new <- ref$predicted.id

ref$celltype.new[ref$predicted.id %in% c("Basal", "Dividing_Basal")] <- "Basal"
ref$celltype.new[ref$predicted.id %in% c("Secretory_Club", "Secretory_Goblet")] <- "Secretory"
ref$celltype.new[ref$predicted.id %in% c("Fibro_alveolar", "Fibro_adventitial", "Fibro_perichondrial",
                                         "Mesothelia", "Fibro_myofibroblast", "Fibro_peribronchial")] <- "Fibroblasts"
#ref$celltype.new[ref$predicted.id %in% c("SMG_Duct", "SMG_Mucous", "SMG_Serous")] <- "SMG"

table(ref$celltype.new)

# Create metadata for COPD condition of the donors
ref$copd_condition <- ifelse(ref$donor %in% c("donor4", "donor5", "donor6"),
                             "copd+", 
                             "copd-")

# Create metadata for infection status according to the set threshold
# Log transformed viral load
ref$virallog10 <- log10(ref$viralLoad)

# Omit NA values
ref$virallog10[is.infinite(ref$virallog10)] <- NA

# Set the threshold
ref$inf_status <- ifelse(ref$viralLoad == 0, "uninfected",
                         ifelse(ref$viralLoad > 0 & ref$viralLoad <= 0.0001, "low",
                                ifelse(ref$viralLoad > 0.0001 & ref$viralLoad <= 0.001, "medium", "high")))
table(ref$inf_status)

# 2. Gene Ontology (GO) Enrichment Analysis
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(openxlsx)

# 2.1 Create a new column for celltype - sample - copd
ref$celltype.treat.copd <- paste(ref$celltype.new, ref$sample, ref$copd_condition, sep = "_")
Idents(ref) <- "celltype.treat.copd"

# 2.1.1 BASAL
exp_gr_basal <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Basal_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Basal_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Basal_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Basal_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Basal_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Basal_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Basal_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Basal_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
library(VennDiagram)
library(ggVennDiagram)

# Define the  threshold
p_val_adj_threshold <- 0.05

# Top genes based on threshold
top_genes <- lapply(exp_gr_basal, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_Basal.xlsx", rowames = TRUE)

options(enrichplot.colours = c("lightblue","red"))

# Gene Set 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_Basal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_Basal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300)

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_Basal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_Basal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300)

# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_Basal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_Basal_GO.png",
       plot=pl, width = 9, height = 12, dpi = 300)

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_Basal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_Basal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300)

# 2.1.2 CILIATED
exp_gr_ciliated <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Ciliated_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Ciliated_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Ciliated_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Ciliated_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Ciliated_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Ciliated_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Ciliated_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Ciliated_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
# Top genes based on threshold
top_genes <- lapply(exp_gr_ciliated, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_Ciliated.xlsx", rowames = TRUE)


# Extract the genes from areas of interest on Venn diagram
# Gene Set 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_Ciliated_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_Ciliated_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300)

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_Ciliated_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_Ciliated_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300)   


# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_Ciliated_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_Ciliated_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300)   

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_Ciliated_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_Ciliated_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300)   


# 2.1.3 DEUTEROSOMAL
exp_gr_deuterosomal <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Deuterosomal_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Deuterosomal_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Deuterosomal_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Deuterosomal_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Deuterosomal_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Deuterosomal_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Deuterosomal_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Deuterosomal_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
# Top genes based on threshold
top_genes <- lapply(exp_gr_deuterosomal, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_Deuterosomal.xlsx", rowames = TRUE)


# Extract the genes from areas of interest on Venn diagram
# Gene Set 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_Deuterosomal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_Deuterosomal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300)   

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_Deuterosomal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_Deuterosomal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_Deuterosomal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_Deuterosomal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_Deuterosomal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_Deuterosomal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 2.1.4 FIBROBLASTS
exp_gr_fibroblasts <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Fibroblasts_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Fibroblasts_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Fibroblasts_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Fibroblasts_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Fibroblasts_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Fibroblasts_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Fibroblasts_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Fibroblasts_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
# Top genes based on threshold
top_genes <- lapply(exp_gr_fibroblasts, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_Fibroblasts.xlsx", rowames = TRUE)


# Extract the genes from areas of interest on Venn diagram
# Gene Set 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_Fibroblasts_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_Fibroblasts_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_Fibroblasts_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_Fibroblasts_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_Fibroblasts_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_Fibroblasts_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_Fibroblasts_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_Fibroblasts_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# 2.1.5 IONOCYTE_n_BRUSH

# 2.1.6 SECRETORY
exp_gr_secretory <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Secretory_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Secretory_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Secretory_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Secretory_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Secretory_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Secretory_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Secretory_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Secretory_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
# Top genes based on threshold
top_genes <- lapply(exp_gr_secretory, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_Secretory.xlsx", rowames = TRUE)


# Extract the genes from areas of interest on Venn diagram
# Gene Set 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_Secretory_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_Secretory_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_Secretory_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_Secretory_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_Secretory_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_Secretory_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_Secretory_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_Secretory_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 2.1.7 SMG_DUCT
exp_gr_duct <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Duct_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Duct_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Duct_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Duct_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Duct_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Duct_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Duct_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Duct_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
# Top genes based on threshold
top_genes <- lapply(exp_gr_duct, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_SMG_Duct.xlsx", rowames = TRUE)

# Extract the genes from areas of interest on Venn diagram
# Gene Set 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_SMG_Duct_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_SMG_Duct_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_SMG_Duct_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_SMG_Duct_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_SMG_Duct_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_SMG_Duct_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_SMG_Duct_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_SMG_Duct_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 2.1.8 SMG_MUCOUS
exp_gr_mucous <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Mucous_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Mucous_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Mucous_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Mucous_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Mucous_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Mucous_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Mucous_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Mucous_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
# Top genes based on threshold
top_genes <- lapply(exp_gr_mucous, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_SMG_Mucous.xlsx", rowames = TRUE)


# Extract the genes from areas of interest on Venn diagram
# Gene Set 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_SMG_Mucous_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_SMG_Mucous_GO.png",
       plot=pl, width = 9, height = 12, dpi = 300) 

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_SMG_Mucous_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_SMG_Mucous_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_SMG_Mucous_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_SMG_Mucous_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_SMG_Mucous_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_SMG_Mucous_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 2.1.9 SMG_SEROUS
exp_gr_serous <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Serous_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Serous_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Serous_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Serous_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Serous_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Serous_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Serous_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "SMG_Serous_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
# Top genes based on threshold
top_genes <- lapply(exp_gr_serous, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_SMG_Serous.xlsx", rowames = TRUE)


# Extract the genes from areas of interest on Venn diagram
# Combination 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_SMG_Serous_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_SMG_Serous_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_SMG_Serous_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_SMG_Serous_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_SMG_Serous_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_SMG_Serous_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_SMG_Serous_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_SMG_Serous_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 2.1.10 SUPRABASAL
exp_gr_suprabasal <- list(
  "cntrl/copd- vs copd+" = FindMarkers(ref,
                                       ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Suprabasal_control_copd-"),
                                       ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Suprabasal_control_copd+"),
                                       test.use = "wilcox"),
  "inf/copd- vs copd+" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Suprabasal_infected_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Suprabasal_infected_copd+"), 
                                     test.use = "wilcox"),
  "copd-/cntrl vs inf" = FindMarkers(ref, 
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Suprabasal_control_copd-"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Suprabasal_infected_copd-"), 
                                     test.use = "wilcox"),
  "copd+/cntrl vs inf" = FindMarkers(ref,
                                     ident.1 = WhichCells(ref, expression = celltype.treat.copd == "Suprabasal_control_copd+"), 
                                     ident.2 = WhichCells(ref, expression = celltype.treat.copd == "Suprabasal_infected_copd+"), 
                                     test.use = "wilcox"))

# Venn diagram for overlaps of interested combinations
# Define the  threshold
p_val_adj_threshold <- 0.05

# Top genes based on threshold
top_genes <- lapply(exp_gr_suprabasal, function(x) rownames(x[x$p_val_adj < p_val_adj_threshold, ]))

# Create the diagram
venn <- venn.diagram(
  x = top_genes,
  filename = NULL,
  col = "black",
  fill = c("lightpink", "lightblue", "lightgreen", "lightyellow"),
  alpha = 0.5,
  cex = 0.6,
  cat.cex = 0.7)

# Calculate the overlaps on the diagram
overlaps <- calculate.overlap(top_genes)

# Custom Labels for regions (overlaps) on the diagram
labels <- list(
  "A" = overlaps$a1,  
  "B" = overlaps$a2,   
  "C" = overlaps$a3,   
  "D" = overlaps$a4,
  "E" = overlaps$a5,
  "F" = overlaps$a6,
  "G" = overlaps$a7,
  "H" = overlaps$a8,
  "I" = overlaps$a9,
  "J" = overlaps$a10,
  "K" = overlaps$a11,
  "L" = overlaps$a12,
  "M" = overlaps$a13,
  "N" = overlaps$a14,
  "O" = overlaps$a15)

# Annotate the diagram with custom labels
grid.newpage()
grid.draw(venn)

grid.text("A", x = 0.4, y = 0.75)
grid.text("B", x = 0.5, y = 0.62)
grid.text("C", x = 0.6, y = 0.75)
grid.text("D", x = 0.32, y = 0.65)
grid.text("E", x = 0.38, y = 0.5)
grid.text("F", x = 0.5, y = 0.5)
grid.text("G", x = 0.62, y = 0.5)
grid.text("H", x = 0.68, y = 0.65)
grid.text("I", x = 0.15, y = 0.5)
grid.text("J", x = 0.3, y = 0.42)
grid.text("K", x = 0.38, y = 0.32)
grid.text("L", x = 0.62, y = 0.32)
grid.text("M", x = 0.7, y = 0.42)
grid.text("N", x = 0.85, y = 0.5)
grid.text("O", x = 0.5, y = 0.15)

write.xlsx(labels, file = "labels_Suprabasal.xlsx", rowames = TRUE)


# Extract the genes from areas of interest on Venn diagram
# Combination 1: 
genes <- labels$B
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "B_Suprabasal_GO.xlsx", rowames = TRUE)

# Bar plot of the enriched genes
pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "B_Suprabasal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
genes <- c(labels$C, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CH_Suprabasal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CH_Suprabasal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
genes <- c(labels$H, labels$N)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "HN_Suprabasal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "HN_Suprabasal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
genes <- c(labels$C, labels$N, labels$H)
GO_results <- enrichGO(gene = genes,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

GO_df <- as.data.frame(GO_results)
write.xlsx(GO_df, file = "CNH_Suprabasal_GO.xlsx", rowames = TRUE)

pl <- barplot(GO_results, showCategory = 20)
ggsave(filename = "CNH_Suprabasal_GO.png",
       plot=pl, width = 9, height = 9, dpi = 300) 
