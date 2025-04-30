# 3. Gene Set Enrichment Analysis (GSEA)
library(goseq)
library(stringr)
library(europepmc)
library(cowplot)
library(enrichplot)

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

# 3.1 BASAL
# Create a new column for celltype - sample - copd
ref$celltype.treat.copd <- paste(ref$celltype.new, ref$sample, ref$copd_condition, sep = "_")
Idents(ref) <- "celltype.treat.copd"

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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_basal$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_basal$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_basal$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_basal$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Basal_GSEA.xlsx", rowNames = TRUE)

# Store the original data
gseEdited <- gse

# Manipulate the data for shorter term names on plots
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

options(enrichplot.colours = c("red","darkgreen"))

# Dot plot
pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Basal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Basal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Basal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Basal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Basal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Basal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Basal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 3.2 CILIATED
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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_ciliated$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_ciliated$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_ciliated$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_ciliated$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Ciliated_GSEA.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Ciliated_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Ciliated_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Ciliated_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Ciliated_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Ciliated_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Ciliated_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Ciliated_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 3.3 DEUTEROSOMAL
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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_deuterosomal$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_deuterosomal$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_deuterosomal$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_deuterosomal$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Deuterosomal_GSEA.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Deuterosomal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Deuterosomal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Deuterosomal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Deuterosomal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Deuterosomal_GSEA.png",
       plot=pl, width = 9, height = 11, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Deuterosomal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Deuterosomal_GSEA.png",
       plot=pl, width = 9, height = 11, dpi = 300) 


# 3.4 FIBROBLASTS
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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_fibroblasts$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_fibroblasts$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_fibroblasts$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_fibroblasts$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Fibroblasts_GSEA.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Fibroblasts_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Fibroblasts_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Fibroblasts_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Fibroblasts_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Fibroblasts_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Fibroblasts_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Fibroblasts_GSEA.png",
       plot=pl, width = 9, height = 12, dpi = 300) 


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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_secretory$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_secretory$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_secretory$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_secretory$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Secretory_GSEA.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Secretory_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Secretory_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Secretory_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Secretory_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Secretory_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Secretory_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Secretory_GSEA.png",
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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_duct$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_duct$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_duct$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_duct$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Duct_GSEA.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Duct_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Duct_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Duct_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Duct_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Duct_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Duct_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Duct_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 3.8 SMG_MUCOUS
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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_mucous$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_mucous$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_mucous$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_mucous$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Mucous_GSEA.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Mucous_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Mucous_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Mucous_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Mucous_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Mucous_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Mucous_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Mucous_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 3.9 SMG_SEROUS
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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_serous$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_serous$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_serous$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_serous$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Serous_GSEA.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Serous_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Serous_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Serous_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Serous_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Serous_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Serous_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Serous_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 3.10 SUPRABASAL
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

# Extract the list of genes, for all four experimental groups, respectively 
A_gr <- as.data.frame(exp_gr_suprabasal$`cntrl/copd- vs copd+`)
B_gr <- as.data.frame(exp_gr_suprabasal$`inf/copd- vs copd+`)
C_gr <- as.data.frame(exp_gr_suprabasal$`copd-/cntrl vs inf`)
D_gr <- as.data.frame(exp_gr_suprabasal$`copd+/cntrl vs inf`)

A_gr$gene <- rownames(A_gr)
B_gr$gene <- rownames(B_gr)
C_gr$gene <- rownames(C_gr)
D_gr$gene <- rownames(D_gr)

# Only keep the genes that have Ensembl ID (gseGO accepts only one type of gene name)
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseGO
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Sort the vector in decreasing order (important for gseGO)
GenelistA <- sort(GenelistA, decreasing = TRUE)
geneListB <- sort(geneListB, decreasing = TRUE)
GeneListC <- sort(GeneListC, decreasing = TRUE)
genelistD <- sort(genelistD, decreasing = TRUE)

# Gene Set 1: 
gse <- gseGO(geneList = GenelistA, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Suprabasal_GSEA.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Suprabasal_GSEA.png",
       plot=pl, width = 9, height = 14, dpi = 300) 

# Gene Set 2: 
gse <- gseGO(geneList = geneListB, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "2_Suprabasal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Suprabasal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
gse <- gseGO(geneList = GeneListC, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "3_Suprabasal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Suprabasal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
gse <- gseGO(geneList = genelistD, 
             ont ="BP", 
             keyType = "ENSEMBL",
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",
             eps = 0,
             pAdjustMethod = "BH")
head(gse)

write.xlsx(gse, file = "4_Suprabasal_GSEA.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Suprabasal_GSEA.png",
       plot=pl, width = 9, height = 9, dpi = 300) 
