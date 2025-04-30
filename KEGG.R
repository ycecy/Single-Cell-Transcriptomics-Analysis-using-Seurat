# Set the directory
setwd("~/ycecy/file")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(stringr)

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

# 4. KEGG Enrichment Analysis
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(openxlsx)
library(pathview)

# 4.1 BASAL
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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.1.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                organism     = 'hsa',
                pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

options(enrichplot.colours = c("pink","darkblue"))

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Basal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.1.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
              organism     = 'hsa',
              minGSSize    = 3,
              maxGSSize    = 800,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH")

head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "1_Basal_gseKEGG.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Basal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.1.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "2_Basal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.1.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
                   organism     = 'hsa',
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH")

head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "2_Basal_gseKEGG.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Basal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.1.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Basal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.1.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
                   organism     = 'hsa',
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH")

head(gse)

# Save the result tables as excel file
write.xlsx(gsekegg, file = "3_Basal_gseKEGG.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gsekeggEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Basal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.1.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways" )
ggsave(filename = "4_Basal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.1.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH")

head(gse)

# Save the result tables as excel file
write.xlsx(gse, file = "4_Basal_gseKEGG.xlsx", rowNames = TRUE)

gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Basal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.2 CILIATED
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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.2.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Ciliated_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.2.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "1_Ciliated_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Ciliated_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.2.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "2_Ciliated_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.2.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "2_Ciliated_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Ciliated_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.2.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Ciliated_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.2.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gsekegg, file = "3_Ciliated_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Ciliated_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.2.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Ciliated_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.2.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "4_Ciliated_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Ciliated_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.3 DEUTEROSOMAL
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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.3.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Deuterosomal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.3.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "1_Deuterosomal_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Deuterosomal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.3.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways" )
ggsave(filename = "2_Deuterosomal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.3.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "2_Deuterosomal_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Deuterosomal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.3.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Deuterosomal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.3.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "3_Deuterosomal_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Deuterosomal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.3.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Deuterosomal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.3.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "4_Deuterosomal_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Deuterosomal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.4.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways") 
ggsave(filename = "1_Fibroblasts_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.4.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "1_Fibroblasts_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Fibroblasts_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.4.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "2_Fibroblasts_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.4.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "2_Fibroblasts_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Fibroblasts_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.4.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Fibroblasts_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.4.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gsekegg, file = "3_Fibroblasts_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gsekeggEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Fibroblasts_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.4.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Fibroblasts_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.4.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "4_Fibroblasts_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Fibroblasts_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# 3.5 IONOCYTE_n_BRUSH

# 3.6 SECRETORY
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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.6.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Secretory_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.6.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "1_Secretory_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Secretory_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.6.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "2_Secretory_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.6.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "2_Secretory_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Secretory_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.6.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Secretory_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.6.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "3_Secretory_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Secretory_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.6.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Secretory_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.6.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "4_Secretory_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Secretory_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# 3.7 SMG_DUCT
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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.7.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Duct_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.7.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "1_Duct_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Duct_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.7.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "2_Duct_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.7.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "2_Duct_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Duct_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.7.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Duct_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.7.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "3_Duct_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Duct_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.7.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Duct_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.7.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "4_Duct_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Duct_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.8 SMG_MUCOUS
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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.7.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Mucous_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.8.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "1_Mucous_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Mucous_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.8.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "2_Mucous_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.8.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "2_Mucous_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Mucous_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.7.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways" )
ggsave(filename = "3_Mucous_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.8.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "3_Mucous_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Mucous_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.8.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Mucous_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.8.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "4_Mucous_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Mucous_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.9 SMG_SEROUS
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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.9.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Serous_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.9.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "1_Serous_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "1_Serous_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.9.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "2_Serous_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.9.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "2_Serous_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Serous_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.9.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 30, title = "Enriched Pathways")
ggsave(filename = "3_Serous_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.9.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gsekegg, file = "3_Serous_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)


pl <- dotplot(gsekeggEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "3_Serous_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.9.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Serous_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.9.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "4_Serous_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Serous_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.10 SUPRABASAL
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

# Only keep the genes that have Ensembl ID
ens.A <- A_gr[str_detect(A_gr$gene, "^ENSG"), ]
ens.B <- B_gr[str_detect(B_gr$gene, "^ENSG"), ]
ens.C <- C_gr[str_detect(C_gr$gene, "^ENSG"), ]
ens.D <- D_gr[str_detect(D_gr$gene, "^ENSG"), ]

# Extract the log fold change of the correspond genes for gseKEGG
GenelistA <- ens.A$avg_log2FC
geneListB <- ens.B$avg_log2FC
GeneListC <- ens.C$avg_log2FC
genelistD <- ens.D$avg_log2FC

names(GenelistA) <- ens.A$gene
names(geneListB) <- ens.B$gene
names(GeneListC) <- ens.C$gene
names(genelistD) <- ens.D$gene

# Gene Set 1: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GenelistA),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
Genelist <- GenelistA[names(gene) %in% names(GenelistA)]
names(Genelist) <- gene

# Sort geneList in decreasing order
Genelist <- sort(Genelist, decreasing = TRUE)

# 4.10.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene       = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Suprabasal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.10.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList = Genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "1_Suprabasal_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "1_Suprabasal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# Gene Set 2: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(geneListB),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
geneList <- geneListB[names(gene) %in% names(geneListB)]
names(geneList) <- gene

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# 4.10.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "2_Suprabasal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.10.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "2_Suprabasal_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "2_Suprabasal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 3: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(GeneListC),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
GeneList <- GeneListC[names(gene) %in% names(GeneListC)]
names(GeneList) <- gene

# Sort geneList in decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

# 4.10.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Suprabasal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.10.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = GeneList,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "3_Suprabasal_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "3_Suprabasal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 

# Gene Set 4: 
# Convert Ensembl IDs to Entrez IDs (required for gseKEGG)
gene <- mapIds(org.Hs.eg.db, 
               keys = names(genelistD),
               column = "ENTREZID", 
               keytype = "ENSEMBL", 
               multiVals = "first")

head(gene)

# Replace Ensembl IDs with Entrez IDs in geneList
genelist <- genelistD[names(gene) %in% names(genelistD)]
names(genelist) <- gene

# Sort geneList in decreasing order
genelist <- sort(genelist, decreasing = TRUE)

# 4.10.1 KEGG Pathway Over-Representation Analysis
enr <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(enr)

enrEdited <- enr
enrEdited@result$Description <- sapply(strsplit(enrEdited@result$Description, " - "), "[", 1)

# Dot plot
pl <- dotplot(enrEdited, showCategory = 20, title = "Enriched Pathways")
ggsave(filename = "4_Suprabasal_enrKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 


# 4.10.2 KEGG Pathway Gene Set Enrichment Analysis
gse <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")

head(gse)

write.xlsx(gse, file = "4_Suprabasal_gseKEGG.xlsx", rowNames = TRUE)
gseEdited <- gse
gseEdited@result$Description <- sapply(strsplit(gseEdited@result$Description, " - "), "[", 1)

pl <- dotplot(gseEdited, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(filename = "4_Suprabasal_gseKEGG.png",
       plot=pl, width = 9, height = 9, dpi = 300) 
