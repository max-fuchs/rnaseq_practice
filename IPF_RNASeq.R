library(GEOquery)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(gplots)
library(pheatmap)
library(EnhancedVolcano)
library(EDASeq)
library(biomaRt)
library(tibble)
library(writexl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(ReactomePA)

setwd(choose.dir()) ##choose the directory where you files should be saved

gse <- getGEO("GSE92592", GSEMatrix = TRUE)

getGEOSuppFiles("GSE92592", makeDirectory = TRUE)

#path to zip files
filename <- "GSE92592/GSE92592_gene.counts.txt.gz"

#reading in count table
counts <- read.table(gzfile(filename), header = TRUE, sep = "\t")

#extracting sample order from gse object
order <- gse[["GSE92592_series_matrix.txt.gz"]]@phenoData@data[["title"]]

#extracting order from count object
current_order <- colnames(counts)

#find indices for reordering
order_indices <- match(order, current_order)

#reordering
counts_ordered <- counts[, order_indices]

#extract coldata from gse
coldata <- gse[["GSE92592_series_matrix.txt.gz"]]@phenoData@data

# Rename the column
colnames(coldata)[colnames(coldata) == "disease state:ch1"] <- "disease"

##Differential gene expression analysis
##Checking annotation
all(rownames(coldata) %in% colnames(counts_ordered))
all(rownames(coldata) == colnames(counts_ordered))

row.names(coldata) <- coldata$title

##Creating DESeq2dataset object

dds <- DESeqDataSetFromMatrix(countData = counts_ordered,
                              colData = coldata,
                              design = ~ disease)

#Pre-filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= nrow(coldata)
dds <- dds[keep,]
nrow(dds)

#rlog transformation
rld <- rlog(dds, blind = F)
head(assay(rld), 3)

#PCA
pcaData <- plotPCA(rld,intgroup = c("disease","title"),returnData = T)
pcaData

percentVar <- round(100*attr(pcaData, 'percentVar'))
ggplot(pcaData) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = disease), size = 5) +
  xlab(paste0('PC1: ', percentVar[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar[2], '% variance')) +
  coord_fixed() +
  geom_text(mapping = aes(x= PC1, y = PC2, label = title), nudge_y = -1)

#Differential expression analysis
dds <- DESeq(dds)
plotDispEsts(dds)
plotRLE(counts(dds))
plotRLE(counts(dds, normalized=T))

#create result objects
IPF_vs_healthy <- results(dds, contrast = c("disease", "Idiopathic Pulmonary Fibrosis", "Control"), pAdjustMethod = "BH")

summary(IPF_vs_healthy)

head(IPF_vs_healthy)

#heatmap of the count matrix
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:200]

pheatmap(assay(rld)[select, ],
         cluster_rows = TRUE,
         show_rownames = FALSE,
         show_colnames = F,
         cluster_cols = F,
         annotation_col = coldata[38],
         color = bluered(100),
         scale = "row")
##plot heatmap degs
IPF_vs_healthy_sig <- as.data.frame(IPF_vs_healthy) %>% dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1.5)

#only keep DEGs
matrix <- assay(rld)[rownames(IPF_vs_healthy_sig)[1:250],]

pheatmap(matrix, scale = 'row',
         show_rownames = F,
         show_colnames = F,
         color = bluered(100),
         fontsize = 8,
         annotation_col = coldata[38],
         cluster_cols = T)

#volcano plot
#names <- read.csv("DEG_symbols_for_volcano.csv", header = F)
#names <- names$V1

EnhancedVolcano(IPF_vs_healthy,
                lab = rownames(IPF_vs_healthy),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Volcano plot: IPF vs Healthy',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0,
                col = c('grey30', 'grey30', 'grey30', 'orange'),
                colAlpha = 0.75,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE)

# Update row names by removing chromosomal info
rownames(IPF_vs_healthy) <- gsub("\\.chr[0-9XY]+", "", rownames(IPF_vs_healthy))
#Extract DEG list
export <- as.data.frame(IPF_vs_healthy)
export$gene_symbol <- row.names(export)
write_xlsx(export, "IPF_vs_healthy.xlsx")


# Assuming 'IPF_vs_healthy' contains your DESeq2 results
# Create a ranked list of genes
geneList <- IPF_vs_healthy$log2FoldChange
names(geneList) <- rownames(IPF_vs_healthy, )

# Sort the gene list based on log2 fold changes
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[abs(geneList) > 3]

# If needed, map your gene symbols to ENTREZ IDs
geneList_mapped <- bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# Keep only those genes that could be mapped
valid_genes <- !is.na(geneList_mapped$ENTREZID)
geneList_filtered <- geneList[names(geneList) %in% geneList_mapped$SYMBOL[valid_genes]]
names(geneList_filtered) <- geneList_mapped$ENTREZID[valid_genes]

x <- enrichPathway(names(geneList_filtered), pvalueCutoff = 0.05, readable=TRUE)

barplot(x, showCategory=20)
