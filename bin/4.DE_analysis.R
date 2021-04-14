library(tximeta)
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(PCAtools)
library(pheatmap)
source("../bin/functions.R")

coldata <- read.delim(here::here("output/salmon_quants/metadata.txt"))
coldata <- coldata %>% dplyr::rename(names = unique_id)
rownames(coldata) <- coldata$names
dir <- file.path(here::here("output/salmon_quants"))
list.files(dir) 
coldata$files <- file.path(dir,paste0(coldata$Sample,"_quant"),"quant.sf")
data.frame(coldata$Sample, file.exists(coldata$files))

se <- tximeta(coldata)
se$Condition <- factor(se$Treatment)
table(se$Condition)
dds <- DESeqDataSet(se, design = ~ Condition)
keep <- rowSums(counts(dds) >= 3) >= 6
dds <- dds[keep, ]
dds <- DESeq(dds)

rld <- rlog(dds, blind = F)
PCA <- pca(assay(rld), scale = T, metadata = coldata) 

res <- results(dds)
res_lfc <- lfcShrink(dds, coef = "Condition_Verafinib_vs_Control",  type = "apeglm")
group <- coldata$Treatment
dds$group <- group
DEG <- as.data.frame(res)
significant <- DEG %>% filter(log2FoldChange > 0 & padj < 0.1 |
                                log2FoldChange < 0 & padj < 0.1)
norm_counts <- counts(dds, normalized = T)
annotation_col <- coldata[, c(2, 3)]