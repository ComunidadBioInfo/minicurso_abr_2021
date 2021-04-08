DEGResults <- function(qlf) {
  ##This function returns a dataframe with all DEG
  ##qlf: Object obatined from the generalized linear model
  qlf <- topTags(qlf, n = Inf)
  qlf <- as.data.frame(qlf) 
  return(qlf)
}

volcanoplotR <- function(dge.obj, logfc, p.adj, type) {
  ##This function adds a new column (T or F) according to the FDR and LFC of each gene in edgeR list of DEG
  ##dge.obj: List with DEG
  ##logFC: logFC threshold used for the differential expression test
  ##p.adj: p.adj or FDR threshold to obtain significant genes
  ##type: Type of the output "edgeR" or "DESeq2"
  ##Updated 5-mar-2021 Rodolfo Chavez
  if(type == "edgeR") {
    volc <- dge.obj %>%
      mutate(condition = ifelse((dge.obj$logFC > logfc) & (dge.obj$FDR < p.adj), "Over-expressed",
                                ifelse((dge.obj$logFC < -logfc) & (dge.obj$FDR < p.adj), "Sub-expressed",
                                       ifelse((dge.obj$logFC > logfc) & (dge.obj$FDR > p.adj), "NS",
                                              ifelse((dge.obj$logFC < -logfc) & (dge.obj$FDR > p.adj), "NS",
                                                     ifelse((dge.obj$logFC < logfc) & (dge.obj$FDR > p.adj), "NS", "NS"))))))
    volcano_plot <- ggplot(volc)+
      geom_point(aes(x = logFC, y = -log10(FDR), color = condition))+
      scale_color_manual(name = "Condition",
                         labels = paste(c("NS", "Over-expressed", "Sub-expressed"), c(sum(volc$condition == "NS"), sum(volc$condition == "Over-expressed"), sum(volc$condition == "Sub-expressed"))),
                         values = c("#6e6d6e","#d84b47","#66c343"))+
      geom_vline(aes(xintercept = logfc), linetype = "dashed")+
      geom_vline(aes(xintercept = -logfc), linetype = "dashed")+
      geom_hline(aes(yintercept = -log10(p.adj)), linetype = "dashed")+
      theme_set(theme_bw())+
      theme(plot.title = element_text(face = "bold", size = 18), 
            axis.title = element_text(size = 18),
            legend.title = element_text(face = "bold", size = 15),
            legend.text = element_text(size = 15), 
            legend.position = "bottom")
  } else {
    volc <- dge.obj %>%
      mutate(condition = ifelse((dge.obj$log2FoldChange > logfc) & (dge.obj$padj < p.adj), "Over-expressed",
                                ifelse((dge.obj$log2FoldChange < -logfc) & (dge.obj$padj < p.adj), "Sub-expressed",
                                       ifelse((dge.obj$log2FoldChange > logfc) & (dge.obj$padj > p.adj), "NS",
                                              ifelse((dge.obj$log2FoldChange < -logfc) & (dge.obj$padj > p.adj), "NS",
                                                     ifelse((dge.obj$log2FoldChange < logfc) & (dge.obj$padj > p.adj), "NS", "NS")))))) %>%
      drop_na()
    volcano_plot <- ggplot(volc)+
      geom_point(aes(x = log2FoldChange, y = -log10(padj), color = condition))+
      scale_color_manual(name = "Condition",
                         labels = paste(c("NS", "Over-expressed", "Sub-expressed"), c(sum(volc$condition == "NS"), sum(volc$condition == "Over-expressed"), sum(volc$condition == "Sub-expressed"))),
                         values = c("#6e6d6e","#d84b47","#66c343"))+
      geom_vline(aes(xintercept = logfc), linetype = "dashed")+
      geom_vline(aes(xintercept = -logfc), linetype = "dashed")+
      geom_hline(aes(yintercept = -log10(p.adj)), linetype = "dashed")+
      theme_set(theme_bw())+
      theme(plot.title = element_text(face = "bold", size = 18), 
            axis.title = element_text(size = 18),
            legend.title = element_text(face = "bold", size = 15),
            legend.text = element_text(size = 15), 
            legend.position = "bottom")
  }
  return(volcano_plot)
}
