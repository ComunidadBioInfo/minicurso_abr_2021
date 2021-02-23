library(tximeta)
library(tximportData)
library(SummarizedExperiment)
library(DESeq2)


dir <- ("../salmon_quants/")
samples <- read.table(file.path(dir,  "metadata.txt"), he = T)
files <- file.path(dir, paste0(samples$Sample, "_quant"), "quant.sf")
data.frame(file = files,  exist = file.exists(files))

coldata <- data.frame(files,  names = samples$Unique_id, condition = samples$Condition, stringsAsFactors = F)
se <- tximeta(coldata)
colData(se)
gse <- summarizeToGene(se)
assayNames(se)
rowRanges(se)
seqinfo(se)
edb <- retrieveDb(se)
class(edb)
gse <- summarizeToGene(se)
rowRanges(gse)

rowRanges(gse)
