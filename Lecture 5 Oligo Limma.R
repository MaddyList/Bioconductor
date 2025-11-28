#SECTION 1 — SETUP AND INSTALLATION

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


bioc_packages <- c("oligo", "GEOquery", "limma")

for (pkg in bioc_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}


cran_packages <- c("readr", "data.table", "dplyr", "ggplot2")

for (pkg in cran_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

library("oligo")
library("GEOquery")
library("limma")
library("DESeq2")
library("edgeR")
library("airway")
library("leukemiasEset")
library("readr")
library("data.table")
library("dplyr")
library("ggplot2")




getGEOSuppFiles("GSE38792")






list.files("GSE38792")



untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")





celfiles <- list.files("GSE38792/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)

rawData

exprs(rawData)[1:4, 1:3]


#Cleaning Sample Names and Metadata
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename



sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames

pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames), 
                               "OSA", "Control")

pData(rawData)
par(mar=c(4,4,2,1))
boxplot(rawData, target="probeset")

normData <- rma(rawData)

normData
par(mar=c(4,4,2,1))
boxplot(normData)



###“Differential Gene Expression Analysis Using limma”
library(oligo)     # normalized data from CEL files
library(limma)     # differential expression


normData       # ExpressionFeatureSet from oligo
pData(normData)  # sample information


expr <- exprs(normData)


group <- factor(pData(normData)$group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design


fit <- lmFit(expr, design)


expression ~ group


contrast.matrix <- makeContrasts(OSA - Control, levels = design)
contrast.matrix


fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)


results <- topTable(fit2, adjust = "BH", number = Inf)
head(results)


# logFC > 0 → upregulated in OSA
#	logFC < 0 → downregulated in OSA
#	adj.P.Val < 0.05 → significant


volcanoplot(fit2)


library(pheatmap)
topgenes <- rownames(results)[1:50]
pheatmap(expr[topgenes, ], scale = "row")








