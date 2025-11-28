
##############################
# SECTION 1 — INSTALL PACKAGES
##############################

# Install BiocManager if missing
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor packages
bioc_packages <- c("minfi", "GEOquery", "R.utils")

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

# CRAN packages
cran_packages <- c("readr", "data.table", "dplyr", "ggplot2")

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# Load libraries
library(minfi)
library(GEOquery)
library(R.utils)
library(dplyr)
library(ggplot2)
library(data.table)
library(readr)


###############################################
# SECTION 2 — DOWNLOAD GEO SUPPLEMENTARY FILES
###############################################

getwd()

# Download IDAT files
getGEOSuppFiles("GSE68777")

# Extract the RAW tar file
untar("GSE68777/GSE68777_RAW.tar",
      exdir = "GSE68777/idat")

# View extracted IDAT files
head(list.files("GSE68777/idat", pattern = "idat"))


###################################################
# SECTION 3 — DECOMPRESS IDAT FILES AND READ THEM
###################################################

# List all .idat.gz files
idatFiles <- list.files("GSE68777/idat",
                        pattern = "idat.gz$",
                        full.names = TRUE)

# Unzip each IDAT file
sapply(idatFiles, gunzip, overwrite = TRUE)

# Read IDAT files into RGChannelSet
rgSet <- read.metharray.exp("GSE68777/idat")
rgSet

# View sample names and metadata
pData(rgSet)
sampleNames(rgSet)


##########################################################
# SECTION 4 — LOAD SERIES MATRIX (METADATA FROM GEO PAGE)
##########################################################
getwd()
# Download metadata
geoMat <- getGEO("GSE68777")

# Decompress matrix file
gunzip("GSE68777/GSE68777_series_matrix.txt.gz",
       overwrite = TRUE)

# Load phenotype data
geoMat <- getGEO(
  filename = "GSE68777/GSE68777_series_matrix.txt",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

# Extract pData
pD.all <- pData(geoMat)

# Select relevant columns
pD <- pD.all[, c("title", "geo_accession",
                 "characteristics_ch1.1",
                 "characteristics_ch1.2")]
pD

# Clean and rename metadata columns
names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex   <- sub("^Sex: ", "", pD$sex)


###############################################################
# SECTION 5 — MATCH GEO METADATA WITH IDAT SAMPLE NAMES
###############################################################

# Fix sample names to match
sampleNames(rgSet) <- sub(".*_5", "5",
                          sampleNames(rgSet))

# Set row names for phenotype table
rownames(pD) <- pD$title

# Reorder phenotype table to match IDAT sample order
pD <- pD[sampleNames(rgSet), ]


########################################################################
# SECTION 6 — INSTALL ILLUMINA 450K ANNOTATION & ARRAY MANIFEST PACKAGES
########################################################################

BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("lumi")

library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(lumi)


##############################################
# SECTION 7 — PREPROCESSING & QUALITY CONTROL
##############################################

# Convert raw to methylation set
mSet <- preprocessRaw(rgSet)

# Calculate detection p-values
detP <- detectionP(rgSet)

# Check failed probes per sample
failed <- colMeans(detP > 0.01)
failed

# Remove probes failing in >10% samples
keep <- rowMeans(detP < 0.01) > 0.9
rgSet <- rgSet[keep, ]

# Normalize (NOOB)
mSet.norm <- preprocessNoob(rgSet)

# Extract Beta and M-values
betaVals <- getBeta(mSet.norm)
mVals <- getM(mSet.norm)

betaVals
mVals


############################
# SECTION 8 — PCA ANALYSIS
############################

# Randomly sample 5000 probes for PCA
beta_small <- betaVals[sample(1:nrow(betaVals), 5000), ]

# PCA on transposed matrix (samples as rows)
pca <- prcomp(t(beta_small), scale = TRUE)

pca

# Plot PCA
par(mar=c(4,4,2,1))
plot(pca$x[,1], pca$x[,2],
     col = as.factor(pD$group),
     pch = 19,
     main = "PCA — Methylation Data")

legend("topright", legend = unique(pD$group),
       col = 1:length(unique(pD$group)), pch = 19)


##########################################################
# SECTION 9 — DIFFERENTIAL METHYLATION ANALYSIS (LIMMA)
##########################################################

library(limma)

# Convert AnnotatedDataFrame to standard data frame
pD <- pD@data

# Encode group and sex as factors
pD$group <- factor(pD$group)
pD$sex   <- factor(pD$sex)

# Design matrix (group comparison)
design <- model.matrix(~ group, data = pD)
design

# Fit linear model
fit <- lmFit(mVals, design)
fit <- eBayes(fit)

fit

# Extract results (treatment vs control)
results <- topTable(fit, coef = 2, number = Inf)
head(results)


#############################################
# SECTION 10 — ANNOTATE CpGs WITH GENOMIC INFO
#############################################

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Add probe IDs as a column
results$CpG_ID <- rownames(results)
ann450k$CpG_ID <- rownames(ann450k)

# Merge annotation with results
results_annot <- merge(results, ann450k,
                       by = "CpG_ID",
                       all.x = TRUE)

head(results_annot)

#############################
# SECTION 11 — VOLCANO PLOT
#############################

volcano <- data.frame(
  logFC = results$logFC,
  negLogP = -log10(results$P.Value),
  significant = results$adj.P.Val < 0.05
)

ggplot(volcano, aes(x = logFC, y = negLogP, color = significant)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Methylated Positions",
       x = "log2 Fold Change (M-values)",
       y = "-log10 P-value") +
  scale_color_manual(values = c("grey", "red"))


########################################
# SECTION 12 — HEATMAP OF TOP 100 CpGs
########################################

install.packages("pheatmap")
library(pheatmap)

# Select top 100 CpGs
top100 <- head(rownames(results), 100)
heat_data <- betaVals[top100, ]

# Sample annotation for heatmap
sample_ann <- data.frame(Group = pD$group)
rownames(sample_ann) <- rownames(pD)

# Plot heatmap
pheatmap(
  heat_data,
  annotation_col = sample_ann,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_rownames = FALSE,
  main = "Top 100 Differentially Methylated CpGs"
)

###############################################################
# END OF FULL PIPELINE
###############################################################




# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68777/
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68777
  
