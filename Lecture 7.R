###############################################################
#   COUNT-BASED RNA-seq ANALYSIS USING edgeR & DESeq2 IN R
###############################################################


##############################
# SECTION 1 — INSTALL PACKAGES
##############################

# BiocManager installs Bioconductor packages.
# This block checks whether BiocManager is present;
# if not, it installs it from CRAN.

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define required Bioconductor packages.
# DESeq2 = differential expression analysis
# edgeR  = count-based analysis using negative binomial modeling
# airway = example RNA-seq dataset

bioc_packages <- c("DESeq2", "edgeR", "airway")

# Install missing packages
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load all libraries
library(DESeq2)
library(edgeR)
library(airway)


##############################
# SECTION 2 — LOAD DATA
##############################

# Load the airway dataset (8 human airway smooth muscle samples).
data(airway)


# Explore the dataset
airway


# This prints the SummarizedExperiment object containing:
# - count matrix
# - sample metadata
# - annotation / ranges


##############################
# EXAMINE RAW COUNT DATA
##############################

# View the first 3 genes × first 3 samples
assay(airway, "counts")[1:3, 1:3]




# The experimental variable — dexamethasone treatment
airway$dex



# “untrt” = untreated
# “trt”   = treated


# Set "untrt" as reference level for comparisons
airway$dex <- relevel(airway$dex, "untrt")

# Show genomic ranges for each gene
granges(airway)



###############################################################
# SECTION 3 — edgeR ANALYSIS
###############################################################

################################################
# STEP 1 — Create edgeR DGEList object
################################################
dge <- DGEList(
  counts = assay(airway, "counts"),   # raw count matrix
  group  = airway$dex                 # treatment groups
)

dge

# Merge sample info into DGEList
dge$samples <- merge(
  dge$samples,
  as.data.frame(colData(airway)),
  by = 0
)

# Add gene names to dge$genes
dge$genes <- data.frame(name = names(rowRanges(airway)))

dge$genes

################################################
# STEP 2 — Normalize counts (TMM normalization)
################################################
# calcNormFactors() scales samples to remove library size bias.
dge <- calcNormFactors(dge)


################################################
# STEP 3 — Create Design Matrix
################################################
# model.matrix() encodes the experimental design for GLM.
design <- model.matrix(~ dge$samples$group)
design


################################################
# STEP 4 — Estimate Dispersion
################################################
# edgeR models overdispersion using negative binomial statistics.
dge <- estimateGLMCommonDisp(dge, design)   # global dispersion
dge <- estimateGLMTagwiseDisp(dge, design)  # gene-specific dispersion


################################################
# STEP 5 — Fit the GLM & Perform Likelihood Ratio Test
################################################
fit <- glmFit(dge, design)         # fit NB-GLM for each gene
lrt <- glmLRT(fit, coef = 2)       # compare treated vs untreated

# View top differentially expressed genes
topTags(lrt)


###############################################################
# SECTION 4 — DESeq2 ANALYSIS
###############################################################

################################################
# STEP 1 — Create DESeq2 Dataset
################################################
# DESeqDataSet = count matrix + metadata + design formula
dds <- DESeqDataSet(
  airway,
  design = ~ dex
)


################################################
# STEP 2 — Run DESeq Pipeline
################################################
# DESeq() performs:
# - size factor estimation
# - dispersion estimation
# - negative binomial model fitting
# - Wald test

dds <- DESeq(dds)


################################################
# STEP 3 — Extract Differential Expression Results
################################################
res <- results(dds)

# Order genes by adjusted p-value (FDR)
res <- res[order(res$padj), ]

head(res)

res

#Dispersion Plot (DESeq2)

par(mar=c(4,4,2,1))
plotDispEsts(dds)

#Sample Distance Heatmap
library(pheatmap)
library(RColorBrewer)

vsd <- vst(dds)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))

#PCA Plot (DESeq2)
plotPCA(vsd, intgroup = "dex")




