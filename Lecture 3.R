#SECTION 1 — Installing and Loading Packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase", "ALL", "hgu95av2.db"))


library(Biobase)
library(ALL)
library(hgu95av2.db)


#SECTION 2 — Loading Example Data
data(ALL)
ALL

#SECTION 3 — Accessing Data Inside ExpressionSet
exprs(ALL)[1:4, 1:4]
sampleNames(ALL)
featureNames(ALL)
pData(ALL)
ALL$sex
ALL$age
pData(ALL)$sex


exprs()                      # experiment measurements  
pData()                      # sample info  
featureNames()               # gene info  
sampleNames()                # sample IDs


#SECTION 4 — Subsetting ExpressionSet
ALL[, 1:5]
ALL[1:10, ]
ALL[1:10, 1:5]
ALL[, c(3, 2, 1)]


#SECTION 5 — Annotation Using hgu95av2.db
ids <- featureNames(ALL)[1:5]
as.list(hgu95av2ENTREZID[ids])


#SECTION 6 — SummarizedExperiment Basics
BiocManager::install("airway")
library(airway)
library(GenomicRanges)

data(airway)
airway


#SECTION 7 — Accessing SummarizedExperiment Data
colData(airway)
colnames(airway)
rownames(airway)
assayNames(airway)
assay(airway, "counts")
rowRanges(airway)[[1]]
start(airway)[1:10]


#SECTION 8 — Genomic Subsetting
gr <- GRanges(seqnames = "1", ranges = IRanges(start = 1, end = 10^7))
subsetByOverlaps(airway, gr)

subsetByOverlaps()


#SECTION 9 — GEOquery Basics
BiocManager::install("GEOquery")
library(GEOquery)

getGEO("GSE11675")
getGEOSuppFiles("GSE11675")

getwd()
gset <- getGEO(filename = "GSE11675_series_matrix.txt.gz")
gset


#SECTION 10 — biomaRt Basics
BiocManager::install("biomaRt")
library(biomaRt)
listMarts()

mart <- useMart("ensembl")

mart <- useMart

listDatasets(mart)
listDatasets(mart)

ensembl <- useDataset("hsapiens_gene_ensembl", mart)

ensembl <- useDatase ("hsapiens_gene_ensemb", mart)

ensembl

listAttributes(ensembl)
listFilters(ensembl)

listAttributes()
listFilters()



#SECTION 11 — S4 Classes and Methods
class(ALL)
isS4(ALL)
isS4
getClass("ExpressionSet")
showMethods("as.data.frame")
getMethod("as.data.frame", "DataFrame")
validObject(ALL)

validObject()
