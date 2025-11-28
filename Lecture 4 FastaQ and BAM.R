
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


bioc_packages <- c("ShortRead", "Rsamtools")

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

library("ShortRead")
library("Rsamtools")
library("leukemiasEset")
library("readr")
library("data.table")
library("dplyr")
library("ggplot2")



fastqDir <- system.file("extdata", "E-MTAB-1147", package = "ShortRead")

fastqPath <- list.files(fastqDir, pattern = ".fastq.gz$", full = TRUE)[1]
reads <- readFastq(fastqPath)
reads


sread(reads)[1:2]     #extracts the DNA sequences stored in the FASTQ object

quality(reads)[1:2]   #extracts the quality part of the FASTQ read

id(reads)[1:2]

as(quality(reads), "matrix")[1:2,1:10]   #Converts ASCII characters to numeric




getwd()

# Define the file name(s)
fastq_file_R1 <- "ERR127302_1.fastq.gz" 

# Read the single-end file (or Read 1 of a paired-end set)
reads_R1 <- readFastq(fastq_file_R1)

# Check the object type
class(reads_R1)

# Review a summary of the reads
reads_R1


sread(reads_R1)[1:2]











##Reading BAM Files with Rsamtools (Bioconductor)

bamPath <- system.file("extdata", "ex1.bam", package = "Rsamtools")
bamFile <- BamFile(bamPath)
bamFile



seqinfo(bamFile)

seqinfo(bamFile)

aln <- scanBam(bamFile)

scanBam()



aln <- aln[[1]]
names(aln)





lapply(aln, function(x) x[1])








yieldSize(bamFile) <- 1000
open(bamFile)



chunk1 <- scanBam(bamFile)
chunk2 <- scanBam(bamFile)

chunk1
chunk2

close(bamFile)


countBam(bamFile)
countBam()




bam_file_name <- "06.4_raw_merged_f4_sort_rmdup_q30_IR_30bp.bam"

bam_file_name
countBam(bam_file_name)










