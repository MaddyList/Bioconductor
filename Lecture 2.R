# BIOCONDUCTOR – INSTALLATION 

# Install BiocManager (Bioconductor installer)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Check Bioconductor version
BiocManager::version()


# INSTALL CORE BIOCONDUCTOR PACKAGES
BiocManager::install(c(
  "Biostrings",         # Work with DNA/RNA/Protein sequences
  "BSgenome",           # Whole genomes
  "BSgenome.Scerevisiae.UCSC.sacCer3",  # Yeast genome
  "GenomicRanges",      # Genomic coordinates
  "IRanges",            # Interval objects
  "AnnotationHub"       # For downloading gene annotation
))


BiocManager::install(c("BSgenome.Scerevisiae.UCSC.sacCer3"))
BiocManager::install(c("AnnotationHub"))



# LOAD LIBRARIES
library(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(GenomicRanges)
library(IRanges)
library(AnnotationHub)


# PART 1 — BIOSTRINGS: Learning DNA/RNA/Protein in R
dna <- DNAString("ACGTACGTNNNN")
dna

reverse(dna)
complement(dna)
reverseComplement(dna)


# 1.2 DNAStringSet: multiple sequences 
dna_set <- DNAStringSet(c("ACGT", "GTCA", "TTTTTT", "GGGG"))
dna_set


# 1.3 Basic operations
reverse(dna)                 # Reverse order
complement(dna)              # A<->T, C<->G
reverseComplement(dna)       # Very important in biology


# 1.4 Translation: DNA → protein
translate(DNAString("ATGTTCGGA"))


# PART 2 — BSGENOME: Working with a FULL Yeast Genome
yeast <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast


# 2.2 Explore genome structure
seqnames(yeast)       # Chromosome names
seqlengths(yeast)     # Chromosome lengths


# 2.3 Extract chromosome sequence
chrI <- yeast$chrI
chrI
subseq(chrI, 1, 50)   # first 50 bases


# PART 3 — Pattern Matching (Finding DNA Motifs)
## “We can search for any DNA pattern in the genome.”

pattern <- DNAString("ACGTACGT")

# Count matches
countPattern(pattern, chrI)

# Find locations
matchPattern(pattern, chrI)


# PART 4 — IRanges & GenomicRanges (Coordinates on a Genome)
## “IRanges = start/end positions.
## GRanges = positions with chromosomes.”


# 4.1 IRanges basics
IRanges(start = c(10, 50, 100), width = 20)


# 4.2 GRanges basics
GRanges(
  seqnames = c("chrI", "chrI"),
  ranges = IRanges(start = c(100, 500), width = 50),
  strand = c("+", "-")
)


# PART 5 — Example: Promoter GC Content
# Promoters are upstream regions of genes; their GC content can affect gene regulation.

# 5.1 Load yeast gene annotation using AnnotationHub
ah <- AnnotationHub()
query(ah, "sacCer3")   # search for yeast genome annotations


# 5.2 Load the yeast gene database (TxDb)
# Replace AH52272 with the TxDb ID that appears in your query() results
txdb <- ah[["AH52272"]]     
genes <- genes(txdb)


# 5.3 Extract promoter regions (2000 bp upstream of each gene)
promoters <- promoters(genes, upstream = 2000, downstream = 0)
promoters <- trim(promoters)

promoters


# 5.4 Compute GC content
promoter_seqs <- getSeq(BSgenome.Scerevisiae.UCSC.sacCer3, promoters)
gc_content <- letterFrequency(promoter_seqs, "GC", as.prob = TRUE)

# View summary of GC content across all promoters
summary(gc_content)


# PART 6 — Coverage & Rle (Read Depth)
##   “Rle compresses repeated values, useful for coverage.”
reads <- GRanges(seqnames = "chrI",
  ranges = IRanges(start = c(1,5,10), width = 5)
)

cov <- coverage(reads)
cov


reads <- GRanges(seqnames = "chrI"), 
ranges = IRanges(start = c(1,5,10), width =5))

cov <- coverage(reads)





# End
sessionInfo()

