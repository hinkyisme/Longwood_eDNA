# This code is for installation of dada2 utilizing Bioconductor, optimal in R v. 4.4.x but currently in 4.3 using Bioconductor 3.18 
# Additionally this code also serves as an example walkthrough of dada2 using an example dataset 

# Nick Duellman 
# Sources: 
# Bioconductor install https://www.bioconductor.org/install/  
# DADA 2 Install http://benjjneb.github.io/dada2/tutorial.html

#Bioconductor Install

# Core packages install


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

# View available packages

BiocManager::available()


# dada2 install 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# load dada2 
BiocManager::install("dada2")
library(dada2); packageVersion("dada2")

# inspect sample data with file path (example paired end fastq files)
path <- "/Users/jamesonhinkle/Desktop/folders/Longwood/Research/eDNA/R/Longwood_eDNA/dada2 supplemental example data - MiSeq_SOP" #change depending on file location
list.files(path)

# sort and match files for forward and reverse files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE)) #dont use Rs for single read
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# String manipulations may need modification depending on file format

# inspecting read quality (forward and reverse read qualities)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#reads must overlap after truncation in order to merge them! 

# Filter and Trim | Create subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, truncLen=c(240),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
#11/11/24, why truncLen of 240?


head(out)

#Note: Speed up downstream computation, tighten maxEE. Relax maxEE if too few reads are passing through, and reduce truncLen to remove
#low quality reads. 

# Error rates 

errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

plotErrors(errF, nominalQ=TRUE) #Estimated error rates (black line) are good fit with observed (points), and error rates drop with increased quality of samples.


# Infer sample sequences - ID and quantify the exact ASV present in a sample

dadaFs <- dada(filtFs, err=errF, multithread = FALSE) 
dadaRs <- dada(filtRs, err=errF, multithread = FALSE) #dont use reverse reads for single read

dadaFs[[1]]
help("dada-class") #explore dada-class diagnostics


# pooling information - dada2 package offers two types dada(.., pool = TRUE) for standard pooling


# Merge forward and reverse for full denoised sequences
# **Dont merge files when working with single read data**
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#Inspect the merger data frame from sample 1 
head(mergers[[1]])

# Construct sequence table

#seqtab <- makeSequenceTable(mergers)
seqtab <- makeSequenceTable(dadaFs)
#dim(seqtab)
dim(seqtab)

# Inspect dist of length sequences

table(nchar(getSequences(seqtab)))

# Remove chimeras


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Check reads that made it through pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)


# Assign taxonomy


