#QA/QC test script 04/07/2026

library(Biostrings)
library(DECIPHER)

# Load your sequences
dna <- readDNAStringSet("mussel_genomes.fasta")

# 1. Remove gaps or unusual characters if they exist
dna <- RemoveGaps(dna)

# 2. Orient all sequences in the same direction
# This is crucial for mitochondria which might be sequenced in different directions
dna <- OrientNucleotides(dna)

# Perform the alignment
aligned_dna <- AlignSeqs(dna)

# Browse the alignment in your browser to check for quality
BrowseSeqs(aligned_dna)

BiocManager::install("ggmsa")
library(ggmsa)

# Plot a specific range (e.g., the first 100 bases)
ggmsa(aligned_dna, start = 1, end = 100, char_width = 0.5, seq_name = TRUE) + 
  geom_seqlogo() + 
  geom_msaBar()

library(ggplot2)

# This will plot the conservation/identity across the whole genome
# Note: For 130 full genomes, this plot may be very dense.
ggmsa(aligned_dna, font = NULL, color = "Clustal") + 
  geom_msaBar() + 
  theme_minimal() +
  labs(title = "Global Conservation Profile of 130 Mitochondrial Genomes",
       x = "Position (bp)", 
       y = "Genomes")

# above may be a little much for my machine (at least for a quick look)
# see below for quick proportion of sequences aligned plot

alignment_matrix <- as.matrix(aligned_dna)
coverage <- apply(alignment_matrix, 1, function(x) sum(x != "-") / length(x))
coverage_df <- data.frame(
  Sequence = names(aligned_dna),
  Coverage = round(coverage * 100, 2)
)
head(coverage_df)
summary(coverage_df$Coverage) # percent coverage alignment

col_occupancy <- apply(alignment_matrix, 2, function(x) mean(x != "-"))
mean(col_occupancy) * 100

plot(col_occupancy, type="l", main="Column Occupancy per Position",
     ylab="Proportion of Sequences Present", xlab="Alignment Position") # plot of that coverage


# Export the aligned sequences
writeXStringSet(aligned_dna, file = "musselgenomes_gapsrm.fasta")
