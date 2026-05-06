library(Biostrings)

# Reference dataset input
ref_set <- DNAStringSet(ref_seq)
names(ref_set) <- "NCBI_Reference_16S"

# Extracted sequences (example)
# Ensure they already have meaningful names
# e.g. names(extracted_16S) <- qc_table$locus

length(ref_set)        # should be 1
length(extracted_16S)  # should be your 130 sequences

validate_extracted_16S <- function(ref_set,
                                   extracted_16S,
                                   output_prefix="16S_validation") {
  
  library(Biostrings)
  library(DECIPHER)
  
  # -------------------------
  # Validate inputs
  # -------------------------
  if (!inherits(ref_set, "DNAStringSet"))
    stop("ref_set must be a DNAStringSet")
  
  if (!inherits(extracted_16S, "DNAStringSet"))
    stop("extracted_16S must be a DNAStringSet")
  
  if (length(ref_set) != 1)
    stop("ref_set must contain exactly one sequence")
  
  if (is.null(names(ref_set)))
    stop("ref_set must be named")
  
  if (is.null(names(extracted_16S)))
    stop("extracted_16S must have names")
  
  # -------------------------
  # Combine sequences
  # -------------------------
  all_seqs <- c(ref_set, extracted_16S)
  
  # -------------------------
  # Multiple alignment
  # -------------------------
  aligned <- AlignSeqs(all_seqs,
                       iterations=2,
                       refinements=1,
                       verbose=FALSE)
  
  aln_matrix <- as.matrix(aligned)
  ref_aln <- aln_matrix[1, ]
  
  # -------------------------
  # Consensus sequence
  # -------------------------
  consensus <- ConsensusSequence(aligned,
                                 threshold=0.5,
                                 ambiguity=FALSE)
  
  # Safe coercion (version-proof)
  if (inherits(consensus, "DNAString")) {
    consensus_set <- DNAStringSet(consensus)
  } else {
    consensus_set <- consensus
  }
  
  names(consensus_set) <- "Consensus_16S"
  
  writeXStringSet(consensus_set,
                  paste0(output_prefix, "_consensus.fasta"))
  
  # -------------------------
  # Plot 1: mismatch heatmap
  # -------------------------
  mismatch_matrix <- sweep(aln_matrix[-1, ], 2, ref_aln, FUN="!=")
  
  pdf(paste0(output_prefix, "_mismatch_heatmap.pdf"),
      width=10, height=6)
  
  image(t(mismatch_matrix),
        col=c("white","red"),
        axes=FALSE,
        main="Mismatch vs Reference (Red = mismatch)")
  
  dev.off()
  
  # -------------------------
  # Plot 2: per-position conservation
  # -------------------------
  conservation <- apply(aln_matrix, 2, function(col) {
    col <- col[col != "-"]
    if (length(col) == 0) return(NA)
    max(prop.table(table(col)))
  })
  
  pdf(paste0(output_prefix, "_conservation_plot.pdf"),
      width=10, height=4)
  
  plot(conservation, type="l",
       ylab="Max Base Frequency",
       xlab="Alignment Position",
       main="Per-Position Conservation")
  
  abline(h=0.5, col="red", lty=2)
  
  dev.off()
  
  # -------------------------
  # Plot 3: Reference coverage
  # -------------------------
  ref_positions <- cumsum(ref_aln != "-")
  
  coverage <- apply(aln_matrix[-1, ], 2, function(col) {
    sum(col != "-")
  })
  
  pdf(paste0(output_prefix, "_reference_coverage.pdf"),
      width=10, height=4)
  
  plot(ref_positions, coverage,
       type="l",
       xlab="Reference 16S Coordinate",
       ylab="Number of Sequences Covering",
       main="Coverage Across Reference 16S")
  
  dev.off()
  
  # -------------------------
  # Return objects
  # -------------------------
  return(list(
    aligned = aligned,
    consensus = consensus_set
  ))
}