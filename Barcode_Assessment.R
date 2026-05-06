library(Biostrings)
library(ggplot2)


# do this first to get "genomes" below from unannotated mitochondrial genomes

extract_genomes_hardened <- function(gb_file, verbose = TRUE) {
  # Load Biostrings
  if(!requireNamespace("Biostrings", quietly = TRUE))
    stop("Biostrings package is required")
  library(Biostrings)
  
  # Read the GenBank file
  gb_lines <- readLines(gb_file, warn = FALSE)
  rec_ends <- which(trimws(gb_lines) == "//")
  
  if (length(rec_ends) == 0) stop("No record separators ('//') found in file.")
  
  genomes <- list()
  qc_report <- data.frame(
    genome = character(),
    locus = character(),
    accession = character(),
    declared_length = numeric(),
    extracted_length = numeric(),
    circular = logical(),
    invalid_chars_removed = numeric(),
    truncated_origin = logical(),
    duplicate = logical(),
    stringsAsFactors = FALSE
  )
  
  seen_names <- character()
  rec_start <- 1
  
  for (i in seq_along(rec_ends)) {
    
    rec_lines <- gb_lines[rec_start:(rec_ends[i] - 1)]
    rec_start <- rec_ends[i] + 1
    
    if (length(rec_lines) == 0) next
    
    # ---- LOCUS parsing ----
    locus_line <- grep("^LOCUS", rec_lines, value = TRUE)
    locus_name <- NA
    declared_length <- NA
    circular_flag <- FALSE
    
    if (length(locus_line) > 0) {
      parts <- strsplit(locus_line, "\\s+")[[1]]
      locus_name <- parts[2]
      # Extract declared length if available
      length_match <- regmatches(locus_line, regexpr("[0-9]+ bp", locus_line))
      if (length(length_match) > 0)
        declared_length <- as.numeric(gsub(" bp", "", length_match))
      circular_flag <- grepl("circular", locus_line, ignore.case = TRUE)
    }
    
    # ---- ACCESSION parsing ----
    accession_line <- grep("^ACCESSION", rec_lines, value = TRUE)
    accession <- if (length(accession_line) > 0) strsplit(accession_line, "\\s+")[[1]][2] else NA
    
    genome_name <- if (!is.na(accession)) accession else locus_name
    if (is.na(genome_name)) genome_name <- paste0("Genome_", i)
    
    # ---- Duplicate detection ----
    duplicate_flag <- genome_name %in% seen_names
    if (duplicate_flag) {
      suffix <- 1
      while (paste0(genome_name, "_", suffix) %in% seen_names) suffix <- suffix + 1
      genome_name <- paste0(genome_name, "_", suffix)
    }
    seen_names <- c(seen_names, genome_name)
    
    # ---- ORIGIN extraction ----
    origin_idx <- grep("^ORIGIN", rec_lines)
    truncated_flag <- FALSE
    
    if (length(origin_idx) == 0) {
      if (verbose) message("No ORIGIN block for record ", genome_name, " -> skipping")
      next
    }
    
    dna_lines <- rec_lines[(origin_idx + 1):length(rec_lines)]
    
    # Detect accidental feature bleed into ORIGIN
    if (any(grepl("FEATURES", dna_lines))) truncated_flag <- TRUE
    
    # ---- Safe DNA cleaning ----
    dna_seq_raw <- paste(dna_lines, collapse = "")
    dna_seq_clean <- gsub("[^ACGTNacgtn]", "", dna_seq_raw)  # removes digits, spaces, punctuation
    invalid_removed <- nchar(dna_seq_raw) - nchar(dna_seq_clean)
    dna_seq <- toupper(dna_seq_clean)
    
    if (nchar(dna_seq) == 0) {
      if (verbose) message("Empty sequence after cleaning for ", genome_name, " -> skipping")
      next
    }
    
    # ---- Length check vs declared ----
    extracted_length <- nchar(dna_seq)
    if (!is.na(declared_length) && abs(declared_length - extracted_length) > 5) {
      truncated_flag <- TRUE
      if (verbose)
        message("Length mismatch for ", genome_name,
                ": declared=", declared_length,
                " extracted=", extracted_length)
    }
    
    # Add sequence to list
    genomes[[genome_name]] <- dna_seq
    
    # Update QC
    qc_report <- rbind(qc_report,
                       data.frame(
                         genome = genome_name,
                         locus = locus_name,
                         accession = accession,
                         declared_length = declared_length,
                         extracted_length = extracted_length,
                         circular = circular_flag,
                         invalid_chars_removed = invalid_removed,
                         truncated_origin = truncated_flag,
                         duplicate = duplicate_flag,
                         stringsAsFactors = FALSE
                       ))
  }
  
  # ---- Build final DNAStringSet ----
  # Remove empty sequences if any
  genomes <- genomes[nchar(genomes) > 0]
  
  dna_set <- DNAStringSet(unlist(genomes))
  names(dna_set) <- names(genomes)
  
  if (verbose) {
    message("Extracted ", length(dna_set), " genomes from file.")
    message("Circular genomes detected: ", sum(qc_report$circular, na.rm = TRUE))
    message("Duplicates handled: ", sum(qc_report$duplicate, na.rm = TRUE))
    message("Sequences with invalid chars removed: ", sum(qc_report$invalid_chars_removed > 0))
    message("Truncated ORIGIN blocks: ", sum(qc_report$truncated_origin))
  }
  
  return(list(
    genomes = dna_set,
    qc = qc_report
  ))
}

result <- extract_genomes_hardened("Mussel_mitogenome.gb")

genomes <- result$genomes    # Clean DNAStringSet
qc_table <- result$qc        # QC table for inspection

# Inspect sequences
genomes[1:3]

# Inspect QC
head(qc_table)
subset(qc_table, truncated_origin | invalid_chars_removed > 0 | duplicate)

library(Biostrings)
library(ggplot2)

# Below is code to align all sequences from unannotated sequences, the alleged 16S
# sequences, and our reference sequence that is annotated from a 2024 paper.  It
# aligns all the sequences and checks if our identified region is actually in the expected 
# region from the extracted sequences based on the whole genomes to the reference sequence.
# Follow all steps in order below to complete, only having completed the step above.

# --- Step 1: Inputs ---
# Clean full genomes from extractor
# result <- extract_genomes_hardened("Mussel_mitogenome.gb")
# genomes <- result$genomes
# qc_table <- result$qc

# Extracted 16S sequences from your extraction function
extracted_16S <- readDNAStringSet("16S_all.fasta")

# Trusted reference 16S (PP103562.1)
reference_16S <- readDNAStringSet("Oligamentina_mito_complete.fasta")
ref_seq <- reference_16S[[1]]  # assume single sequence

# --- Step 2: Validation Function ---
library(Biostrings)
library(parallel)
library(ggplot2)

library(Biostrings)
library(parallel)
library(ggplot2)

validate_16S_pipeline_with_locus <- function(genomes, extracted_16S, ref_seq,
                                             locus_names = NULL,
                                             identity_threshold = 97,
                                             truncation_threshold = 50,
                                             window = 20,
                                             verbose = TRUE,
                                             ncores = 2) {
  
  genome_names <- names(genomes)
  
  if (!is.null(locus_names)) {
    if (length(locus_names) != length(genome_names)) {
      stop("locus_names must have the same length as genomes")
    }
    names(genomes) <- locus_names
    genome_names <- locus_names
  }
  
  # Helper for one genome
  validate_one_genome <- function(gname) {
    genome_seq <- genomes[[gname]]
    genome_len <- nchar(as.character(genome_seq))
    
    # Align forward and reverse
    aln_fwd <- pairwiseAlignment(ref_seq, genome_seq, type="local")
    aln_rev <- pairwiseAlignment(ref_seq, reverseComplement(genome_seq), type="local")
    
    pid_fwd <- pid(aln_fwd)
    pid_rev <- pid(aln_rev)
    
    if (pid_fwd >= pid_rev) {
      aln_best <- aln_fwd
      orientation <- "forward"
      subj_seq <- genome_seq
    } else {
      aln_best <- aln_rev
      orientation <- "reverse"
      subj_seq <- reverseComplement(genome_seq)
    }
    
    ref_start <- start(subject(aln_best))
    ref_end <- end(subject(aln_best))
    ref_width <- nchar(as.character(pattern(aln_best)))
    pid_val <- pid(aln_best)
    
    # Extracted 16S comparison
    if (gname %in% names(extracted_16S)) {
      extracted_seq <- extracted_16S[[gname]]
      extracted_len <- nchar(as.character(extracted_seq))
      trunc_bp <- ref_width - extracted_len
      flag <- ifelse(abs(trunc_bp) > truncation_threshold, "TRUNCATED", "OK")
    } else {
      extracted_len <- NA
      trunc_bp <- NA
      flag <- "MISSING"
    }
    
    if (pid_val < identity_threshold) {
      flag <- ifelse(flag == "OK", "LOW_IDENTITY", paste(flag, "LOW_IDENTITY", sep=";"))
    }
    
    # --- Alignment window for inspection ---
    win_start <- max(1, ref_start - window)
    win_end <- min(genome_len, ref_end + window)
    alignment_window <- as.character(subseq(subj_seq, start = win_start, end = win_end))
    
    return(data.frame(
      locus = gname,
      genome_length = genome_len,
      extracted_length = extracted_len,
      ref_aln_start = ref_start,
      ref_aln_end = ref_end,
      ref_aln_width = ref_width,
      percent_identity = pid_val,
      truncation_bp = trunc_bp,
      orientation = orientation,
      flag = flag,
      alignment_window = alignment_window,
      stringsAsFactors = FALSE
    ))
  }
  
  # --- Parallel execution ---
  results_list <- mclapply(genome_names, validate_one_genome,
                           mc.cores = ncores)
  
  validation_results <- do.call(rbind, results_list)
  
  if (verbose) {
    message("Validation complete for ", nrow(validation_results), " genomes.")
    message("Flags summary:")
    print(table(validation_results$flag))
  }
  
  return(validation_results)
}

# --- Step 3: Run the pipeline ---
# Run on 2 cores
# Run pipeline
validation_results_locus <- validate_16S_pipeline_with_locus(
  genomes = genomes,
  extracted_16S = extracted_16S,
  ref_seq = ref_seq,
  locus_names = qc_table$locus,   # vector of locus names from your GB file
  ncores = 2,
  verbose = TRUE
)


# --- Step 4: Diagnostic plots ---
# Histogram of 16S start positions by locus
ggplot(validation_results_locus, aes(x = ref_aln_start)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of 16S Start Positions", x = "Start Position", y = "Count")

# Truncation bar plot
ggplot(validation_results_locus, aes(x = locus, y = truncation_bp, fill = flag)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "16S Truncation Across Loci", x = "Locus", y = "Truncation bp")

# Percent identity scatter plot
ggplot(validation_results_locus, aes(x = locus, y = percent_identity, color = flag)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "16S Percent Identity vs Locus", x = "Locus", y = "Percent Identity")

# Quick way to inspect alignment window for a specific locus
validation_results_locus$alignment_window[5]  # first locus


# Assuming validation_results_locus from
# validate_16S_pipeline_with_locus()

# Distinguish MISSING loci clearly
validation_results_locus$flag <- factor(validation_results_locus$flag,
                                        levels = c("OK", "TRUNCATED", "LOW_IDENTITY", "MISSING"))

# Percent Identity Scatter Plot with MISSING highlighted
ggplot(validation_results_locus, aes(x = locus, y = percent_identity, color = flag, shape = flag)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("OK" = "green", 
                                "TRUNCATED" = "orange", 
                                "LOW_IDENTITY" = "red", 
                                "MISSING" = "black")) +
  scale_shape_manual(values = c("OK" = 16, 
                                "TRUNCATED" = 17, 
                                "LOW_IDENTITY" = 18, 
                                "MISSING" = 4)) + # cross for MISSING
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "16S Percent Identity by Locus",
       x = "Locus",
       y = "Percent Identity",
       color = "QC Flag",
       shape = "QC Flag")

##### MULTIPLE ALIGNMENT OF ALL SEQUENCES TO SEE WHAT WE HAVE #####

library(Biostrings)
library(DECIPHER)
library(ggplot2)

# ---- Step 1: Combine sequences ----
# Include reference + extracted sequences
# Wrap the reference sequence in a DNAStringSet
ref_set <- DNAStringSet(ref_seq)

# Set the name afterwards
names(ref_set) <- "Reference_16S"

# Combine with extracted 16S sequences
all_seqs <- c(ref_set, extracted_16S)
names(all_seqs)[1] <- "Reference_16S"

# Make a DNAStringSet
all_seqs_set <- DNAStringSet(all_seqs)

# ---- Step 2: Perform multiple sequence alignment using DECIPHER ----
aligned_seqs <- AlignSeqs(all_seqs_set, iterations = 2, refinements = 1, verbose = FALSE)

# ---- Step 3: Compute consensus sequence ----
consensus_seq <- ConsensusSequence(aligned_seqs, threshold = 0.5, ambiguity = FALSE)

# Convert to DNAStringSet for export
consensus_fasta <- DNAStringSet(consensus_seq)
names(consensus_seq) <- "Consensus_16S"

# ---- Step 4: Export consensus FASTA for BLAST ----
writeXStringSet(consensus_fasta, "consensus_16S_DECIPHER.fasta")

