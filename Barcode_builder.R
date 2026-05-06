# ----------------
# Script: Extract 16S rRNA from a single GenBank flat file
# ----------------
library(ape)
library(stringr)
library(dplyr)
library(purrr)

# To check for missing/necessary metadata, use the following:
readLines("ecomp_sequence.gb")

# --- User params ---
gb_file <- "ecomp_sequence.gb"      # <-- replace with your file
out_prefix <- "custom_16S_ref"
min_len <- 200
# -------------------

# this function below pulls all 130 sequences from the genbank flat file and gives you a 
# fasta sequence for each that we can turn into a barcode

extract_16S_multi_to_fasta <- function(gb_file,
                                       out_fasta = "16S_all.fasta",
                                       wrap = 80,
                                       verbose = TRUE) {
  # helper: reverse complement
  rc <- function(seq) {
    seq <- toupper(seq)
    seq_chars <- rev(strsplit(seq, "")[[1]])
    seq_rc <- paste(seq_chars, collapse = "")
    # complement A<->T, C<->G; leave N and others unchanged
    chartr("ACGT", "TGCA", seq_rc)
  }
  
  # read file
  gb_lines <- readLines(gb_file, warn = FALSE)
  # find record end lines (exactly "//")
  rec_ends <- which(gb_lines == "//")
  if (length(rec_ends) == 0) stop("No record separators ('//') found in file.")
  
  seq_list <- list()
  rec_start <- 1
  rec_idx <- 0
  
  for (end_idx in rec_ends) {
    rec_idx <- rec_idx + 1
    rec_lines <- gb_lines[rec_start:(end_idx - 1)]
    rec_start <- end_idx + 1
    
    if (length(rec_lines) == 0) next
    
    # LOCUS & ACCESSION for naming
    locus_line <- grep("^LOCUS", rec_lines, value = TRUE)
    accession_line <- grep("^ACCESSION", rec_lines, value = TRUE)
    locus_name <- NA
    if (length(locus_line) > 0) {
      locus_name <- strsplit(locus_line, "\\s+")[[1]][2]
    } else if (length(accession_line) > 0) {
      locus_name <- strsplit(accession_line, "\\s+")[[1]][2]
    } else {
      locus_name <- paste0("Genome_", rec_idx)
    }
    
    # ORIGIN: build raw sequence string (remove numbers/spaces)
    origin_idx <- grep("^ORIGIN", rec_lines)
    if (length(origin_idx) == 0) {
      if (verbose) message("Record ", rec_idx, " has no ORIGIN -> skipped")
      next
    }
    dna_lines <- rec_lines[(origin_idx + 1):length(rec_lines)]
    dna_seq <- gsub("[0-9\\s]", "", paste(dna_lines, collapse = ""))
    dna_seq <- toupper(dna_seq)
    if (nchar(dna_seq) == 0) next
    
    # FEATURES block
    feat_idx <- grep("^FEATURES", rec_lines)
    if (length(feat_idx) == 0) {
      if (verbose) message("Record ", rec_idx, " has no FEATURES -> skipped")
      next
    }
    feat_lines <- rec_lines[(feat_idx + 1):(origin_idx - 1)]
    # identify feature entry starts: lines with non-space at column 6 (i.e., starts with 5 spaces then non-space)
    entry_start <- which(grepl("^\\s{5}\\S", feat_lines, perl = TRUE))
    if (length(entry_start) == 0) next
    
    # build full feature entries (each entry = multiple wrapped lines)
    entry_ranges <- c(entry_start, length(feat_lines) + 1)
    feature_entries <- vector("list", length(entry_start))
    for (k in seq_along(entry_start)) {
      s <- entry_start[k]
      e <- entry_ranges[k + 1] - 1
      entry_text <- paste(feat_lines[s:e], collapse = " ")
      feature_entries[[k]] <- entry_text
    }
    
    # search entries for 16S-like features (check both location and qualifiers)
    hits <- which(sapply(feature_entries, function(x) {
      grepl("16S|16 s|16_S|16S ribosomal|16S rRNA|16S_RNA|16s|16S ribosomal RNA", x, ignore.case = TRUE)
    }))
    if (length(hits) == 0) {
      # also look for generic 'rRNA' with product containing 16S in qualifiers
      hits <- which(sapply(feature_entries, function(x) {
        grepl("rRNA", x, ignore.case = TRUE) && grepl("16S|16s", x, ignore.case = TRUE)
      }))
    }
    if (length(hits) == 0) {
      if (verbose) message("Record ", rec_idx, " no candidate 16S features found -> skipped")
      next
    }
    
    # For each hit, extract coordinates and sequence
    seq_count_in_rec <- 0
    for (h in hits) {
      entry <- feature_entries[[h]]
      # find coordinate patterns like 123..456
      coords <- regmatches(entry, gregexpr("[0-9]+\\.{2}[0-9]+", entry, perl = TRUE))[[1]]
      if (length(coords) == 0) next
      
      # check complement
      is_complement <- grepl("complement\\s*\\(", entry, ignore.case = TRUE)
      
      # extract parts (support join/multiple ranges)
      parts <- sapply(coords, function(coord) {
        se <- as.numeric(unlist(strsplit(coord, "\\.\\.")))
        if (length(se) != 2) return("")
        if (se[1] < 1) se[1] <- 1
        if (se[2] > nchar(dna_seq)) se[2] <- nchar(dna_seq)
        substr(dna_seq, se[1], se[2])
      }, USE.NAMES = FALSE)
      seq_string <- paste(parts, collapse = "")
      seq_string <- gsub("\\s+", "", seq_string)
      if (nchar(seq_string) == 0) next
      
      if (is_complement) seq_string <- rc(seq_string)
      
      # name: if multiple 16S in the same record, index them
      seq_count_in_rec <- seq_count_in_rec + 1
      seq_name <- if (seq_count_in_rec == 1) locus_name else paste0(locus_name, "_16Scopy", seq_count_in_rec)
      # ensure unique names across file
      if (seq_name %in% names(seq_list)) {
        suffix <- 1
        while (paste0(seq_name, "_", suffix) %in% names(seq_list)) suffix <- suffix + 1
        seq_name <- paste0(seq_name, "_", suffix)
      }
      seq_list[[seq_name]] <- seq_string
    } # end hits loop
  } # end records loop
  
  if (length(seq_list) == 0) stop("No 16S sequences extracted from file.")
  
  # write FASTA (variable-length sequences allowed)
  con <- file(out_fasta, "w")
  for (nm in names(seq_list)) {
    writeLines(paste0(">", nm), con = con)
    wrapped <- strwrap(seq_list[[nm]], width = wrap)
    writeLines(wrapped, con = con)
  }
  close(con)
  
  if (verbose) message("Wrote ", length(seq_list), " sequences to ", out_fasta)
  return(invisible(seq_list))
}

# extract and write all 16S to FASTA
seqs <- extract_16S_multi_to_fasta("Mussel_mitogenome.gb", out_fasta = "16S_all.fasta", wrap = 80)

# check how many were extracted
length(seqs)   # should be 130 if everything found
names(seqs)[1:5]  # inspect first few names


# QA/QC of alignment

alignment_matrix <- as.matrix(aln)
coverage <- apply(alignment_matrix, 1, function(x) sum(x != "-") / length(x))
coverage_df <- data.frame(
  Sequence = names(aln),
  Coverage = round(coverage * 100, 2)
)
head(coverage_df)
summary(coverage_df$Coverage) # percent coverage alignment

col_occupancy <- apply(alignment_matrix, 2, function(x) mean(x != "-"))
mean(col_occupancy) * 100

plot(col_occupancy, type="l", main="Column Occupancy per Position",
     ylab="Proportion of Sequences Present", xlab="Alignment Position") # plot of that coverage
