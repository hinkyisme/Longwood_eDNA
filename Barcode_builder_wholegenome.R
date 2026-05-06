#Barcode_builder V2
#J.Hinkle
#Nov 11, 2025

#This script/function instead of grepling the 16S region specifically just takes a relatively poorly annotated 
#genbank flat file and grabs the necessary annotations (e.g. Locus, Accession, etc.) and puts it into
#one nice fasta file for alighnment.


extract_genomes_to_fasta <- function(gb_file, output_fasta = "genomes.fasta") {
  # Read file
  gb_text <- readLines(gb_file, warn = FALSE)
  
  # Split records by GenBank separator
  records <- strsplit(paste(gb_text, collapse = "\n"), "//\\s*\n", perl = TRUE)[[1]]
  records <- records[nzchar(records)]
  
  # Split records into lines
  rec_lines <- strsplit(records, "\n")
  
  # Extract names
  extract_name <- function(rec) {
    locus_line <- grep("^LOCUS", rec, value = TRUE)
    acc_line   <- grep("^ACCESSION", rec, value = TRUE)
    if (length(locus_line) > 0) {
      strsplit(locus_line[1], "\\s+")[[1]][2]
    } else if (length(acc_line) > 0) {
      strsplit(acc_line[1], "\\s+")[[1]][2]
    } else {
      paste0("Genome_", sample(1e6, 1))
    }
  }
  
  # Extract sequences
  extract_seq <- function(rec) {
    start <- grep("^[[:space:]]*[Oo][Rr][Ii][Gg][Ii][Nn]", rec)
    if (length(start) == 0) return(NA_character_)
    seq_lines <- rec[(start + 1):length(rec)]
    seq <- gsub("[^ATGCatgc]", "", paste(seq_lines, collapse = ""))
    if (nchar(seq) == 0) return(NA_character_)
    toupper(seq)
  }
  
  # Apply both vectorized
  seqs  <- vapply(rec_lines, extract_seq, character(1))
  names <- vapply(rec_lines, extract_name, character(1))
  
  # Filter out invalids
  valid <- nzchar(seqs)
  seqs  <- seqs[valid]
  names <- names[valid]
  
  # Build FASTA lines
  fasta_lines <- unlist(
    mapply(function(nm, sq) c(paste0(">", nm), strwrap(sq, width = 80)),
           names, seqs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  )
  
  # Write file
  writeLines(fasta_lines, output_fasta)
  message("✅ Extracted ", length(seqs), " genomes to ", output_fasta)
  
  invisible(output_fasta)
}

extract_genomes_to_fasta("Mussel_mitogenome_full.gb", "mussel_genomes.fasta")

##### Plotting below

library(DECIPHER)

# Read in your FASTA file
seqs <- readDNAStringSet("mussel_genomes.fasta")

# Align the sequences
alignment <- AlignSeqs(seqs, processors = parallel::detectCores())

# Save the alignment (optional)
WriteXStringSet(alignment, "mussel_genomes_aligned.fasta")

#Native DECIPHER plotting/visualization in html
BrowseSeqs(alignment)

  #If we want the E.comp genome as a reference
file.copy("reference.fasta", "mussel_genomes.fasta", overwrite = FALSE)
write("\n", "mussel_genomes.fasta", append = TRUE)
file.append("mussel_genomes.fasta", "mussel_others.fasta")

#QA/QC of alignment

alignment_matrix <- as.matrix(alignment)
coverage <- apply(alignment_matrix, 1, function(x) sum(x != "-") / length(x))
coverage_df <- data.frame(
  Sequence = names(alignment),
  Coverage = round(coverage * 100, 2)
)
head(coverage_df)
summary(coverage_df$Coverage) # percent coverage alignment

col_occupancy <- apply(alignment_matrix, 2, function(x) mean(x != "-"))
mean(col_occupancy) * 100

plot(col_occupancy, type="l", main="Column Occupancy per Position",
     ylab="Proportion of Sequences Present", xlab="Alignment Position") # plot of that coverage






