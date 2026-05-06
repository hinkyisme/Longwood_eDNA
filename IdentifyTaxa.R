# --- 1. LOAD THE DATA ---
# Read the 130 genomes from the single GenBank file
# (Note: readDNAStringSet works on .gb files by grabbing the 'ORIGIN' sequence)
all_genomes <- readDNAStringSet("mussel_genomes.fasta")

# 2 Load the 9 .seq test files
setwd("/Users/jamesonhinkle/Desktop/folders/Longwood/Research/eDNA/R/Longwood_eDNA/Sanger_BarcodeAlign_test_GoodRead/")
test_files <- list.files(pattern = "\\.seq$", full.names = TRUE)
dna_list <- lapply(test_files, function(f) {
  tryCatch({
    readDNAStringSet(f)
  }, error = function(e) {
    stop(paste("Error reading file:", f, "\nDetails:", e$message))
  })
})
setwd("/Users/jamesonhinkle/Desktop/folders/Longwood/Research/eDNA/R/Longwood_eDNA/")
test_seqs <- do.call(c, dna_list)
clean_file_names <- basename(as.character(test_files))
clean_file_names <- gsub("\\.seq$", "", clean_file_names, ignore.case = TRUE)
names(test_seqs) <- clean_file_names


# --- 3. MANUAL 16S EXTRACTION ---

f_primer <- "GGGTCACCAACTCCGCTAAC" 
r_primer <- "GGGCAGCTAAGGCTGGAAAA"

# NOTE: these are the primers from A22 in "Selected 16s primers".
# They also only work on whole genomes with gaps still contained.
# We likely need a better primer.

# Find the Forward Primer positions
# (Allowing 2 mismatches manually)
f_hits <- vmatchPattern(f_primer, all_genomes, max.mismatch = 2)

# Find the Reverse Primer positions 
# (We look for the reverse complement of the R primer)
r_hits <- vmatchPattern(reverseComplement(DNAString(r_primer)), all_genomes, max.mismatch = 2)

# Extract only the sequences where BOTH primers were found
valid_idx <- which(elementNROWS(f_hits) > 0 & elementNROWS(r_hits) > 0)

# Create the training set by cutting between the primers
training_list <- lapply(valid_idx, function(i) {
  start_pos <- end(f_hits[[i]][1]) + 1  # Start right after F primer
  end_pos <- start(r_hits[[i]][1]) - 1   # End right before R primer
  
  # Only extract if it's a reasonable length (prevents circular genome errors)
  if ((end_pos - start_pos) > 400 & (end_pos - start_pos) < 2500) {
    return(subseq(all_genomes[[i]], start = start_pos, end = end_pos))
  } else {
    return(NULL)
  }
})

# Combine into a clean DNAStringSet
# A. Combine the list into a DNAStringSet (plural)
# Using 'DNAStringSet' around the 'unlist' ensures it's the right class
training_set <- DNAStringSet(unlist(training_list))

# B. Extract the corresponding names from your original 130 genomes
# We filter the names to match only the sequences that survived the extraction
original_names <- names(all_genomes)[valid_idx][!sapply(training_list, is.null)]

# C. Assign the names
# This should now work without the 'S4 object' error
names(training_set) <- original_names

# D. Final Sanity Check
print(training_set)

# Check success
cat("Extracted", length(training_set), "valid 16S barcodes.\n")

# --- 5. TRAIN THE TAXONOMIC CLASSIFIER (Memory Optimized) ---

# 1. Clean and format labels as before
species_names <- gsub("^(\\w+\\s\\w+).*", "\\1", names(training_set))
train_labels <- paste0("Root;Mitochondria;", species_names, ";")

# 2. Build the model with memory-saving parameters
# maxKmerSize = 8 (default is 8, reducing it significantly drops memory usage)

trained_model <- LearnTaxa(training_set, 
                           train_labels, 
                           K = 8)

# 3. Clean up R's memory immediately after training
gc()

# --- 6. ASSIGN TAXONOMY TO TEST FILES ---
# strand = "both" is critical if test files might be in reverse orientation
assignments <- IdTaxa(test_seqs, 
                      trained_model, 
                      strand = "both", 
                      threshold = 50)

# --- 7. FORMAT AND VIEW RESULTS ---
results <- data.frame(
  Test_File = names(test_seqs),
  Taxonomy = sapply(assignments, function(x) paste(x$taxon, collapse = " > ")),
  Confidence = sapply(assignments, function(x) min(x$confidence))
)

print(results)

# --- 8. COMBINE AND ALIGN ---
# Join your references and test sequences into one set
all_for_alignment <- c(training_set, test_seqs)

# Perform the multiple sequence alignment (MSA)
# This uses the DECIPHER 'AlignSeqs' function
alignment <- AlignSeqs(all_for_alignment)

# --- 9. VIEW IN BROWSER ---
# This opens an interactive window
# highlight = 0 shows all bases
# highlight = 1 shows only the differences compared to the first sequence
BrowseSeqs(alignment, highlight = 1)

# --- 10. EXTRACT RESULTS ---
# This function pulls the most specific (last) name that passed your 50% threshold
final_verdict <- sapply(assignments, function(x) {
  # The 'taxon' vector contains: "Root", "Mitochondria", "Genus species"
  # We take the last element in that vector
  tail(x$taxon, n = 1)
})

# --- 11. CREATE A SUMMARY TABLE ---
# Combine the filename, the final assigned taxon, and the confidence
summary_table <- data.frame(
  Sample_Name = names(test_seqs),
  Assigned_Taxon = final_verdict,
  Confidence_Score = sapply(assignments, function(x) {
    # We take the confidence of that final assigned level
    tail(x$confidence, n = 1)
  })
)

# Sort by confidence so you see the most certain hits first
summary_table <- summary_table[order(summary_table$Confidence_Score, decreasing = TRUE), ]

print(summary_table)
write.csv(summary_table, "Barcode_test_summary.csv")
