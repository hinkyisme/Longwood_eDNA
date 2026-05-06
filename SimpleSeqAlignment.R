library(ape)

file_list <- list.files(path = "/Users/jamesonhinkle/Desktop/folders/Longwood/Research/eDNA/R/Longwood_eDNA/Sanger_BarcodeAlign_test_GoodRead", 
                        pattern = "\\.seq$", 
                        full.names = TRUE, 
                        ignore.case = TRUE)

# --- SAFETY CHECK ---
if (length(file_list) == 0) {
  stop(paste("No .seq files found in:", seq_folder, 
             "\nCheck your path and file extensions!"))
} else {
  message(paste("Found", length(file_list), "files. Starting import..."))
}
# ---------------------

# 3. Read and Combine
# We'll use a loop to make it easier to debug which specific file might fail
dna_list <- lapply(file_list, function(f) {
  tryCatch({
    readDNAStringSet(f)
  }, error = function(e) {
    stop(paste("Error reading file:", f, "\nDetails:", e$message))
  })
})

dna_set <- do.call(c, dna_list)

names(dna_set) <- gsub("\\.seq$", "", basename(file_list), ignore.case = TRUE)

alignment <- AlignSeqs(dna_set)

BrowseSeqs(alignment)

# get consensus sequence

consensus_dna <- ConsensusSequence(alignment, 
                                   threshold = 0.6, 
                                   noConsensusChar = "N")

writeXStringSet(consensus_dna, file = "Nine_Sanger_results_ForBlast.fasta")

BiocManager::install("annotate")
library(annotate)

blast_result <- blastSequences(consensus_dna, as="data.frame") # R port into BLAST
