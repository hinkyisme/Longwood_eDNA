# ============================================================
# server.R  –  eDNA Wiki - Build Barcodes, Primers, and Identify Taxa.
# ============================================================

library(shiny)
library(DECIPHER)
library(Biostrings)
library(DT)
library(tools)
library(ape)
library(stringr)
library(httr)
library(rvest)

# ============================================================
# PRIMER DESIGN ENGINE
# All functions are self-contained and dependency-free beyond base R.
# ============================================================

# --- Nearest-Neighbor Tm (SantaLucia 1998, 50 mM Na+, 250 nM primer) ---
calc_tm <- function(seq) {
  seq <- toupper(seq)
  n   <- nchar(seq)
  if (n < 2) return(NA_real_)

  # NN parameters: dH (kcal/mol), dS (cal/mol/K)
  nn_dH <- c(AA=7.9,  AT=7.2,  AC=8.4,  AG=7.8,
             TA=7.2,  TT=7.9,  TC=8.2,  TG=8.5,
             CA=8.5,  CT=7.8,  CC=8.0,  CG=10.6,
             GA=8.2,  GT=8.4,  GC=9.8,  GG=8.0)
  nn_dS <- c(AA=22.2, AT=20.4, AC=22.4, AG=21.0,
             TA=21.3, TT=22.2, TC=22.2, TG=22.7,
             CA=22.7, CT=21.0, CC=19.9, CG=27.2,
             GA=22.2, GT=22.4, GC=24.4, GG=19.9)

  dH_sum <- 0; dS_sum <- 0
  for (i in seq_len(n - 1)) {
    dinuc <- substr(seq, i, i + 1)
    dH_sum <- dH_sum + (nn_dH[dinuc] %||% 8.0)   # fallback for ambiguous bases
    dS_sum <- dS_sum + (nn_dS[dinuc] %||% 22.0)
  }
  # initiation corrections (terminal GC vs AT)
  init_dH <- ifelse(substr(seq,1,1) %in% c("G","C"), 0.1, 2.3)
  init_dS <- ifelse(substr(seq,1,1) %in% c("G","C"), -2.8, 4.1)
  init_dH <- init_dH + ifelse(substr(seq,n,n) %in% c("G","C"), 0.1, 2.3)
  init_dS <- init_dS + ifelse(substr(seq,n,n) %in% c("G","C"), -2.8, 4.1)

  dH_total <- (dH_sum + init_dH) * 1000   # convert to cal/mol
  dS_total  <- dS_sum + init_dS
  R  <- 1.987        # cal / mol / K
  CT <- 250e-9       # primer concentration (M)
  Tm <- dH_total / (dS_total + R * log(CT / 4)) - 273.15
  round(Tm, 1)
}

# Null-coalescing helper (base R doesn't have %||%)
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

# --- GC content ---
gc_content <- function(seq) {
  seq <- toupper(seq)
  gc  <- str_count(seq, "[GC]")
  round(gc / nchar(seq) * 100, 1)
}

# --- Reverse complement ---
rev_comp <- function(seq) {
  chartr("ACGT", "TGCA", paste(rev(strsplit(toupper(seq), "")[[1]]), collapse=""))
}

# --- 3'-end GC clamp check (last 5 bases should have 1–3 GC) ---
gc_clamp_ok <- function(seq) {
  tail5 <- substr(seq, nchar(seq)-4, nchar(seq))
  n_gc  <- str_count(tail5, "[GC]")
  n_gc >= 1 && n_gc <= 3
}

# --- Homopolymer run check (no run > 4 identical bases) ---
no_long_runs <- function(seq) {
  !grepl("(A{5,}|T{5,}|G{5,}|C{5,})", toupper(seq))
}

# --- Simple hairpin / self-dimer penalty (3'-end 6-mer vs full primer) ---
hairpin_penalty <- function(seq) {
  seq <- toupper(seq)
  n   <- nchar(seq)
  if (n < 12) return(0)
  tail6 <- substr(seq, n-5, n)
  tail6_rc <- rev_comp(tail6)
  # count complementary base pairs in any 6-window of the primer
  max_bp <- 0
  for (i in seq_len(n - 5)) {
    window <- substr(seq, i, i+5)
    bp <- sum(strsplit(window,"")[[1]] == strsplit(tail6_rc,"")[[1]])
    if (bp > max_bp) max_bp <- bp
  }
  max_bp   # 0–6; penalise if > 3
}

# --- Score a single candidate primer (lower = better) ---
score_primer <- function(seq) {
  tm  <- calc_tm(seq)
  gc  <- gc_content(seq)
  penalty <- 0

  # Tm range 55–65°C ideal
  if (is.na(tm))          penalty <- penalty + 50
  else if (tm < 50)       penalty <- penalty + (50 - tm) * 2
  else if (tm > 68)       penalty <- penalty + (tm - 68) * 2
  else if (tm < 55)       penalty <- penalty + (55 - tm)
  else if (tm > 65)       penalty <- penalty + (tm - 65)

  # GC range 40–60%
  if (gc < 40)            penalty <- penalty + (40 - gc) * 0.5
  else if (gc > 60)       penalty <- penalty + (gc - 60) * 0.5

  # GC clamp
  if (!gc_clamp_ok(seq))  penalty <- penalty + 5

  # Homopolymer runs
  if (!no_long_runs(seq)) penalty <- penalty + 10

  # Hairpin / self-dimer
  hp <- hairpin_penalty(seq)
  if (hp > 3)             penalty <- penalty + (hp - 3) * 4

  list(tm = tm, gc = gc, penalty = penalty)
}

# --- Enumerate all candidate primers in one strand ---
# Returns a data.frame of candidates passing hard filters
enumerate_candidates <- function(seq, kmin = 18, kmax = 24,
                                 tm_min = 50, tm_max = 70,
                                 gc_min = 35, gc_max = 65) {
  seq  <- toupper(seq)
  n    <- nchar(seq)
  rows <- list()

  for (k in kmin:kmax) {
    for (start in seq_len(n - k + 1)) {
      subseq <- substr(seq, start, start + k - 1)
      # Quick pre-filter: skip if ambiguous bases > 0
      if (grepl("[^ACGT]", subseq)) next
      gc  <- gc_content(subseq)
      if (gc < gc_min || gc > gc_max) next
      tm  <- calc_tm(subseq)
      if (is.na(tm) || tm < tm_min || tm > tm_max) next
      if (!no_long_runs(subseq)) next
      sc  <- score_primer(subseq)
      rows[[length(rows)+1]] <- data.frame(
        start   = start,
        end     = start + k - 1,
        length  = k,
        seq     = subseq,
        tm      = sc$tm,
        gc      = sc$gc,
        penalty = sc$penalty,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# --- Pair forward + reverse candidates into primer pairs ---
design_primer_pairs <- function(barcode_seq,
                                kmin = 18, kmax = 24,
                                amp_min = 100, amp_max = 600,
                                tm_delta_max = 5,
                                n_pairs = 20) {

  barcode_seq <- toupper(barcode_seq)
  n           <- nchar(barcode_seq)

  # Forward candidates on sense strand
  fwd_cands <- enumerate_candidates(barcode_seq, kmin, kmax)

  # Reverse candidates on antisense strand (enumerate on RC, then map coords back)
  rc_seq    <- rev_comp(barcode_seq)
  rev_cands <- enumerate_candidates(rc_seq, kmin, kmax)
  if (!is.null(rev_cands)) {
    # Convert RC coords back to sense-strand coords
    rev_cands$bind_end   <- n - rev_cands$start + 1
    rev_cands$bind_start <- n - rev_cands$end   + 1
  }

  if (is.null(fwd_cands) || is.null(rev_cands)) {
    return(data.frame(
      Pair = integer(0), Fwd_Seq = character(0), Rev_Seq = character(0),
      Fwd_Tm = numeric(0), Rev_Tm = numeric(0),
      Fwd_GC = numeric(0), Rev_GC = numeric(0),
      Amplicon_bp = integer(0), Pair_Score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # Keep top-100 candidates per strand by penalty to limit O(n²) pairing
  fwd_top <- head(fwd_cands[order(fwd_cands$penalty), ], 100)
  rev_top <- head(rev_cands[order(rev_cands$penalty), ], 100)

  pairs <- list()
  for (i in seq_len(nrow(fwd_top))) {
    for (j in seq_len(nrow(rev_top))) {
      f <- fwd_top[i, ]
      r <- rev_top[j, ]
      amp <- r$bind_end - f$start + 1
      if (amp < amp_min || amp > amp_max) next
      if (abs(f$tm - r$tm) > tm_delta_max) next
      # Combined score: sum of individual penalties + Tm delta bonus
      pair_score <- f$penalty + r$penalty + abs(f$tm - r$tm)
      pairs[[length(pairs)+1]] <- data.frame(
        Fwd_Seq     = f$seq,
        Rev_Seq     = r$seq,
        Fwd_Start   = f$start,
        Fwd_End     = f$end,
        Rev_Bind_Start = r$bind_start,
        Rev_Bind_End   = r$bind_end,
        Fwd_Tm      = f$tm,
        Rev_Tm      = r$tm,
        Fwd_GC      = f$gc,
        Rev_GC      = r$gc,
        Amplicon_bp = as.integer(amp),
        Pair_Score  = round(pair_score, 2),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(pairs) == 0) {
    return(data.frame(
      Pair = integer(0), Fwd_Seq = character(0), Rev_Seq = character(0),
      Fwd_Tm = numeric(0), Rev_Tm = numeric(0),
      Fwd_GC = numeric(0), Rev_GC = numeric(0),
      Amplicon_bp = integer(0), Pair_Score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  result <- do.call(rbind, pairs)
  result <- result[order(result$Pair_Score), ]

  # Deduplicate: drop pairs too similar to a higher-ranked pair (offset < 5 bp)
  keep <- rep(TRUE, nrow(result))
  for (i in seq_len(nrow(result) - 1)) {
    if (!keep[i]) next
    for (j in (i+1):nrow(result)) {
      if (!keep[j]) next
      if (abs(result$Fwd_Start[i] - result$Fwd_Start[j]) < 5 &&
          abs(result$Rev_Bind_End[i] - result$Rev_Bind_End[j]) < 5) {
        keep[j] <- FALSE
      }
    }
  }
  result <- result[keep, ]
  result <- head(result, n_pairs)
  result$Pair <- seq_len(nrow(result))
  result[, c("Pair","Fwd_Seq","Rev_Seq","Fwd_Tm","Rev_Tm",
             "Fwd_GC","Rev_GC","Amplicon_bp","Pair_Score",
             "Fwd_Start","Fwd_End","Rev_Bind_Start","Rev_Bind_End")]
}

# --- Build HTML sequence viewer with primer pair highlighted ---
highlight_barcode <- function(barcode_seq, fwd_start, fwd_end,
                              rev_bind_start, rev_bind_end,
                              line_width = 60) {
  seq  <- toupper(barcode_seq)
  n    <- nchar(seq)
  chars <- strsplit(seq, "")[[1]]

  styled <- character(n)
  for (i in seq_len(n)) {
    base <- chars[i]
    col  <- switch(base, A="#2196F3", T="#FF5722", G="#4CAF50", C="#FF9800", "#999")
    cls  <- ""
    if (i >= fwd_start && i <= fwd_end)           cls <- " fwd-primer"
    if (i >= rev_bind_start && i <= rev_bind_end)  cls <- " rev-primer"
    styled[i] <- sprintf('<span style="color:%s" class="base%s">%s</span>', col, cls, base)
  }

  # Wrap into lines with position numbers
  lines <- character(0)
  pos   <- 1
  while (pos <= n) {
    chunk_end <- min(pos + line_width - 1, n)
    chunk     <- paste(styled[pos:chunk_end], collapse="")
    lines <- c(lines, sprintf('<span class="pos-num">%5d</span> %s', pos, chunk))
    pos   <- chunk_end + 1
  }

  paste0(
    '<div class="seq-viewer">',
    '<div class="seq-legend">',
    '<span class="leg fwd-leg">&#9632; Forward primer</span>',
    '<span class="leg rev-leg">&#9632; Reverse primer (binding site)</span>',
    '</div>',
    '<pre class="seq-pre">', paste(lines, collapse="\n"), '</pre>',
    '</div>'
  )
}

# ============================================================
# GENBANK 16S EXTRACTOR (from Barcode_builder.R, adapted for Shiny)
# ============================================================

extract_16S_from_gb <- function(gb_file, verbose = FALSE) {
  rc <- function(seq) chartr("ACGT","TGCA", paste(rev(strsplit(toupper(seq),"")[[1]]),collapse=""))

  gb_lines <- readLines(gb_file, warn = FALSE)
  rec_ends <- which(gb_lines == "//")
  if (length(rec_ends) == 0) stop("No '//' record separators found.")

  seq_list <- list()
  rec_start <- 1

  for (rec_idx in seq_along(rec_ends)) {
    end_idx   <- rec_ends[rec_idx]
    rec_lines <- gb_lines[rec_start:(end_idx - 1)]
    rec_start <- end_idx + 1
    if (length(rec_lines) == 0) next

    locus_line     <- grep("^LOCUS",     rec_lines, value=TRUE)
    accession_line <- grep("^ACCESSION", rec_lines, value=TRUE)
    locus_name <- if (length(locus_line)     > 0) strsplit(locus_line,     "\\s+")[[1]][2] else
                  if (length(accession_line) > 0) strsplit(accession_line, "\\s+")[[1]][2] else
                  paste0("Genome_", rec_idx)

    origin_idx <- grep("^ORIGIN", rec_lines)
    if (length(origin_idx) == 0) next
    dna_seq <- toupper(gsub("[0-9\\s]","", paste(rec_lines[(origin_idx+1):length(rec_lines)], collapse="")))
    if (nchar(dna_seq) == 0) next

    feat_idx <- grep("^FEATURES", rec_lines)
    if (length(feat_idx) == 0) next
    feat_lines   <- rec_lines[(feat_idx+1):(origin_idx-1)]
    entry_start  <- which(grepl("^\\s{5}\\S", feat_lines, perl=TRUE))
    if (length(entry_start) == 0) next
    entry_ranges <- c(entry_start, length(feat_lines)+1)
    feature_entries <- lapply(seq_along(entry_start), function(k) {
      paste(feat_lines[entry_start[k]:(entry_ranges[k+1]-1)], collapse=" ")
    })

    hits <- which(sapply(feature_entries, function(x)
      grepl("16S|16S ribosomal|16S rRNA", x, ignore.case=TRUE)))
    if (length(hits)==0)
      hits <- which(sapply(feature_entries, function(x)
        grepl("rRNA",x,ignore.case=TRUE) && grepl("16S|16s",x,ignore.case=TRUE)))
    if (length(hits)==0) next

    seq_count <- 0
    for (h in hits) {
      entry  <- feature_entries[[h]]
      coords <- regmatches(entry, gregexpr("[0-9]+\\.{2}[0-9]+", entry, perl=TRUE))[[1]]
      if (length(coords)==0) next
      is_comp <- grepl("complement\\s*\\(", entry, ignore.case=TRUE)
      parts   <- sapply(coords, function(coord) {
        se <- as.numeric(strsplit(coord,"\\.\\.")[[1]])
        if (length(se)!=2) return("")
        se[1] <- max(se[1],1); se[2] <- min(se[2], nchar(dna_seq))
        substr(dna_seq, se[1], se[2])
      }, USE.NAMES=FALSE)
      seq_string <- gsub("\\s+","", paste(parts, collapse=""))
      if (nchar(seq_string)==0) next
      if (is_comp) seq_string <- rc(seq_string)
      seq_count <- seq_count + 1
      nm <- if (seq_count==1) locus_name else paste0(locus_name,"_16Scopy",seq_count)
      while (nm %in% names(seq_list)) nm <- paste0(nm,"_dup")
      seq_list[[nm]] <- seq_string
    }
  }
  seq_list
}

# ============================================================
# PRIMER-BLAST NCBI INTEGRATION
#
# Uses NCBI's Primer-BLAST CGI endpoint (same one the web form posts to).
# Step 1 – POST job: returns a job_key / request_id in the response URL.
# Step 2 – Poll GET until status is no longer "Job is running".
# Step 3 – Parse the results HTML table with rvest.
#
# NCBI fair-use: max 3 requests/sec; we add a 1 s delay before polling.
# The organism field defaults to "all" (taxid 1) so results reflect
# genome-wide specificity; users can override via the UI.
# ============================================================

PRIMER_BLAST_URL <- "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"

# POST a Primer-BLAST job and return the polling URL + job_key
primer_blast_submit <- function(fwd_seq, rev_seq,
                                organism   = "Mytilus",   # genus or taxon name
                                db         = "nt",
                                max_hits   = 20) {

  body <- list(
    INPUT_SEQUENCE   = "",          # no template – we supply the primers directly
    PRIMER5_INPUT    = fwd_seq,
    PRIMER3_INPUT    = rev_seq,
    SEARCHMODE       = "1",         # "1" = Check primers, no design
    PRIMER_PRODUCT_MIN = "50",
    PRIMER_PRODUCT_MAX = "2000",
    ORGANISM         = organism,
    PRIMER_SPECIFICITY_DATABASE = db,
    HITSIZE          = as.character(max_hits),
    NUM_TARGETS_WITH_PRIMERS = "20",
    NUM_TARGETS      = "20",
    SHOW_SSEQ        = "on",
    PRIMER_LEFT_TM   = "60",
    PRIMER_RIGHT_TM  = "60",
    CMD              = "web"
  )

  resp <- POST(PRIMER_BLAST_URL,
               body   = body,
               encode = "form",
               add_headers(
                 `User-Agent` = "R/Shiny eDNA-pipeline (educational use)"
               ),
               timeout(30))

  if (http_error(resp)) {
    stop(paste("Primer-BLAST submission failed:", status_code(resp)))
  }

  # The redirect URL contains job_key=<ID>
  final_url <- resp$url
  if (is.null(final_url) || !grepl("job_key=", final_url)) {
    # Try to extract from response body
    body_text <- content(resp, "text", encoding = "UTF-8")
    m <- regmatches(body_text, regexpr("job_key=[A-Za-z0-9_]+", body_text))
    if (length(m) == 0) stop("Could not extract Primer-BLAST job_key from response.")
    final_url <- paste0(PRIMER_BLAST_URL, "?", m, "&CMD=web")
  }
  final_url
}

# Poll the result URL until the job finishes (or timeout)
# Returns the raw HTML content string
primer_blast_poll <- function(result_url,
                              max_wait_sec = 120,
                              poll_interval = 5) {
  Sys.sleep(1)   # brief courtesy pause before first poll
  elapsed <- 0

  repeat {
    resp <- GET(result_url,
                add_headers(`User-Agent` = "R/Shiny eDNA-pipeline (educational use)"),
                timeout(20))
    if (http_error(resp)) stop(paste("Poll request failed:", status_code(resp)))
    html_text <- content(resp, "text", encoding = "UTF-8")

    # NCBI signals "still running" with this phrase
    still_running <- grepl("Job is running|job_status.*RUNNING|Please wait", html_text, ignore.case = TRUE)

    if (!still_running) return(list(html = html_text, url = result_url))

    elapsed <- elapsed + poll_interval
    if (elapsed >= max_wait_sec) stop("Primer-BLAST job timed out after ", max_wait_sec, " seconds.")
    Sys.sleep(poll_interval)
  }
}

# Parse the finished HTML page into a tidy data.frame of BLAST hits
primer_blast_parse <- function(html_text) {
  page <- tryCatch(read_html(html_text), error = function(e) NULL)
  if (is.null(page)) return(NULL)

  # ── Check for "no hits" message ──────────────────────────────────────────
  no_hits_msg <- html_text |>
    (\(h) grepl("No hits found|0 targets|could not find any", h, ignore.case = TRUE))()
  if (no_hits_msg) {
    return(data.frame(
      Result = "No off-target hits found — primers appear specific.",
      stringsAsFactors = FALSE
    ))
  }

  # ── Try to parse the primer pair / product table ──────────────────────────
  # NCBI renders results in a <div class="prPairInfo"> or a plain <table>
  tables <- html_table(page, fill = TRUE)

  # Filter to tables that look like BLAST hit tables (have accession-like content)
  hit_tables <- Filter(function(tbl) {
    ncol(tbl) >= 3 && nrow(tbl) >= 1 &&
      any(grepl("NM_|NR_|XM_|NC_|[A-Z]{1,2}[0-9]{5,}", as.character(tbl), ignore.case = FALSE))
  }, tables)

  if (length(hit_tables) == 0) {
    # Fall back: return a simple "check URL" message
    return(data.frame(
      Result = "Results received — please view the full report via the NCBI link above.",
      stringsAsFactors = FALSE
    ))
  }

  # Take the largest hit table and clean it up
  tbl <- hit_tables[[which.max(sapply(hit_tables, nrow))]]
  # Drop entirely-NA columns
  tbl <- tbl[, colSums(is.na(tbl)) < nrow(tbl)]
  # Standardise column names
  colnames(tbl) <- make.names(colnames(tbl), unique = TRUE)
  tbl
}

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # --- Utility: read any supported DNA format ---
  smart_read <- function(file_path) {
    ext <- tolower(file_ext(file_path))
    if (ext %in% c("fastq","fq")) {
      return(as(readQualityScaledDNAStringSet(file_path), "DNAStringSet"))
    } else if (ext == "seq") {
      lines <- readLines(file_path)
      if (any(grepl("^>", lines))) return(readDNAStringSet(file_path))
      return(DNAStringSet(paste(lines, collapse="")))
    } else {
      return(readDNAStringSet(file_path))
    }
  }

  # ── TAXONOMY PIPELINE ───────────────────────────────────────
  results <- eventReactive(input$run, {
    req(input$ref_fasta, input$test_files)
    withProgress(message = "Initializing Pipeline", value = 0, {
      incProgress(0.2, detail = "Importing reference genomes...")
      ref_dna <- smart_read(input$ref_fasta$datapath)

      incProgress(0.2, detail = "Extracting 16S barcode regions...")
      extracted_ref <- TrimDNA(ref_dna,
                               left  = input$fwd,
                               right = input$rev,
                               type  = "sequences",
                               max.mismatch = input$mismatch)

      incProgress(0.3, detail = "Training taxonomic classifier...")
      train_labels  <- paste0("Root;Mitochondria;", names(extracted_ref))
      trained_model <- LearnTaxa(extracted_ref, train_labels, maxKmerSize = 7)

      incProgress(0.2, detail = "Processing test sequences...")
      test_dna <- smart_read(input$test_files$datapath[1])
      ids      <- IdTaxa(test_dna, trained_model, strand="both", threshold=input$threshold)

      incProgress(0.1, detail = "Cleaning up memory...")
      gc()
      list(ids=ids, test_dna=test_dna)
    })
  })

  output$results_table <- renderDT({
    data <- results()
    taxa <- sapply(data$ids, function(x) paste(x$taxon,      collapse=" > "))
    conf <- sapply(data$ids, function(x) paste(round(x$confidence,1), collapse=" | "))
    df   <- data.frame(Sequence_ID=names(data$test_dna), Taxonomy_Path=taxa, Confidence_Scores=conf)
    datatable(df, options=list(pageLength=10, scrollX=TRUE))
  })

  # ── BARCODE BUILDER ─────────────────────────────────────────

  # Reactive: extract 16S sequences from uploaded GenBank file
  barcode_seqs <- eventReactive(input$run_barcode, {
    req(input$gb_file)
    withProgress(message="Extracting 16S barcodes...", value=0.3, {
      seqs <- tryCatch(
        extract_16S_from_gb(input$gb_file$datapath, verbose=FALSE),
        error = function(e) { showNotification(paste("Extraction error:", e$message), type="error"); NULL }
      )
      incProgress(0.7)
      seqs
    })
  })

  # Reactive: compute consensus barcode from extracted sequences
  consensus_barcode <- reactive({
    seqs <- barcode_seqs()
    req(seqs, length(seqs) > 0)
    dna_set <- DNAStringSet(unlist(seqs))
    if (length(dna_set) == 1) return(as.character(dna_set[[1]]))
    # Build consensus via DECIPHER alignment + consensusString
    aligned   <- AlignSeqs(dna_set, verbose=FALSE)
    cons_raw  <- consensusString(aligned, ambiguityMap="N", threshold=0.6)
    # Strip gap characters
    gsub("-","", cons_raw)
  })

  # Summary box: extraction stats
  output$barcode_summary <- renderUI({
    seqs <- barcode_seqs()
    if (is.null(seqs) || length(seqs)==0) {
      return(tags$p("No sequences extracted yet.", style="color:#888;"))
    }
    lens <- nchar(unlist(seqs))
    tagList(
      tags$p(tags$b(length(seqs)), " 16S sequences extracted"),
      tags$p("Length range: ", min(lens), "–", max(lens), " bp"),
      tags$p("Consensus length: ", nchar(consensus_barcode()), " bp")
    )
  })

  # Table of extracted barcode names
  output$barcode_table <- renderDT({
    seqs <- barcode_seqs()
    req(seqs)
    df <- data.frame(
      Name   = names(seqs),
      Length = nchar(unlist(seqs)),
      stringsAsFactors = FALSE
    )
    datatable(df, options=list(pageLength=10, scrollX=TRUE, dom="ftp"), rownames=FALSE)
  })

  # Download extracted barcodes as FASTA
  output$download_barcodes <- downloadHandler(
    filename = function() paste0("16S_barcodes_", Sys.Date(), ".fasta"),
    content  = function(file) {
      seqs <- barcode_seqs()
      req(seqs)
      lines <- unlist(lapply(names(seqs), function(nm)
        c(paste0(">", nm), strwrap(seqs[[nm]], width=80))
      ))
      writeLines(lines, file)
    }
  )

  # ── PRIMER DESIGNER ─────────────────────────────────────────

  # Reactive: run primer design on consensus barcode (or user-pasted sequence)
  primer_results <- eventReactive(input$run_primers, {
    # Priority: pasted sequence > consensus from GenBank upload
    target_seq <- trimws(input$custom_barcode)
    if (nchar(target_seq) < 50) {
      target_seq <- tryCatch(consensus_barcode(), error=function(e) NULL)
    }
    validate(need(!is.null(target_seq) && nchar(target_seq) >= 50,
      "Please upload a GenBank file and run the Barcode Builder, or paste a barcode sequence (≥50 bp) above."))

    withProgress(message="Designing primers...", value=0.1, {
      incProgress(0.4, detail="Enumerating candidates...")
      pairs <- design_primer_pairs(
        barcode_seq  = target_seq,
        kmin         = input$primer_kmin,
        kmax         = input$primer_kmax,
        amp_min      = input$amp_min,
        amp_max      = input$amp_max,
        tm_delta_max = input$tm_delta,
        n_pairs      = 25
      )
      incProgress(0.5, detail="Ranking pairs...")
      list(pairs=pairs, target_seq=target_seq)
    })
  })

  # Primer results table
  output$primer_table <- renderDT({
    data <- primer_results()
    df   <- data$pairs
    if (nrow(df)==0) {
      return(datatable(data.frame(Message="No primer pairs found with current settings. Try relaxing amplicon size or Tm delta."),
                       options=list(dom="t"), rownames=FALSE))
    }
    # Display subset of columns
    display <- df[, c("Pair","Fwd_Seq","Rev_Seq","Fwd_Tm","Rev_Tm","Fwd_GC","Rev_GC","Amplicon_bp","Pair_Score")]
    colnames(display) <- c("Pair","Forward Primer (5'→3')","Reverse Primer (5'→3')",
                           "Fwd Tm (°C)","Rev Tm (°C)","Fwd GC%","Rev GC%","Amplicon (bp)","Score")
    datatable(
      display,
      selection = "single",
      rownames  = FALSE,
      options   = list(pageLength=10, scrollX=TRUE,
                       columnDefs=list(list(className="dt-center", targets=c(3,4,5,6,7,8)))),
      class     = "stripe hover compact"
    ) |>
      formatRound(columns=c("Fwd Tm (°C)","Rev Tm (°C)","Fwd GC%","Rev GC%","Score"), digits=1)
  })

  # Sequence viewer: highlight selected primer pair on barcode
  output$seq_viewer <- renderUI({
    data  <- primer_results()
    df    <- data$pairs
    s_row <- input$primer_table_rows_selected
    if (is.null(s_row) || nrow(df)==0) {
      return(tags$p("Select a primer pair in the table above to highlight it on the barcode sequence.",
                    style="color:#888; font-style:italic;"))
    }
    pair <- df[s_row, ]
    html <- highlight_barcode(
      data$target_seq,
      fwd_start      = pair$Fwd_Start,
      fwd_end        = pair$Fwd_End,
      rev_bind_start = pair$Rev_Bind_Start,
      rev_bind_end   = pair$Rev_Bind_End
    )
    HTML(html)
  })

  # Primer detail card for selected pair
  output$primer_detail <- renderUI({
    data  <- primer_results()
    df    <- data$pairs
    s_row <- input$primer_table_rows_selected
    if (is.null(s_row) || nrow(df)==0) return(NULL)
    p <- df[s_row, ]
    tagList(
      tags$hr(),
      tags$h5("Selected Pair Detail"),
      tags$table(class="table table-sm table-bordered", style="max-width:600px;",
        tags$thead(tags$tr(tags$th(""), tags$th("Forward"), tags$th("Reverse"))),
        tags$tbody(
          tags$tr(tags$td("Sequence (5'→3')"),
                  tags$td(tags$code(p$Fwd_Seq)),
                  tags$td(tags$code(p$Rev_Seq))),
          tags$tr(tags$td("Length (bp)"),
                  tags$td(nchar(p$Fwd_Seq)), tags$td(nchar(p$Rev_Seq))),
          tags$tr(tags$td("Tm (°C)"),
                  tags$td(p$Fwd_Tm), tags$td(p$Rev_Tm)),
          tags$tr(tags$td("GC%"),
                  tags$td(p$Fwd_GC), tags$td(p$Rev_GC)),
          tags$tr(tags$td(tags$b("Amplicon")),
                  tags$td(colspan="2", paste(p$Amplicon_bp, "bp")))
        )
      )
    )
  })

  # ── PRIMER-BLAST CHECK ──────────────────────────────────────

  # Reactive value to hold PrimerBLAST results (keyed by pair index)
  blast_result <- reactiveVal(NULL)
  blast_result_url <- reactiveVal(NULL)

  # Fire when user clicks "Check via Primer-BLAST"
  observeEvent(input$run_primerblast, {
    data  <- primer_results()
    df    <- data$pairs
    s_row <- input$primer_table_rows_selected

    if (is.null(s_row) || nrow(df) == 0) {
      showNotification("Please select a primer pair from the table first.", type = "warning")
      return()
    }

    pair <- df[s_row, ]

    # Validate sequences look like DNA
    if (!grepl("^[ACGT]+$", pair$Fwd_Seq, ignore.case = TRUE) ||
        !grepl("^[ACGT]+$", pair$Rev_Seq, ignore.case = TRUE)) {
      showNotification("Primer sequences contain ambiguous bases — cannot submit to Primer-BLAST.", type = "error")
      return()
    }

    blast_result(NULL)      # clear previous result
    blast_result_url(NULL)

    withProgress(message = "Submitting to NCBI Primer-BLAST…", value = 0.1, {

      incProgress(0.2, detail = "Posting job…")
      result_url <- tryCatch(
        primer_blast_submit(
          fwd_seq  = pair$Fwd_Seq,
          rev_seq  = pair$Rev_Seq,
          organism = input$blast_organism,
          db       = input$blast_db
        ),
        error = function(e) {
          showNotification(paste("Submission error:", e$message), type = "error", duration = 10)
          NULL
        }
      )

      if (is.null(result_url)) return()
      blast_result_url(result_url)

      incProgress(0.4, detail = "Waiting for NCBI results (may take 30–90 s)…")
      raw <- tryCatch(
        primer_blast_poll(result_url, max_wait_sec = 180, poll_interval = 6),
        error = function(e) {
          showNotification(paste("Poll error:", e$message), type = "error", duration = 10)
          NULL
        }
      )

      if (is.null(raw)) return()

      incProgress(0.3, detail = "Parsing results…")
      parsed <- primer_blast_parse(raw$html)
      blast_result(list(table = parsed, url = raw$url, pair = pair))
    })
  })

  # Render the PrimerBLAST results panel
  output$blast_result_panel <- renderUI({
    res <- blast_result()
    if (is.null(res)) return(
      tags$p("Select a primer pair and click Check via Primer-BLAST.",
             style = "color:#888; font-style:italic;")
    )

    pair <- res$pair
    tagList(
      # Direct link to NCBI results
      tags$div(style = "margin-bottom:10px;",
        tags$a(href   = res$url,
               target = "_blank",
               class  = "btn btn-sm btn-info",
               tags$span(class="glyphicon glyphicon-new-window"),
               " View Full Results on NCBI Primer-BLAST")
      ),
      # Which pair was checked
      tags$p(tags$b("Checked pair: "),
             tags$code(pair$Fwd_Seq), " / ", tags$code(pair$Rev_Seq)),
      # Hits table or message
      if (!is.null(res$table) && "Result" %in% colnames(res$table)) {
        tags$div(class = "alert alert-success", res$table$Result[1])
      } else if (!is.null(res$table) && nrow(res$table) > 0) {
        tagList(
          tags$p(tags$b(nrow(res$table)), " potential amplification product(s) found across the selected database:"),
          DTOutput("blast_hits_table")
        )
      } else {
        tags$div(class = "alert alert-warning",
                 "Results were received but could not be parsed automatically. Use the NCBI link above.")
      }
    )
  })

  output$blast_hits_table <- renderDT({
    res <- blast_result()
    req(res, !is.null(res$table), !"Result" %in% colnames(res$table))
    datatable(res$table,
              options  = list(pageLength = 10, scrollX = TRUE, dom = "ftp"),
              rownames = FALSE,
              class    = "stripe hover compact")
  })

  # ── DOWNLOAD BUTTON ──────────────────────────────────────────

  # Render download button only after primers have been computed
  output$download_primers_btn <- renderUI({
    req(primer_results())
    downloadButton("download_primers", "Download Primer CSV",
                   class = "btn-sm btn-default", style = "width:100%;")
  })

  # Download primer pairs as CSV
  output$download_primers <- downloadHandler(
    filename = function() paste0("primer_pairs_", Sys.Date(), ".csv"),
    content  = function(file) {
      data <- primer_results()
      write.csv(data$pairs[, c("Pair","Fwd_Seq","Rev_Seq","Fwd_Tm","Rev_Tm",
                               "Fwd_GC","Rev_GC","Amplicon_bp","Pair_Score")],
                file, row.names=FALSE)
    }
  )

}
