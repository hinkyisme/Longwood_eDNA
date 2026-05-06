# ============================================================
# server.R  —  eDNA Barcode Builder & Identifier
# ============================================================

library(shiny)
library(DECIPHER)
library(Biostrings)
library(DT)
library(tools)
library(ape)

# ── Helper: reverse complement (base R, no Bioconductor) ──────
rc_seq <- function(seq) {
  seq <- toupper(seq)
  seq_chars <- rev(strsplit(seq, "")[[1]])
  chartr("ACGT", "TGCA", paste(seq_chars, collapse = ""))
}

# ── Helper: smart multi-format reader ─────────────────────────
smart_read <- function(file_path) {
  ext <- tolower(file_ext(file_path))
  if (ext %in% c("fastq", "fq")) {
    return(as(readQualityScaledDNAStringSet(file_path), "DNAStringSet"))
  } else if (ext == "seq") {
    lines <- readLines(file_path, warn = FALSE)
    if (any(grepl("^>", lines))) return(readDNAStringSet(file_path))
    dna <- DNAStringSet(paste(grep("^[ACGTNacgtn]+", lines, value = TRUE), collapse = ""))
    return(dna)
  } else {
    return(readDNAStringSet(file_path))
  }
}

# ── Helper: extract 16S from GenBank flat file ────────────────
extract_16S_from_gb <- function(gb_file, verbose_fn = message) {
  gb_lines  <- readLines(gb_file, warn = FALSE)
  rec_ends  <- which(gb_lines == "//")
  if (length(rec_ends) == 0) stop("No record separators ('//') found in file.")

  seq_list  <- list()
  rec_start <- 1
  rec_idx   <- 0

  for (end_idx in rec_ends) {
    rec_idx   <- rec_idx + 1
    rec_lines <- gb_lines[rec_start:(end_idx - 1)]
    rec_start <- end_idx + 1
    if (length(rec_lines) == 0) next

    # Accession / locus name
    locus_line     <- grep("^LOCUS",     rec_lines, value = TRUE)
    accession_line <- grep("^ACCESSION", rec_lines, value = TRUE)
    locus_name <- if (length(locus_line) > 0) {
      strsplit(locus_line, "\\s+")[[1]][2]
    } else if (length(accession_line) > 0) {
      strsplit(accession_line, "\\s+")[[1]][2]
    } else {
      paste0("Genome_", rec_idx)
    }

    # Organism name (for better labels)
    org_lines <- grep("/organism=", rec_lines, value = TRUE)
    org_name  <- if (length(org_lines) > 0) {
      gsub('.*?/organism="([^"]+)".*', "\\1", org_lines[1])
    } else {
      locus_name
    }

    # Build raw sequence from ORIGIN
    origin_idx <- grep("^ORIGIN", rec_lines)
    if (length(origin_idx) == 0) next
    dna_lines <- rec_lines[(origin_idx + 1):length(rec_lines)]
    dna_seq   <- toupper(gsub("[0-9\\s /]", "", paste(dna_lines, collapse = "")))
    if (nchar(dna_seq) == 0) next

    # Features block
    feat_idx <- grep("^FEATURES", rec_lines)
    if (length(feat_idx) == 0) next
    feat_lines  <- rec_lines[(feat_idx + 1):(origin_idx - 1)]
    entry_start <- which(grepl("^\\s{5}\\S", feat_lines, perl = TRUE))
    if (length(entry_start) == 0) next

    entry_ranges  <- c(entry_start, length(feat_lines) + 1)
    feature_entries <- lapply(seq_along(entry_start), function(k) {
      paste(feat_lines[entry_start[k]:(entry_ranges[k + 1] - 1)], collapse = " ")
    })

    # Match 16S rRNA features
    hits <- which(sapply(feature_entries, function(x)
      grepl("16S|16s|16S ribosomal|16S rRNA", x, ignore.case = TRUE)))
    if (length(hits) == 0) {
      hits <- which(sapply(feature_entries, function(x)
        grepl("rRNA", x, ignore.case = TRUE) && grepl("16S|16s", x, ignore.case = TRUE)))
    }
    if (length(hits) == 0) {
      verbose_fn(paste("Record", rec_idx, "(", locus_name, ") — no 16S feature found"))
      next
    }

    seq_count <- 0
    for (h in hits) {
      entry  <- feature_entries[[h]]
      coords <- regmatches(entry, gregexpr("[0-9]+\\.{2}[0-9]+", entry, perl = TRUE))[[1]]
      if (length(coords) == 0) next

      is_comp <- grepl("complement\\s*\\(", entry, ignore.case = TRUE)

      parts <- sapply(coords, function(coord) {
        se <- as.numeric(unlist(strsplit(coord, "\\.\\.")))
        if (length(se) != 2) return("")
        se[1] <- max(1, se[1]); se[2] <- min(nchar(dna_seq), se[2])
        substr(dna_seq, se[1], se[2])
      }, USE.NAMES = FALSE)

      seq_string <- gsub("\\s+", "", paste(parts, collapse = ""))
      if (nchar(seq_string) == 0) next
      if (is_comp) seq_string <- rc_seq(seq_string)

      seq_count <- seq_count + 1
      seq_name  <- if (seq_count == 1) org_name else paste0(org_name, "_copy", seq_count)
      # uniquify
      base <- seq_name; sfx <- 1
      while (seq_name %in% names(seq_list)) { seq_name <- paste0(base, "_", sfx); sfx <- sfx + 1 }
      seq_list[[seq_name]] <- seq_string
    }
  }

  if (length(seq_list) == 0) stop("No 16S sequences extracted from GenBank file.")
  DNAStringSet(seq_list)
}

# ── SERVER ────────────────────────────────────────────────────
server <- function(input, output, session) {

  log_msgs <- reactiveVal(character(0))
  append_log <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste(...))
    log_msgs(c(log_msgs(), msg))
  }

  pipeline_state <- reactiveVal("idle")   # idle | running | done | error

  # ── Status bar ──────────────────────────────────────────────
  output$status_ui <- renderUI({
    s <- pipeline_state()
    dot_class <- switch(s, idle = "status-dot", running = "status-dot running", done = "status-dot ready", error = "status-dot")
    label <- switch(s, idle = "Awaiting pipeline execution", running = "Pipeline running…",
                    done = "Pipeline complete", error = "Pipeline error — see Log")
    div(class = "status-bar", span(class = dot_class), label)
  })

  # ── Reactive: full pipeline ──────────────────────────────────
  pipeline_result <- eventReactive(input$run, {
    log_msgs(character(0))
    pipeline_state("running")

    tryCatch({
      withProgress(message = "Running pipeline", value = 0, {

        # ---- STEP 1: Load / extract reference barcodes --------
        incProgress(0.15, detail = "Loading reference sequences…")
        append_log("=== STEP 1: Load reference sequences ===")

        mode <- isolate(input$input_mode)
        if (is.null(mode)) mode <- "genbank"

        ref_barcodes <- NULL
        extraction_info <- list(n = 0, names = character(0), lengths = integer(0))

        if (mode == "genbank") {
          req(input$gb_file)
          append_log("Mode: GenBank flat file →", input$gb_file$name)
          ref_barcodes <- extract_16S_from_gb(input$gb_file$datapath,
                                              verbose_fn = append_log)
          append_log("Extracted", length(ref_barcodes), "16S sequences from GenBank file.")
        } else {
          req(input$ref_fasta)
          append_log("Mode: FASTA/FASTQ reference →", input$ref_fasta$name)
          ref_dna <- smart_read(input$ref_fasta$datapath)
          append_log("Loaded", length(ref_dna), "reference sequences.")

          # primer-trim to isolate barcode region
          incProgress(0.1, detail = "Trimming primers from reference…")
          append_log("=== STEP 1b: Primer trimming ===")
          fwd_p <- isolate(input$fwd); rev_p <- isolate(input$rev)
          append_log("Forward:", fwd_p)
          append_log("Reverse:", rev_p)
          append_log("Max mismatches:", isolate(input$mismatch))

          f_hits <- vmatchPattern(fwd_p, ref_dna, max.mismatch = isolate(input$mismatch))
          r_hits <- vmatchPattern(reverseComplement(DNAString(rev_p)), ref_dna,
                                  max.mismatch = isolate(input$mismatch))

          valid_idx <- which(elementNROWS(f_hits) > 0 & elementNROWS(r_hits) > 0)
          append_log("Both primers found in", length(valid_idx), "of", length(ref_dna), "sequences.")

          if (length(valid_idx) == 0) stop("No reference sequences had both primers. Check primer sequences or increase mismatches.")

          barcode_list <- lapply(valid_idx, function(i) {
            s <- end(f_hits[[i]][1]) + 1
            e <- start(r_hits[[i]][1]) - 1
            if ((e - s) > 150 & (e - s) < 3000) return(subseq(ref_dna[[i]], start = s, end = e))
            return(NULL)
          })
          keep     <- !sapply(barcode_list, is.null)
          ref_barcodes <- DNAStringSet(unlist(barcode_list[keep]))
          names(ref_barcodes) <- names(ref_dna)[valid_idx[keep]]
          append_log("Retained", length(ref_barcodes), "barcodes after length filter.")
        }

        extraction_info <- list(
          n       = length(ref_barcodes),
          names   = head(names(ref_barcodes), 10),
          lengths = as.integer(width(ref_barcodes))
        )

        # ---- STEP 2: Train IDTAXA model -----------------------
        incProgress(0.25, detail = "Training taxonomic classifier…")
        append_log("=== STEP 2: Train IdTaxa model ===")
        append_log("K-mer size:", isolate(input$kmer))

        species_raw   <- gsub("^(\\S+\\s+\\S+).*", "\\1", names(ref_barcodes))
        train_labels  <- paste0("Root;Mitochondria;", species_raw, ";")
        trained_model <- LearnTaxa(ref_barcodes, train_labels,
                                   maxKmerSize = isolate(input$kmer))
        append_log("Model trained on", length(ref_barcodes), "sequences.")

        # ---- STEP 3: Load + classify test sequences -----------
        incProgress(0.25, detail = "Classifying test sequences…")
        append_log("=== STEP 3: Classify test sequences ===")
        req(input$test_files)

        test_list <- lapply(seq_len(nrow(input$test_files)), function(i) {
          append_log("  Reading:", input$test_files$name[i])
          seqs <- smart_read(input$test_files$datapath[i])
          if (is.null(names(seqs)) || all(names(seqs) == "")) {
            names(seqs) <- paste0(tools::file_path_sans_ext(input$test_files$name[i]),
                                  "_", seq_along(seqs))
          }
          seqs
        })
        test_dna <- do.call(c, test_list)
        append_log("Total test sequences:", length(test_dna))

        ids <- IdTaxa(test_dna, trained_model,
                      strand    = "both",
                      threshold = isolate(input$threshold))
        append_log("Classification complete.")

        # ---- STEP 4: Alignment (ref + test) -------------------
        incProgress(0.15, detail = "Aligning all sequences…")
        append_log("=== STEP 4: Multiple sequence alignment ===")
        all_seqs  <- c(ref_barcodes, test_dna)
        alignment <- AlignSeqs(all_seqs, verbose = FALSE)
        append_log("Alignment dimensions:", length(alignment), "sequences ×",
                   unique(width(alignment)), "columns")

        # Coverage stats
        aln_mat    <- as.matrix(alignment)
        cov_per_seq <- apply(aln_mat, 1, function(x) sum(x != "-") / length(x))
        col_occ     <- apply(aln_mat, 2, function(x) mean(x != "-"))

        gc()
        pipeline_state("done")
        append_log("Pipeline finished successfully.")

        list(
          ref_barcodes    = ref_barcodes,
          extraction_info = extraction_info,
          trained_model   = trained_model,
          test_dna        = test_dna,
          ids             = ids,
          alignment       = alignment,
          cov_per_seq     = cov_per_seq,
          col_occ         = col_occ
        )
      })
    }, error = function(e) {
      pipeline_state("error")
      append_log("ERROR:", conditionMessage(e))
      NULL
    })
  })

  # ── Tab: Barcode Extraction ──────────────────────────────────
  output$extraction_stats_ui <- renderUI({
    res <- pipeline_result(); req(res)
    ei  <- res$extraction_info
    div(class = "stat-row",
      div(class = "stat-box", div(class = "stat-val", ei$n), div(class = "stat-lbl", "Barcodes extracted")),
      div(class = "stat-box", div(class = "stat-val", round(mean(ei$lengths))), div(class = "stat-lbl", "Mean length (bp)")),
      div(class = "stat-box", div(class = "stat-val", min(ei$lengths)), div(class = "stat-lbl", "Min length (bp)")),
      div(class = "stat-box", div(class = "stat-val", max(ei$lengths)), div(class = "stat-lbl", "Max length (bp)"))
    )
  })

  output$extraction_log <- renderPrint({
    res <- pipeline_result(); req(res)
    ei  <- res$extraction_info
    cat("── Extracted Barcode Names (first 10) ──\n\n")
    cat(paste(seq_along(ei$names), ei$names, sep = ". "), sep = "\n")
    cat("\n── Length Distribution ──\n\n")
    print(summary(ei$lengths))
  })

  output$dl_barcode_ui <- renderUI({
    req(pipeline_result())
    downloadButton("dl_barcodes", "⬇  Download Barcodes FASTA", class = "dl-btn")
  })

  output$dl_barcodes <- downloadHandler(
    filename = function() paste0("barcodes_", Sys.Date(), ".fasta"),
    content  = function(file) {
      res <- pipeline_result(); req(res)
      writeXStringSet(res$ref_barcodes, filepath = file)
    }
  )

  # ── Tab: Taxonomy Results ────────────────────────────────────
  results_df <- reactive({
    res <- pipeline_result(); req(res)
    final <- sapply(res$ids, function(x) tail(x$taxon, 1))
    conf  <- sapply(res$ids, function(x) tail(x$confidence, 1))
    full  <- sapply(res$ids, function(x) paste(x$taxon, collapse = " > "))
    cscores <- sapply(res$ids, function(x) paste(round(x$confidence, 1), collapse = " | "))

    data.frame(
      Sample        = names(res$test_dna),
      Assigned_Taxon = final,
      Confidence    = round(conf, 1),
      Full_Path     = full,
      All_Scores    = cscores,
      stringsAsFactors = FALSE
    )
  })

  output$taxonomy_stats_ui <- renderUI({
    df  <- results_df()
    n   <- nrow(df)
    hi  <- sum(df$Confidence >= 80, na.rm = TRUE)
    med <- sum(df$Confidence >= isolate(input$threshold) & df$Confidence < 80, na.rm = TRUE)
    lo  <- sum(df$Confidence < isolate(input$threshold), na.rm = TRUE)
    div(class = "stat-row",
      div(class = "stat-box", div(class = "stat-val", n),  div(class = "stat-lbl", "Total queries")),
      div(class = "stat-box", div(class = "stat-val", style = "color:var(--accent);", hi), div(class = "stat-lbl", "High confidence (≥80%)")),
      div(class = "stat-box", div(class = "stat-val", style = "color:var(--accent2);", med), div(class = "stat-lbl", "Moderate confidence")),
      div(class = "stat-box", div(class = "stat-val", style = "color:var(--warn);", lo), div(class = "stat-lbl", "Below threshold"))
    )
  })

  output$results_table <- renderDT({
    df <- results_df()
    datatable(df,
              rownames  = FALSE,
              selection = "none",
              options   = list(
                pageLength = 15,
                scrollX    = TRUE,
                dom        = "frtip",
                columnDefs = list(list(className = "dt-left", targets = "_all"))
              )
    ) %>%
      formatStyle("Confidence",
                  background = styleInterval(
                    c(isolate(input$threshold), 80),
                    c("rgba(240,136,62,0.18)", "rgba(88,166,255,0.18)", "rgba(57,211,83,0.18)")
                  )
      )
  })

  output$dl_results_ui <- renderUI({
    req(pipeline_result())
    tagList(
      downloadButton("dl_csv", "⬇  Download CSV",  class = "dl-btn"),
      downloadButton("dl_aln", "⬇  Download Alignment FASTA", class = "dl-btn")
    )
  })

  output$dl_csv <- downloadHandler(
    filename = function() paste0("taxonomy_results_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(results_df(), file, row.names = FALSE)
  )

  output$dl_aln <- downloadHandler(
    filename = function() paste0("alignment_", Sys.Date(), ".fasta"),
    content  = function(file) {
      res <- pipeline_result(); req(res)
      writeXStringSet(res$alignment, filepath = file)
    }
  )

  # ── Tab: Alignment QC ───────────────────────────────────────
  output$align_summary <- renderPrint({
    res <- pipeline_result(); req(res)
    cat("── Alignment Summary ──────────────────────────────\n\n")
    cat("Sequences aligned  :", length(res$alignment), "\n")
    cat("Alignment length   :", unique(width(res$alignment)), "columns\n\n")
    cat("── Per-sequence coverage ──────────────────────────\n\n")
    cov_df <- data.frame(
      Sequence = names(res$cov_per_seq),
      Coverage_pct = round(res$cov_per_seq * 100, 2)
    )
    print(summary(cov_df$Coverage_pct))
    cat("\n── Column occupancy ───────────────────────────────\n\n")
    cat("Mean col. occupancy :", round(mean(res$col_occ) * 100, 2), "%\n")
    cat("Min  col. occupancy :", round(min(res$col_occ)  * 100, 2), "%\n")
    cat("Max  col. occupancy :", round(max(res$col_occ)  * 100, 2), "%\n")
  })

  output$coverage_plot <- renderPlot({
    res <- pipeline_result(); req(res)
    par(bg = "#f8f9fa", fg = "#54595d", col.axis = "#54595d",
        col.lab = "#202122", col.main = "#202122",
        family = "serif", mar = c(4, 4, 3, 2))
    plot(res$col_occ, type = "l", col = "#3366cc", lwd = 1.8,
         main = "Column Occupancy per Alignment Position",
         ylab = "Proportion of Sequences Present",
         xlab = "Alignment Position",
         ylim = c(0, 1), axes = FALSE)
    axis(1, col = "#a2a9b1", col.ticks = "#a2a9b1")
    axis(2, col = "#a2a9b1", col.ticks = "#a2a9b1", las = 1)
    abline(h = mean(res$col_occ), col = "#c00", lty = 2, lwd = 1)
    legend("bottomright", legend = "Mean occupancy", col = "#c00",
           lty = 2, lwd = 1, bty = "n", text.col = "#54595d", cex = 0.9)
    box(col = "#a2a9b1")
  }, bg = "#f8f9fa")

  # ── Tab: Log ────────────────────────────────────────────────
  output$process_log <- renderText({
    paste(log_msgs(), collapse = "\n")
  })
}
