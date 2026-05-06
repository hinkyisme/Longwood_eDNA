# eDNA Pipeline — Two-Page Shiny App
# server.R

library(shiny)
library(DECIPHER)
library(Biostrings)
library(DT)
library(tools)
library(ape)
library(stringr)
library(dplyr)
library(purrr)

server <- function(input, output, session) {

  # ── Active page state ──────────────────────────────────────────────────────
  active_page <- reactiveVal("builder")   # "builder" | "identify"

  observeEvent(input$go_builder, {
    active_page("builder")
    session$sendCustomMessage("setActiveNav", "builder")
  })

  observeEvent(input$go_identify, {
    active_page("identify")
    session$sendCustomMessage("setActiveNav", "identify")
  })

  # ── Shared helpers ─────────────────────────────────────────────────────────
  smart_read <- function(file_path) {
    ext <- tolower(file_ext(file_path))
    if (ext %in% c("fastq", "fq")) {
      return(as(readQualityScaledDNAStringSet(file_path), "DNAStringSet"))
    } else if (ext == "seq") {
      lines <- readLines(file_path)
      if (any(grepl("^>", lines))) return(readDNAStringSet(file_path))
      return(DNAStringSet(paste(lines, collapse = "")))
    } else {
      return(readDNAStringSet(file_path))
    }
  }

  # ============================================================
  #  PAGE 1 — BARCODE BUILDER  (Barcode_builder.R logic)
  # ============================================================

  output$page_builder <- renderUI({
    req(active_page() == "builder")
    fluidRow(
      column(4,
        div(class = "card",
          h4("1. Input GenBank File"),
          fileInput("gb_file", "GenBank flat file (.gb)",
                    accept = c(".gb", ".gbk", ".genbank")),
          div(class = "status-box",
            "Upload a multi-record GenBank file. Each record will be searched
             for 16S rRNA features and extracted to FASTA."
          )
        ),
        div(class = "card",
          h4("2. Extraction Settings"),
          numericInput("min_len", "Min. sequence length (bp)", value = 200, min = 50, step = 50),
          textInput("out_prefix", "Output file prefix", value = "custom_16S_ref")
        ),
        div(class = "card",
          h4("3. Run"),
          actionButton("run_builder", "Extract 16S Barcodes", class = "run-btn"),
          hr(),
          uiOutput("builder_dl_ui")
        )
      ),
      column(8,
        div(class = "card",
          h4("Extraction Log"),
          verbatimTextOutput("builder_log")
        ),
        div(class = "card",
          h4("Extracted Sequences Preview"),
          DTOutput("builder_table")
        ),
        div(class = "card",
          h4("Alignment Coverage"),
          plotOutput("builder_coverage_plot", height = "260px")
        )
      )
    )
  })

  # Reactive: extract 16S from uploaded GenBank
  builder_result <- eventReactive(input$run_builder, {
    req(input$gb_file)

    log_msgs <- character(0)
    log <- function(...) { log_msgs <<- c(log_msgs, paste0(...)) }

    withProgress(message = "Extracting 16S sequences", value = 0, {

      incProgress(0.1, detail = "Reading GenBank file...")
      gb_lines <- readLines(input$gb_file$datapath, warn = FALSE)

      # reverse complement helper
      rc <- function(seq) {
        seq <- toupper(seq)
        chartr("ACGT", "TGCA", paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
      }

      rec_ends  <- which(gb_lines == "//")
      if (length(rec_ends) == 0) stop("No '//' record separators found in file.")
      log("Records found: ", length(rec_ends))

      seq_list <- list()
      rec_start <- 1
      rec_idx   <- 0

      incProgress(0.2, detail = "Parsing records...")
      for (end_idx in rec_ends) {
        rec_idx   <- rec_idx + 1
        rec_lines <- gb_lines[rec_start:(end_idx - 1)]
        rec_start <- end_idx + 1
        if (length(rec_lines) == 0) next

        # Name
        locus_line     <- grep("^LOCUS",     rec_lines, value = TRUE)
        accession_line <- grep("^ACCESSION", rec_lines, value = TRUE)
        locus_name <- if (length(locus_line) > 0) {
          strsplit(locus_line, "\\s+")[[1]][2]
        } else if (length(accession_line) > 0) {
          strsplit(accession_line, "\\s+")[[1]][2]
        } else {
          paste0("Genome_", rec_idx)
        }

        # ORIGIN
        origin_idx <- grep("^ORIGIN", rec_lines)
        if (length(origin_idx) == 0) { log("Record ", rec_idx, " — no ORIGIN, skipped"); next }
        dna_lines <- rec_lines[(origin_idx + 1):length(rec_lines)]
        dna_seq   <- toupper(gsub("[0-9\\s]", "", paste(dna_lines, collapse = "")))
        if (nchar(dna_seq) == 0) next

        # FEATURES
        feat_idx <- grep("^FEATURES", rec_lines)
        if (length(feat_idx) == 0) { log("Record ", rec_idx, " — no FEATURES, skipped"); next }
        feat_lines  <- rec_lines[(feat_idx + 1):(origin_idx - 1)]
        entry_start <- which(grepl("^\\s{5}\\S", feat_lines, perl = TRUE))
        if (length(entry_start) == 0) next

        entry_ranges   <- c(entry_start, length(feat_lines) + 1)
        feature_entries <- lapply(seq_along(entry_start), function(k) {
          paste(feat_lines[entry_start[k]:(entry_ranges[k + 1] - 1)], collapse = " ")
        })

        # 16S hits
        hits <- which(sapply(feature_entries, function(x)
          grepl("16S|16S ribosomal|16S rRNA|16S_RNA|16s", x, ignore.case = TRUE)))
        if (length(hits) == 0)
          hits <- which(sapply(feature_entries, function(x)
            grepl("rRNA", x, ignore.case = TRUE) && grepl("16S|16s", x, ignore.case = TRUE)))
        if (length(hits) == 0) { log("Record ", rec_idx, " — no 16S features"); next }

        seq_count_in_rec <- 0
        for (h in hits) {
          entry  <- feature_entries[[h]]
          coords <- regmatches(entry, gregexpr("[0-9]+\\.{2}[0-9]+", entry, perl = TRUE))[[1]]
          if (length(coords) == 0) next
          is_complement <- grepl("complement\\s*\\(", entry, ignore.case = TRUE)
          parts <- sapply(coords, function(coord) {
            se <- as.numeric(unlist(strsplit(coord, "\\.\\.")))
            if (length(se) != 2) return("")
            se[1] <- max(1, se[1]); se[2] <- min(nchar(dna_seq), se[2])
            substr(dna_seq, se[1], se[2])
          }, USE.NAMES = FALSE)
          seq_string <- gsub("\\s+", "", paste(parts, collapse = ""))
          if (nchar(seq_string) < input$min_len) next
          if (is_complement) seq_string <- rc(seq_string)
          seq_count_in_rec <- seq_count_in_rec + 1
          nm <- if (seq_count_in_rec == 1) locus_name else paste0(locus_name, "_copy", seq_count_in_rec)
          while (nm %in% names(seq_list)) nm <- paste0(nm, "_dup")
          seq_list[[nm]] <- seq_string
        }
      }

      if (length(seq_list) == 0) stop("No 16S sequences passed filters.")
      log("Extracted sequences: ", length(seq_list))

      incProgress(0.4, detail = "Writing FASTA...")
      out_fasta <- tempfile(fileext = ".fasta")
      con <- file(out_fasta, "w")
      for (nm in names(seq_list)) {
        writeLines(paste0(">", nm), con = con)
        writeLines(strwrap(seq_list[[nm]], width = 80), con = con)
      }
      close(con)

      incProgress(0.2, detail = "Aligning for QC...")
      dna_set   <- DNAStringSet(seq_list)
      aln       <- tryCatch(AlignSeqs(dna_set, verbose = FALSE), error = function(e) NULL)
      cov_df    <- NULL
      col_occ   <- NULL
      if (!is.null(aln)) {
        mat      <- as.matrix(aln)
        coverage <- apply(mat, 1, function(x) sum(x != "-") / length(x)) * 100
        cov_df   <- data.frame(Name     = names(aln),
                               Length   = nchar(as.character(aln)),
                               Coverage = round(coverage, 2))
        col_occ  <- apply(mat, 2, function(x) mean(x != "-"))
        log("Mean column occupancy: ", round(mean(col_occ) * 100, 1), "%")
      }

      incProgress(0.1, detail = "Done.")
      list(seq_list = seq_list, out_fasta = out_fasta,
           cov_df = cov_df, col_occ = col_occ, log = log_msgs)
    })
  })

  output$builder_log <- renderText({
    req(builder_result())
    paste(builder_result()$log, collapse = "\n")
  })

  output$builder_table <- renderDT({
    req(builder_result())
    df <- data.frame(
      Name   = names(builder_result()$seq_list),
      Length = nchar(unlist(builder_result()$seq_list)),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$builder_coverage_plot <- renderPlot({
    req(builder_result(), !is.null(builder_result()$col_occ))
    occ <- builder_result()$col_occ
    par(mar = c(4, 4, 2, 1), bg = "#fff", col.axis = "#1a2e44", col.lab = "#1a2e44")
    plot(occ, type = "l", col = "#7ecfb3", lwd = 2,
         main = "Column Occupancy per Alignment Position",
         ylab = "Proportion Present", xlab = "Position",
         ylim = c(0, 1), frame.plot = FALSE)
    abline(h = mean(occ), col = "#1a2e44", lty = 2, lwd = 1.2)
    legend("bottomright", legend = paste0("Mean: ", round(mean(occ)*100, 1), "%"),
           lty = 2, col = "#1a2e44", bty = "n", cex = .85)
  })

  output$builder_dl_ui <- renderUI({
    req(builder_result())
    downloadButton("dl_fasta", "Download FASTA", class = "btn btn-success dl-btn")
  })

  output$dl_fasta <- downloadHandler(
    filename = function() paste0(input$out_prefix, "_16S.fasta"),
    content  = function(file) file.copy(builder_result()$out_fasta, file)
  )

  # ============================================================
  #  PAGE 2 — IDENTIFY TAXA  (IdentifyTaxa.R logic)
  # ============================================================

  output$page_identify <- renderUI({
    req(active_page() == "identify")
    fluidRow(
      column(4,
        div(class = "card",
          h4("1. Reference Barcodes"),
          fileInput("ref_fasta", "Reference FASTA (16S barcodes)",
                    accept = c(".fasta", ".fa", ".seq", ".fastq", ".fq")),
          div(class = "status-box",
            "Use the FASTA exported from the Barcode Builder, or any
             pre-built 16S reference library."
          )
        ),
        div(class = "card",
          h4("2. Test Sequences"),
          fileInput("test_files", "Query sequences (multiple OK)",
                    multiple = TRUE,
                    accept   = c(".fasta", ".fa", ".seq", ".fastq", ".fq"))
        ),
        div(class = "card",
          h4("3. Primer Trimming"),
          textInput("fwd", "Forward Primer (5'→3')", "GGGTCACCAACTCCGCTAAC"),
          textInput("rev", "Reverse Primer  (5'→3')", "GGGCAGCTAAGGCTGGAAAA"),
          sliderInput("mismatch", "Max Mismatches", min = 0, max = 5, value = 2)
        ),
        div(class = "card",
          h4("4. Classifier Settings"),
          numericInput("threshold", "Confidence Threshold (%)",
                       value = 50, min = 0, max = 100),
          actionButton("run_identify", "Identify Taxa", class = "run-btn"),
          hr(),
          uiOutput("identify_dl_ui")
        )
      ),
      column(8,
        div(class = "card",
          h4("Taxonomy Results"),
          DTOutput("results_table")
        ),
        div(class = "card",
          h4("Summary"),
          verbatimTextOutput("identify_summary")
        ),
        div(class = "card",
          h4("Alignment Browser"),
          verbatimTextOutput("align_summary")
        )
      )
    )
  })

  identify_result <- eventReactive(input$run_identify, {
    req(input$ref_fasta, input$test_files)

    withProgress(message = "Running Taxonomic Pipeline", value = 0, {

      # 1. Load reference
      incProgress(0.15, detail = "Loading reference barcodes...")
      ref_dna <- smart_read(input$ref_fasta$datapath)

      # 2. Trim primers
      incProgress(0.15, detail = "Trimming primers from reference...")
      extracted_ref <- TrimDNA(ref_dna,
                               left          = input$fwd,
                               right         = input$rev,
                               type          = "sequences",
                               max.mismatch  = input$mismatch)
      # Drop zero-length sequences after trimming
      keep <- width(extracted_ref) > 0
      extracted_ref <- extracted_ref[keep]

      # 3. Train
      incProgress(0.25, detail = "Training classifier (LearnTaxa)...")
      species_names <- gsub("^(\\w+\\s+\\w+).*", "\\1", names(extracted_ref))
      train_labels  <- paste0("Root;Mitochondria;", species_names, ";")
      trained_model <- LearnTaxa(extracted_ref, train_labels, maxKmerSize = 7)

      # 4. Load test sequences (all files)
      incProgress(0.15, detail = "Loading test sequences...")
      test_list <- lapply(seq_len(nrow(input$test_files)), function(i) {
        tryCatch(smart_read(input$test_files$datapath[i]),
                 error = function(e) NULL)
      })
      test_list <- Filter(Negate(is.null), test_list)
      test_dna  <- do.call(c, test_list)
      # Use original filenames (without extension) as sequence names
      clean_names <- gsub("\\.[^.]+$", "", input$test_files$name)
      if (length(clean_names) == length(test_dna)) names(test_dna) <- clean_names

      # 5. Classify
      incProgress(0.2, detail = "Classifying with IdTaxa...")
      ids <- IdTaxa(test_dna, trained_model, strand = "both",
                    threshold = input$threshold)

      # 6. Align reference + test
      incProgress(0.05, detail = "Aligning sequences...")
      all_seqs  <- c(extracted_ref, test_dna)
      alignment <- tryCatch(AlignSeqs(all_seqs, verbose = FALSE), error = function(e) NULL)

      incProgress(0.05, detail = "Finalising...")
      gc()

      list(ids = ids, test_dna = test_dna, alignment = alignment)
    })
  })

  output$results_table <- renderDT({
    req(identify_result())
    data <- identify_result()
    final_taxon <- sapply(data$ids, function(x) tail(x$taxon,      n = 1))
    final_conf  <- sapply(data$ids, function(x) tail(x$confidence, n = 1))
    full_path   <- sapply(data$ids, function(x) paste(x$taxon, collapse = " > "))
    df <- data.frame(
      Sample_Name     = names(data$test_dna),
      Assigned_Taxon  = final_taxon,
      Confidence      = round(final_conf, 1),
      Full_Path       = full_path,
      stringsAsFactors = FALSE
    )
    df <- df[order(df$Confidence, decreasing = TRUE), ]
    datatable(df, options = list(pageLength = 10, scrollX = TRUE)) |>
      formatStyle("Confidence",
                  background = styleColorBar(c(0, 100), "#7ecfb3"),
                  backgroundSize = "100% 80%",
                  backgroundRepeat = "no-repeat",
                  backgroundPosition = "center")
  })

  output$identify_summary <- renderText({
    req(identify_result())
    data   <- identify_result()
    confs  <- sapply(data$ids, function(x) tail(x$confidence, n = 1))
    taxa   <- sapply(data$ids, function(x) tail(x$taxon, n = 1))
    passed <- sum(confs >= input$threshold)
    paste0(
      "Sequences classified:  ", length(data$ids), "\n",
      "Above threshold (>=", input$threshold, "%):  ", passed, "\n",
      "Unique taxa assigned:  ", length(unique(taxa)), "\n",
      "Median confidence:     ", round(median(confs), 1), "%\n"
    )
  })

  output$align_summary <- renderText({
    req(identify_result())
    aln <- identify_result()$alignment
    if (is.null(aln)) return("Alignment could not be computed.")
    mat <- as.matrix(aln)
    occ <- mean(apply(mat, 2, function(x) mean(x != "-"))) * 100
    paste0(
      "Sequences aligned:   ", nrow(mat), "\n",
      "Alignment length:    ", ncol(mat), " bp\n",
      "Mean col. occupancy: ", round(occ, 1), "%\n\n",
      "(Use BrowseSeqs() locally for interactive viewing)"
    )
  })

  # Download CSV
  output$identify_dl_ui <- renderUI({
    req(identify_result())
    downloadButton("dl_csv", "Download CSV", class = "btn btn-success dl-btn")
  })

  output$dl_csv <- downloadHandler(
    filename = function() paste0("barcode_results_", Sys.Date(), ".csv"),
    content  = function(file) {
      data <- identify_result()
      df <- data.frame(
        Sample_Name    = names(data$test_dna),
        Assigned_Taxon = sapply(data$ids, function(x) tail(x$taxon, n = 1)),
        Confidence     = sapply(data$ids, function(x) round(tail(x$confidence, n = 1), 1)),
        Full_Path      = sapply(data$ids, function(x) paste(x$taxon, collapse = " > ")),
        stringsAsFactors = FALSE
      )
      write.csv(df[order(df$Confidence, decreasing = TRUE), ], file, row.names = FALSE)
    }
  )
}
