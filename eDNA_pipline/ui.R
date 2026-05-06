# ============================================================
# ui.R  —  eDNA Barcode Builder & Identifier
# ============================================================

library(shiny)
library(DT)
library(tools)
library(shinycssloaders)

options(shiny.maxRequestSize = 500 * 1024^2)

ui <- fluidPage(

  # ── HEAD: fonts + global CSS ──────────────────────────────
  tags$head(
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Linux+Libertine+O&family=Source+Serif+4:wght@400;600&family=Source+Code+Pro:wght@400;600&display=swap",
      rel  = "stylesheet"
    ),
    tags$style(HTML("

      /* ── Wikipedia-inspired Palette ─────────────────────── */
      :root {
        --bg:        #f8f9fa;
        --surface:   #ffffff;
        --border:    #a2a9b1;
        --border-lt: #eaecf0;
        --accent:    #3366cc;
        --accent2:   #0645ad;
        --warn:      #c00;
        --text:      #202122;
        --muted:     #54595d;
        --subtle:    #72777d;
        --highlight: #eaf3fb;
        --font-head: 'Source Serif 4', 'Linux Libertine O', Georgia, 'Times New Roman', serif;
        --font-body: -apple-system, 'Linux Libertine O', Georgia, serif;
        --font-mono: 'Source Code Pro', 'Courier New', monospace;
        --radius:    2px;
      }

      /* ── Reset ──────────────────────────────────────────── */
      *, *::before, *::after { box-sizing: border-box; }
      html, body { height: 100%; }
      body {
        background: var(--bg);
        color: var(--text);
        font-family: var(--font-body);
        font-size: 14px;
        margin: 0;
      }

      /* ── Scrollbar ─────────────────────────────────────── */
      ::-webkit-scrollbar { width: 8px; height: 8px; }
      ::-webkit-scrollbar-track { background: var(--bg); }
      ::-webkit-scrollbar-thumb { background: var(--border); border-radius: 2px; }

      /* ── Banner ─────────────────────────────────────────── */
      .app-banner {
        background: var(--surface);
        border-bottom: 1px solid var(--border);
        padding: 14px 24px 12px;
        display: flex;
        align-items: center;
        gap: 14px;
      }
      .app-banner .logo-box {
        width: 42px; height: 42px;
        border: 1px solid var(--border-lt);
        border-radius: 50%;
        display: flex; align-items: center; justify-content: center;
        flex-shrink: 0;
        background: var(--bg);
      }
      .app-banner .logo-box svg { width: 26px; height: 26px; fill: var(--accent); }
      .app-banner h1 {
        font-family: var(--font-head);
        font-size: 20px;
        font-weight: 600;
        color: var(--text);
        margin: 0 0 2px;
        letter-spacing: 0;
        border-bottom: none;
      }
      .app-banner p {
        margin: 0;
        color: var(--subtle);
        font-size: 12px;
        font-style: italic;
        letter-spacing: 0;
        text-transform: none;
      }

      /* ── Layout ─────────────────────────────────────────── */
      .main-grid {
        display: grid;
        grid-template-columns: 310px 1fr;
        height: calc(100vh - 74px);
        overflow: hidden;
      }
      .left-panel {
        background: var(--bg);
        border-right: 1px solid var(--border-lt);
        overflow-y: auto;
        padding: 16px 14px;
      }
      .right-panel {
        background: var(--surface);
        overflow-y: auto;
        padding: 20px 28px;
      }

      /* ── Section cards ───────────────────────────────────── */
      .card {
        background: var(--surface);
        border: 1px solid var(--border-lt);
        border-radius: var(--radius);
        padding: 12px 14px;
        margin-bottom: 12px;
      }
      .card-title {
        font-family: var(--font-head);
        font-size: 13px;
        font-weight: 600;
        text-transform: none;
        letter-spacing: 0;
        color: var(--text);
        margin: 0 0 10px;
        padding-bottom: 6px;
        border-bottom: 1px solid var(--border-lt);
        display: flex;
        align-items: center;
        gap: 8px;
      }
      .card-title .step-badge {
        background: var(--accent);
        color: #fff;
        border-radius: 2px;
        padding: 1px 6px;
        font-size: 10px;
        font-weight: 600;
        font-family: var(--font-mono);
        line-height: 16px;
      }

      /* ── Inputs ──────────────────────────────────────────── */
      label { color: var(--muted); font-size: 12px; }
      .form-control, .selectize-input {
        background: var(--surface) !important;
        border: 1px solid var(--border) !important;
        color: var(--text) !important;
        border-radius: var(--radius) !important;
        font-family: var(--font-body) !important;
        font-size: 13px !important;
      }
      .form-control:focus, .selectize-input.focus {
        border-color: var(--accent) !important;
        box-shadow: 0 0 0 2px rgba(51,102,204,0.15) !important;
        outline: none !important;
      }
      input[type='file'] { color: var(--muted); }
      .btn-file {
        background: var(--bg) !important;
        border: 1px solid var(--border) !important;
        color: var(--text) !important;
        font-family: var(--font-body) !important;
      }
      .irs--shiny .irs-bar, .irs--shiny .irs-handle { background: var(--accent); }
      .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
        background: var(--accent); color: #fff;
      }

      /* ── Run button ──────────────────────────────────────── */
      #run {
        width: 100%;
        background: var(--accent) !important;
        border: 1px solid #2a4b8d !important;
        border-radius: var(--radius) !important;
        color: #fff !important;
        font-family: var(--font-body) !important;
        font-weight: 600 !important;
        font-size: 13px !important;
        letter-spacing: 0 !important;
        text-transform: none !important;
        padding: 10px !important;
        margin-top: 6px !important;
        cursor: pointer;
        transition: background 0.15s ease;
      }
      #run:hover { background: #2a4b8d !important; }
      #run:active { background: #1e3a6e !important; }

      /* ── Download buttons ────────────────────────────────── */
      .dl-btn {
        background: var(--bg) !important;
        border: 1px solid var(--border) !important;
        color: var(--accent2) !important;
        font-family: var(--font-body) !important;
        font-size: 12px !important;
        border-radius: var(--radius) !important;
        padding: 5px 10px !important;
        margin-right: 6px;
        margin-bottom: 6px;
        cursor: pointer;
        transition: background 0.1s;
        text-decoration: none;
      }
      .dl-btn:hover { background: var(--highlight) !important; color: var(--accent2) !important; }

      /* ── Tab bar ─────────────────────────────────────────── */
      .nav-tabs {
        border-bottom: 1px solid var(--border) !important;
        margin-bottom: 18px;
      }
      .nav-tabs > li > a {
        font-family: var(--font-body) !important;
        font-size: 13px !important;
        font-weight: 400 !important;
        text-transform: none !important;
        letter-spacing: 0 !important;
        color: var(--accent2) !important;
        background: var(--bg) !important;
        border: 1px solid transparent !important;
        border-bottom: none !important;
        padding: 7px 14px !important;
        border-radius: 2px 2px 0 0 !important;
        transition: background 0.1s;
      }
      .nav-tabs > li.active > a {
        color: var(--text) !important;
        background: var(--surface) !important;
        border-color: var(--border) var(--border) var(--surface) !important;
        font-weight: 600 !important;
      }
      .nav-tabs > li > a:hover {
        background: var(--highlight) !important;
        color: var(--accent2) !important;
        border-color: var(--border-lt) !important;
      }

      /* ── Status badge ────────────────────────────────────── */
      .status-bar {
        font-size: 12px;
        color: var(--muted);
        margin-bottom: 14px;
        display: flex;
        align-items: center;
        gap: 8px;
        font-style: italic;
      }
      .status-dot {
        width: 7px; height: 7px;
        border-radius: 50%;
        background: var(--border);
        flex-shrink: 0;
      }
      .status-dot.ready { background: #14866d; box-shadow: 0 0 5px #14866d55; }
      .status-dot.running { background: var(--warn); animation: pulse 1s ease-in-out infinite; }
      @keyframes pulse {
        0%, 100% { opacity: 1; }
        50% { opacity: 0.3; }
      }

      /* ── Verbatim output ─────────────────────────────────── */
      pre {
        background: var(--bg) !important;
        border: 1px solid var(--border-lt) !important;
        color: var(--text) !important;
        font-family: var(--font-mono) !important;
        font-size: 12px !important;
        border-radius: var(--radius) !important;
        padding: 14px !important;
        white-space: pre-wrap !important;
        word-break: break-all;
      }

      /* ── DT table overrides ──────────────────────────────── */
      .dataTables_wrapper { color: var(--text) !important; }
      table.dataTable thead th {
        background: var(--bg) !important;
        color: var(--text) !important;
        font-family: var(--font-body) !important;
        font-size: 13px !important;
        font-weight: 600 !important;
        text-transform: none !important;
        letter-spacing: 0 !important;
        border-bottom: 2px solid var(--border) !important;
      }
      table.dataTable tbody tr { background: var(--surface) !important; }
      table.dataTable tbody tr:hover td { background: var(--highlight) !important; }
      table.dataTable tbody tr:nth-child(even) td { background: var(--bg) !important; }
      table.dataTable tbody tr:nth-child(even):hover td { background: var(--highlight) !important; }
      table.dataTable tbody td { border-top: 1px solid var(--border-lt) !important; color: var(--text) !important; font-size: 13px !important; }
      .dataTables_info, .dataTables_length, .dataTables_filter,
      .dataTables_paginate { color: var(--muted) !important; font-size: 12px !important; }
      .dataTables_paginate .paginate_button { color: var(--accent2) !important; }
      .dataTables_paginate .paginate_button.current {
        background: var(--bg) !important;
        border: 1px solid var(--border) !important;
        color: var(--text) !important;
        border-radius: 2px !important;
      }
      .dataTables_filter input {
        background: var(--surface) !important;
        border: 1px solid var(--border) !important;
        color: var(--text) !important;
        border-radius: var(--radius) !important;
        padding: 4px 8px !important;
      }

      /* ── Summary stat boxes ──────────────────────────────── */
      .stat-row { display: flex; gap: 10px; margin-bottom: 16px; flex-wrap: wrap; }
      .stat-box {
        flex: 1 1 100px;
        background: var(--bg);
        border: 1px solid var(--border-lt);
        border-radius: var(--radius);
        padding: 10px 14px;
      }
      .stat-box .stat-val {
        font-family: var(--font-head);
        font-size: 24px;
        font-weight: 600;
        color: var(--accent);
        line-height: 1;
        margin-bottom: 4px;
      }
      .stat-box .stat-lbl {
        font-size: 11px;
        color: var(--subtle);
        text-transform: none;
        letter-spacing: 0;
      }

      /* ── Primer tag ──────────────────────────────────────── */
      .primer-tag {
        display: inline-block;
        background: var(--bg);
        border: 1px solid var(--border-lt);
        border-radius: 2px;
        padding: 2px 7px;
        font-size: 11px;
        font-family: var(--font-mono);
        color: var(--text);
        word-break: break-all;
        margin-top: 2px;
      }

      /* ── Divider ─────────────────────────────────────────── */
      hr { border: none; border-top: 1px solid var(--border-lt); margin: 12px 0; }

      /* ── Mode toggle ─────────────────────────────────────── */
      .mode-toggle {
        display: flex;
        background: var(--bg);
        border: 1px solid var(--border);
        border-radius: var(--radius);
        overflow: hidden;
        margin-bottom: 12px;
      }
      .mode-toggle .mode-btn {
        flex: 1;
        text-align: center;
        padding: 7px 4px;
        font-family: var(--font-body);
        font-size: 12px;
        font-weight: 400;
        text-transform: none;
        letter-spacing: 0;
        color: var(--accent2);
        cursor: pointer;
        transition: background 0.1s, color 0.1s;
        user-select: none;
      }
      .mode-toggle .mode-btn.active {
        background: var(--accent);
        color: #fff;
        font-weight: 600;
      }

      /* misc */
      .shiny-progress-container { z-index: 9999; }
      .shiny-progress .progress-bar { background-color: var(--accent) !important; }
      .file-badge {
        display: inline-block;
        background: var(--highlight);
        border: 1px solid #c8ccd1;
        border-radius: 2px;
        color: var(--accent2);
        padding: 1px 6px;
        font-size: 10px;
        margin-left: 6px;
        vertical-align: middle;
      }
      small { color: var(--subtle); font-style: italic; }
    "))
  ),

  # ── Banner ─────────────────────────────────────────────────
  div(class = "app-banner",
    div(class = "logo-box",
      tags$svg(xmlns = "http://www.w3.org/2000/svg", viewBox = "0 0 24 24",
        tags$path(d = "M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-1 14H9V8h2v8zm4 0h-2V8h2v8z")
      )
    ),
    div(
      tags$h1("eDNA Barcode Builder"),
      tags$p("Non-model organism • Mitochondrial & Genomic DNA • 16S / COI / Custom loci")
    )
  ),

  # ── Main grid ──────────────────────────────────────────────
  div(class = "main-grid",

    # ── LEFT PANEL ──────────────────────────────────────────
    div(class = "left-panel",

      # Step 1 — Input mode
      div(class = "card",
        div(class = "card-title",
          span(class = "step-badge", "01"),
          "Input Mode"
        ),
        div(class = "mode-toggle",
          div(id = "mode_genbank_btn", class = "mode-btn active",
              onclick = "Shiny.setInputValue('input_mode','genbank'); document.getElementById('mode_genbank_btn').classList.add('active'); document.getElementById('mode_fasta_btn').classList.remove('active'); document.getElementById('panel_genbank').style.display='block'; document.getElementById('panel_fasta').style.display='none';",
              "GenBank (.gb)"
          ),
          div(id = "mode_fasta_btn", class = "mode-btn",
              onclick = "Shiny.setInputValue('input_mode','fasta'); document.getElementById('mode_fasta_btn').classList.add('active'); document.getElementById('mode_genbank_btn').classList.remove('active'); document.getElementById('panel_genbank').style.display='none'; document.getElementById('panel_fasta').style.display='block';",
              "FASTA / FASTQ"
          )
        ),

        # GenBank panel
        div(id = "panel_genbank",
          fileInput("gb_file", "GenBank Flat File (.gb)",
                    accept = c(".gb", ".gbk", ".genbank")),
          tags$small(style = "color:var(--muted);",
            "Multi-record files supported. 16S rRNA features are auto-extracted.")
        ),

        # FASTA panel (hidden initially)
        div(id = "panel_fasta", style = "display:none;",
          fileInput("ref_fasta", "Reference Genomes (FASTA / FASTQ / .seq)",
                    accept = c(".fasta", ".fa", ".seq", ".fastq", ".fq")),
          tags$small(style = "color:var(--muted);",
            "Full mitogenomes or pre-extracted amplicons.")
        )
      ),

      # Step 2 — PCR Primers
      div(class = "card",
        div(class = "card-title",
          span(class = "step-badge", "02"),
          "PCR Primer Trimming"
        ),
        textInput("fwd", "Forward Primer (5'→3')", "GGGTCACCAACTCCGCTAAC",
                  placeholder = "e.g. GGTCAACAAATCATAAAGATATTGG"),
        textInput("rev", "Reverse Primer (5'→3')", "GGGCAGCTAAGGCTGGAAAA",
                  placeholder = "e.g. TAAACTTCAGGGTGACCAAAAAATCA"),
        sliderInput("mismatch", "Max Primer Mismatches", min = 0, max = 5, value = 2, ticks = FALSE),
        tags$small(style = "color:var(--muted);",
          "Default primers target 16S in Mytilus / mussel mitogenomes.")
      ),

      # Step 3 — Classifier
      div(class = "card",
        div(class = "card-title",
          span(class = "step-badge", "03"),
          "Classifier Settings"
        ),
        numericInput("threshold", "Confidence Threshold (%)", value = 50, min = 0, max = 100),
        numericInput("kmer", "K-mer Size", value = 8, min = 4, max = 12),
        tags$small(style = "color:var(--muted);", "Lower K-mer = less memory. Recommended: 8.")
      ),

      # Step 4 — Test sequences
      div(class = "card",
        div(class = "card-title",
          span(class = "step-badge", "04"),
          "Test Sequences"
        ),
        fileInput("test_files", NULL,
                  multiple = TRUE,
                  accept = c(".fasta", ".fa", ".seq", ".fastq", ".fq"),
                  placeholder = "FASTA / FASTQ / .seq"),
        tags$small(style = "color:var(--muted);",
          "Sanger reads or short amplicons to identify.")
      ),

      # Run button
      actionButton("run", "▶  Execute Pipeline")
    ),

    # ── RIGHT PANEL ─────────────────────────────────────────
    div(class = "right-panel",

      # Status bar
      uiOutput("status_ui"),

      tabsetPanel(id = "main_tabs",

        # ── Tab A: Extraction QC ─────────────────────────────
        tabPanel("Barcode Extraction",
          br(),
          uiOutput("extraction_stats_ui"),
          withSpinner(verbatimTextOutput("extraction_log"), type = 4, color = "#39d353"),
          br(),
          uiOutput("dl_barcode_ui")
        ),

        # ── Tab B: Taxonomy Results ──────────────────────────
        tabPanel("Taxonomy Results",
          br(),
          uiOutput("taxonomy_stats_ui"),
          withSpinner(DTOutput("results_table"), type = 4, color = "#39d353"),
          br(),
          uiOutput("dl_results_ui")
        ),

        # ── Tab C: Alignment Summary ─────────────────────────
        tabPanel("Alignment QC",
          br(),
          withSpinner(verbatimTextOutput("align_summary"), type = 4, color = "#39d353"),
          plotOutput("coverage_plot", height = "260px")
        ),

        # ── Tab D: Log ───────────────────────────────────────
        tabPanel("Log",
          br(),
          verbatimTextOutput("process_log")
        )
      )
    ) # right-panel
  ) # main-grid
)
