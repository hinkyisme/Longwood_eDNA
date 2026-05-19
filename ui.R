#
# ui.R  –  eDNA Wiki - Build Barcodes, Primers, and Identify Taxa.
#

library(shiny)
library(DECIPHER)
library(Biostrings)
library(DT)
library(tools)
library(httr)
library(rvest)

# Increase upload limit to 500 MB
options(shiny.maxRequestSize = 500 * 1024^2)

# ── CSS for sequence viewer and general polish ──────────────────────────────
seq_viewer_css <- "
  .seq-viewer {
    background: #1e1e2e;
    border-radius: 6px;
    padding: 14px 16px;
    margin-top: 10px;
    overflow-x: auto;
  }
  .seq-pre {
    font-family: 'Courier New', Courier, monospace;
    font-size: 13px;
    line-height: 1.7;
    color: #cdd6f4;
    margin: 0;
    white-space: pre;
  }
  .pos-num {
    color: #6c7086;
    user-select: none;
    font-size: 11px;
  }
  .base.fwd-primer {
    background-color: rgba(137,220,235,0.25);
    border-bottom: 2px solid #89dceb;
    font-weight: bold;
  }
  .base.rev-primer {
    background-color: rgba(250,179,135,0.25);
    border-bottom: 2px solid #fab387;
    font-weight: bold;
  }
  .seq-legend {
    margin-bottom: 8px;
    font-size: 12px;
    color: #a6adc8;
  }
  .seq-legend .leg { margin-right: 18px; }
  .fwd-leg { color: #89dceb; }
  .rev-leg { color: #fab387; }
  /* Primer score badge */
  .score-good  { color: #2e7d32; font-weight: 600; }
  .score-ok    { color: #e65100; font-weight: 600; }
  .score-poor  { color: #c62828; font-weight: 600; }
"

# ── UI definition ────────────────────────────────────────────────────────────
ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  tags$head(tags$style(HTML(seq_viewer_css))),
  titlePanel("eDNA Wiki - Build Barcodes, Primers, and Identify Taxa."),

  sidebarLayout(

    # ── SIDEBAR ──────────────────────────────────────────────
    sidebarPanel(
      width = 3,

      # ---------- Taxonomy Pipeline controls ----------
      conditionalPanel(
        condition = "input.main_tabs == 'Taxonomy Table' ||
                     input.main_tabs == 'Alignment Summary' ||
                     input.main_tabs == 'Log'",

        h4("1. Upload Genomic Data"),
        fileInput("ref_fasta", "Reference Genomes (FASTA / FASTQ / .seq)",
                  accept = c(".fasta",".fa",".seq",".fastq",".fq")),
        fileInput("test_files", "Test Datasets (multiple files allowed)",
                  multiple = TRUE,
                  accept   = c(".fasta",".fa",".seq",".fastq",".fq")),

        hr(),
        h4("2. PCR Primer Trimming"),
        textInput("fwd", "Forward Primer (5'→3')", "GGTCAACAAATCATAAAGATATTGG"),
        textInput("rev", "Reverse Primer (5'→3')", "TAAACTTCAGGGTGACCAAAAAATCA"),
        sliderInput("mismatch", "Max Primer Mismatch", min=0, max=5, value=2),

        hr(),
        h4("3. Classifier Settings"),
        numericInput("threshold", "Confidence Threshold (%)", value=50, min=0, max=100),

        actionButton("run", "Execute Pipeline",
                     class="btn-primary btn-lg", width="100%")
      ),

      # ---------- Barcode Builder controls ----------
      conditionalPanel(
        condition = "input.main_tabs == 'Barcode Builder'",

        h4("1. Upload GenBank File"),
        fileInput("gb_file", "GenBank Flat File (.gb / .gbk)",
                  accept = c(".gb",".gbk",".genbank")),

        hr(),
        actionButton("run_barcode", "Extract 16S Barcodes",
                     class="btn-success btn-lg", width="100%")
      ),

      # ---------- Primer Designer controls ----------
      conditionalPanel(
        condition = "input.main_tabs == 'Primer Designer'",

        h4("1. Target Sequence"),
        tags$p(tags$small(tags$em(
          "Run the Barcode Builder first to auto-populate the consensus, ",
          "or paste a sequence below."
        ))),
        textAreaInput("custom_barcode",
                      "Barcode Sequence (optional override)",
                      placeholder = "Paste a FASTA or raw DNA sequence here…",
                      rows = 4, resize = "vertical"),

        hr(),
        h4("2. Primer Length"),
        fluidRow(
          column(6, numericInput("primer_kmin", "Min (bp)", value=18, min=15, max=30)),
          column(6, numericInput("primer_kmax", "Max (bp)", value=24, min=15, max=35))
        ),

        hr(),
        h4("3. Amplicon Size"),
        fluidRow(
          column(6, numericInput("amp_min", "Min (bp)", value=100, min=50,  max=2000)),
          column(6, numericInput("amp_max", "Max (bp)", value=600, min=100, max=3000))
        ),

        hr(),
        h4("4. Tm Matching"),
        sliderInput("tm_delta", "Max Fwd/Rev Tm Difference (°C)",
                    min=1, max=10, value=5, step=0.5),

        hr(),
        actionButton("run_primers", "Design Primers",
                     class="btn-warning btn-lg", width="100%"),

        hr(),
        h4("5. Primer-BLAST Specificity Check"),
        tags$p(tags$small(tags$em(
          "Select a pair in the table, configure below, then submit to NCBI."
        ))),
        textInput("blast_organism",
                  "Organism / Taxon Name",
                  value = "Mytilus",
                  placeholder = "e.g. Mytilus, Bivalvia, all"),
        selectInput("blast_db", "BLAST Database",
                    choices = c(
                      "Nucleotide collection (nt)" = "nt",
                      "RefSeq RNA"                 = "refseq_rna",
                      "RefSeq Genomic"             = "refseq_genomic",
                      "Whole-genome shotgun (wgs)" = "wgs"
                    ),
                    selected = "nt"),
        actionButton("run_primerblast",
                     "Check Selected Pair via Primer-BLAST",
                     class = "btn-info btn-lg", width = "100%"),

        br(), br(),
        uiOutput("download_primers_btn")
      )
    ),

    # ── MAIN PANEL ───────────────────────────────────────────
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",

        # ---- Taxonomy results ----
        tabPanel("Taxonomy Table",
          br(),
          DTOutput("results_table")
        ),

        tabPanel("Alignment Summary",
          br(),
          verbatimTextOutput("align_summary")
        ),

        tabPanel("Log",
          br(),
          verbatimTextOutput("process_log")
        ),

        # ---- Barcode Builder ----
        tabPanel("Barcode Builder",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h5("Extraction Summary"),
                uiOutput("barcode_summary"),
                br(),
                downloadButton("download_barcodes", "Download FASTA",
                               class="btn-sm btn-default")
              )
            ),
            column(8,
              h5("Extracted 16S Sequences"),
              DTOutput("barcode_table")
            )
          )
        ),

        # ---- Primer Designer ----
        tabPanel("Primer Designer",
          br(),
          fluidRow(
            column(12,
              h5("Ranked Primer Pairs"),
              tags$p(tags$small(
                "Pairs are scored by a combined penalty: Tm deviation from 60 °C, ",
                "GC% deviation from 50%, GC-clamp, homopolymer runs, ",
                "and hairpin/self-dimer risk. Lower score = better pair. ",
                "Click a row to highlight the primers on the barcode sequence below."
              )),
              DTOutput("primer_table"),
              uiOutput("primer_detail")
            )
          ),
          br(),
          fluidRow(
            column(12,
              h5("Barcode Sequence Viewer"),
              uiOutput("seq_viewer")
            )
          ),
          br(),
          fluidRow(
            column(12,
              h5(tags$span(class="glyphicon glyphicon-search"), " Primer-BLAST Specificity Results"),
              tags$p(tags$small(
                "After designing primers, select a pair and click ",
                tags$b("Check Selected Pair via Primer-BLAST"),
                " in the sidebar. Results are fetched live from NCBI and may take 30–90 seconds."
              )),
              tags$div(
                style = "border: 1px solid #dce0e6; border-radius:4px; padding:14px; background:#fafbfc; min-height:80px;",
                uiOutput("blast_result_panel")
              )
            )
          )
        )
      )
    )
  )
)

# ── Download button (only rendered after primers are run) ───────────────────
# (Defined inside server via uiOutput("download_primers_btn") above)
