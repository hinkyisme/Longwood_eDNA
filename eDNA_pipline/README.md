# eDNA Barcode Builder & Identifier — Shiny App

A complete pipeline for building and applying molecular barcodes for non-model organisms from mitochondrial or genomic data.

---

## Required R Packages

```r
install.packages(c("shiny", "DT", "tools", "shinycssloaders", "shinythemes", "ape"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DECIPHER", "Biostrings"))
```

---

## How to Run

Place both `ui.R` and `server.R` in the same folder and run:

```r
shiny::runApp("path/to/shiny_app/")
```

Or open either file in RStudio and click **Run App**.

---

## Pipeline Overview

### Input Mode 1 — GenBank Flat File (.gb)
Upload a GenBank flat file (single or multi-record). The app automatically:
- Parses all records
- Finds features annotated as `16S rRNA`, `16S ribosomal RNA`, or similar
- Extracts the exact nucleotide coordinates
- Handles `complement()` strand features
- Uses the `/organism=` qualifier as the sequence name

### Input Mode 2 — FASTA / FASTQ
Upload a FASTA/FASTQ file of full mitogenomes or pre-assembled sequences. The app:
- Applies the forward and reverse primers (with mismatch tolerance) to trim out the barcode region
- Filters by length (150–3000 bp)

### Steps (both modes)
| Step | Description |
|------|-------------|
| 1 | Load & extract reference barcodes |
| 2 | Train IDTAXA classifier (DECIPHER) |
| 3 | Load & classify test sequences |
| 4 | Multiple sequence alignment (AlignSeqs) + QC |

---

## Tabs

| Tab | Contents |
|-----|----------|
| **Barcode Extraction** | Count, length stats, first 10 names; download extracted FASTA |
| **Taxonomy Results** | DT table with confidence-colored cells; download CSV & alignment |
| **Alignment QC** | Column occupancy plot, summary statistics |
| **Log** | Timestamped step-by-step execution log |

---

## Default Primers (16S — Mytilus / mussel)
- Forward: `GGGTCACCAACTCCGCTAAC`
- Reverse:  `GGGCAGCTAAGGCTGGAAAA`

These can be replaced with any locus primers (COI, 12S, cytB, etc.).

---

## Notes
- Max upload size: 500 MB
- Memory: `gc()` is called after training; lower K-mer (e.g. 6–7) reduces RAM usage for large reference sets
- The alignment step can be slow for >200 sequences; consider pre-filtering references
