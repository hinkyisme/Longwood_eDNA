# eDNA Pipeline — Two-Page Shiny App
# ui.R

library(shiny)
library(DECIPHER)
library(Biostrings)
library(DT)
library(tools)
library(shinythemes)

options(shiny.maxRequestSize = 500 * 1024^2)

ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),

  # ---------- Custom CSS ----------
  tags$head(
    tags$script(HTML("
      Shiny.addCustomMessageHandler('setActiveNav', function(page) {
        document.getElementById('go_builder').classList.toggle('active', page === 'builder');
        document.getElementById('go_identify').classList.toggle('active', page === 'identify');
      });
    ")),
    tags$style(HTML("
    /* ── Page shell ── */
    body { background: #f4f6f9; font-family: 'Georgia', serif; }

    /* ── Top nav bar ── */
    .edna-nav {
      display: flex;
      align-items: center;
      background: #1a2e44;
      padding: 0 28px;
      height: 58px;
      margin-bottom: 28px;
      box-shadow: 0 2px 8px rgba(0,0,0,.25);
    }
    .edna-nav .brand {
      color: #7ecfb3;
      font-size: 1.15rem;
      font-weight: bold;
      letter-spacing: .06em;
      margin-right: 36px;
      white-space: nowrap;
    }
    .edna-nav .nav-btn {
      background: transparent;
      border: none;
      color: #c8d8e8;
      font-size: .95rem;
      padding: 6px 18px;
      margin-right: 4px;
      border-radius: 4px;
      cursor: pointer;
      transition: background .15s, color .15s;
      letter-spacing: .03em;
    }
    .edna-nav .nav-btn:hover  { background: #243f5c; color: #fff; }
    .edna-nav .nav-btn.active { background: #7ecfb3; color: #1a2e44; font-weight: bold; }

    /* ── Section cards ── */
    .card {
      background: #fff;
      border-radius: 8px;
      padding: 20px 24px;
      margin-bottom: 20px;
      box-shadow: 0 1px 4px rgba(0,0,0,.08);
    }
    .card h4 {
      margin-top: 0;
      color: #1a2e44;
      border-bottom: 2px solid #7ecfb3;
      padding-bottom: 6px;
      font-size: 1rem;
      letter-spacing: .04em;
    }

    /* ── Run button ── */
    .run-btn {
      background: #1a2e44 !important;
      color: #7ecfb3 !important;
      border: none !important;
      font-size: 1rem !important;
      letter-spacing: .05em !important;
      border-radius: 6px !important;
      width: 100% !important;
      padding: 10px !important;
      margin-top: 6px !important;
      transition: background .2s !important;
    }
    .run-btn:hover { background: #243f5c !important; }

    /* ── Download button ── */
    .dl-btn { margin-top: 10px; width: 100%; }

    /* ── Status box ── */
    .status-box {
      background: #eaf7f2;
      border-left: 4px solid #7ecfb3;
      border-radius: 4px;
      padding: 10px 14px;
      font-size: .88rem;
      color: #1a2e44;
      margin-bottom: 12px;
    }

    hr { border-color: #dde3eb; margin: 14px 0; }
  "))),

  # ---------- Nav bar ----------
  div(class = "edna-nav",
    span(class = "brand", "🧬 eDNA Pipeline"),
    actionButton("go_builder",  "Barcode Builder",  class = "nav-btn active"),
    actionButton("go_identify", "Identify Taxa",    class = "nav-btn")
  ),

  # ---------- Page: Barcode Builder ----------
  uiOutput("page_builder"),

  # ---------- Page: Identify Taxa ----------
  uiOutput("page_identify")
)
