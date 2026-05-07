library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(dplyr)
library(readr)

# R/ directory is auto-sourced by Shiny:
#   read_multichannel.R  →  read_multichannel()
#   helpers.R            →  is_patient(), ntc_background()
#   results_tert.R       →  process_tert()

# ── Shared plot helper ────────────────────────────────────────────────────────
fa_lollipop <- function(r) {
  r <- r |>
    mutate(
      FA_plot = ifelse(is.na(FA), 0, FA),
      Sample  = reorder(Sample, FA_plot)
    )
  run_label <- paste(unique(r$`Run name`), collapse = ", ")
  pal <- c("Positive" = "#dc3545", "Negative" = "#198754", "Invalid" = "#adb5bd")

  ggplot(r, aes(x = FA_plot, y = Sample, colour = result)) +
    geom_segment(
      aes(x = 0, xend = FA_plot, yend = Sample),
      colour = "grey80", linewidth = 0.5
    ) +
    geom_point(aes(shape = valid), size = 8) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 4),
      labels = c("TRUE" = "Valid", "FALSE" = "Invalid"),
      name   = NULL
    ) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               colour = "grey40", linewidth = 0.6) +
    scale_colour_manual(values = pal, name = "Result") +
    scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x       = "Fractional Abundance (FA %)",
      y       = NULL,
      title   = run_label,
      caption = "Dashed line = FA 0.5% mutation threshold.  \u00d7 = invalid well."
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- page_fluid(
  theme = bs_theme(bootswatch = "flatly") |>
    bs_add_rules("
      .app-header { background-color: #2c3e50; color: #fff;
                    padding: 14px 24px 12px; margin-bottom: 0; }
      .app-header h4 { margin: 0; font-weight: 600; letter-spacing: 0.03em; }
      .nav-tabs { padding: 0 20px; background: #fff;
                  border-bottom: 2px solid #dee2e6; }
      .nav-tabs .nav-link { color: #555; font-weight: 500; }
      .nav-tabs .nav-link.active {
        color: #2c3e50 !important; font-weight: 700;
        border-bottom: 3px solid #2c3e50 !important;
        border-top: none; border-left: none; border-right: none;
        background: transparent;
      }
    "),

  # Title bar
  div(class = "app-header", tags$h4("SIGRID \u00b7 Uromonitor")),

  # Assay tabs
  navset_tab(
    id = "assay_tab",

    # ── Guide tab ─────────────────────────────────────────────────────────────
    nav_panel(
      "\u2139\ufe0f  Guide",
      div(
        style = "max-width:780px; margin: 36px auto; padding: 0 16px;",

        h3("Welcome to SIGRID \u00b7 Uromonitor"),
        p("This app processes digital PCR results from the ",
          tags$strong("Uromonitor\u00ae"), " kit and returns per-sample mutation calls,
          fractional abundance values, and QC flags for the SIGRID study."),

        tags$hr(),

        h5("\U0001f4c2 Required input file"),
        p("Each assay module expects a single ",
          tags$strong("channel_based_multichannel CSV"), " exported from the
          QuantStudio\u2122 Absolute Q\u2122 software after a completed run."),
        p("The file should have the following columns:"),
        tags$table(
          class = "table table-sm table-bordered",
          style = "font-size:0.9em;",
          tags$thead(tags$tr(
            tags$th("Column"), tags$th("Example value"), tags$th("Notes")
          )),
          tags$tbody(
            tags$tr(tags$td("Run"),     tags$td("TERT_146_plate3_reruns_260414"), tags$td("Plate/run identifier")),
            tags$tr(tags$td("Sample"),  tags$td("UM_2025_22_2"),                  tags$td("Patient ID or NTC / PC")),
            tags$tr(tags$td("Well"),    tags$td("A3"),                            tags$td("")),
            tags$tr(tags$td("Channels"),tags$td("FAM  /  VIC  /  FAM VIC"),      tags$td("One row per channel per sample")),
            tags$tr(tags$td("Count"),   tags$td("137"),                           tags$td("Positive droplet count for that channel")),
            tags$tr(tags$td("Total"),   tags$td("20383"),                         tags$td("Total accepted droplets")  )
          )
        ),

        tags$hr(),

        h5("\U0001f9ea Validity criteria (QuantStudio Absolute Q\u2122)"),
        tags$ul(
          tags$li(tags$strong("Valid well:"), " Total events \u2265 20,000 AND IC positive events \u2265 10"),
          tags$li(tags$strong("Mutation positive:"), " NTC-subtracted FAM positives \u2265 5 AND FA \u2265 0.5%"),
          tags$li(tags$strong("FA ="), " FAM\u2082 / (FAM\u2082 + IC\u2082) \u00d7 100  (raw positives, not NTC-subtracted)")
        ),

        tags$hr(),

        h5("\U0001f4cb How to use"),
        tags$ol(
          tags$li("Click the ", tags$strong("TERT-146"), " tab (or FGFR3 when available)."),
          tags$li("Click ", tags$strong("Browse"), " in the sidebar and select your ", tags$em("channel_based_multichannel.csv"), " file."),
          tags$li("Results appear immediately across three sub-tabs:
            Summary boxes, Results table, FA Plot, and Controls."),
          tags$li("Use ", tags$strong("Download Results"), " /", tags$strong(" Download Controls"),
                  " to export the processed data as CSV files.")
        ),

        tags$hr(),

        h5("\u26a0\ufe0f Notes"),
        tags$ul(
          tags$li("Controls: NTC requires FAM = 0; PC requires valid IC (VIC \u2265 10)."),
          tags$li("Samples flagged as ", tags$strong("Invalid"), " have insufficient droplets or IC — they should be repeated per kit guidelines."),
          tags$li("Results follow the manufacturer\u2019s algorithm exactly (Uromonitor\u00ae Instructions v3.5, section 9.3.4).")
        ),

        tags$hr(),
        p(class = "text-muted", style = "font-size:0.85em;",
          "Uromonitor\u00ae is a CE-IVD product by Infogene Lda., Coimbra, Portugal. ",
          "This app is for research use within the SIGRID study only.")
      )
    ),

    # ── TERT-146 tab ──────────────────────────────────────────────────────────
    nav_panel(
      "TERT-146",
      layout_sidebar(
        sidebar = sidebar(
          width = 270,
          fileInput(
            "tert_file", "Upload multichannel CSV",
            accept = ".csv", buttonLabel = "Browse",
            placeholder = "No file selected"
          ),
          hr(),
          uiOutput("tert_dl_ui")
        ),
        uiOutput("tert_summary_row"),
        br(),
        navset_card_tab(
          nav_panel("Results",  DTOutput("tert_tbl_results")),
          nav_panel("FA Plot",  plotOutput("tert_plot_fa", height = "520px")),
          nav_panel("Controls", DTOutput("tert_tbl_controls"))
        )
      )
    ),

    # ── FGFR3 248/249 tab (placeholder) ───────────────────────────────────────
    nav_panel(
      "FGFR3 248/249",
      layout_sidebar(
        sidebar = sidebar(
          width = 270,
          fileInput(
            "fgfr3_file", "Upload multichannel CSV",
            accept = ".csv", buttonLabel = "Browse",
            placeholder = "No file selected"
          )
        ),
        card(
          card_body(
            class = "text-center py-5",
            h4("FGFR3 248/249 module — coming soon"),
            p("This analysis module is under development.", class = "text-muted")
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── TERT reactives ──────────────────────────────────────────────────────────
  tert_processed <- reactive({
    req(input$tert_file)
    process_tert(input$tert_file$datapath)
  })

  tert_results  <- reactive(tert_processed()$results)
  tert_controls <- reactive(tert_processed()$controls)

  # Summary boxes
  output$tert_summary_row <- renderUI({
    req(tert_results())
    r <- tert_results()
    layout_columns(
      value_box("Samples",  nrow(r),                        theme = "secondary"),
      value_box("Positive", sum(r$result == "Positive"),    theme = "danger"),
      value_box("Negative", sum(r$result == "Negative"),    theme = "success"),
      value_box("Invalid",  sum(r$result == "Invalid"),     theme = "warning"),
      col_widths = c(3, 3, 3, 3)
    )
  })

  # Results table
  output$tert_tbl_results <- renderDT({
    req(tert_results())
    tert_results() |>
      mutate(FA = round(FA, 3)) |>
      datatable(rownames = FALSE, selection = "none",
                options = list(pageLength = 25, dom = "ftp")) |>
      formatStyle("result",
        backgroundColor = styleEqual(
          c("Positive", "Negative", "Invalid"),
          c("#f8d7da",  "#d1e7dd",  "#fff3cd")
        ))
  })

  # FA plot
  output$tert_plot_fa <- renderPlot({
    req(tert_results())
    fa_lollipop(tert_results())
  })

  # Controls table
  output$tert_tbl_controls <- renderDT({
    req(tert_controls())
    tert_controls() |>
      datatable(rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE)) |>
      formatStyle("control_pass",
        backgroundColor = styleEqual(c(TRUE, FALSE), c("#d1e7dd", "#f8d7da")))
  })

  # Download buttons
  output$tert_dl_ui <- renderUI({
    req(tert_results())
    tagList(
      downloadButton("tert_dl_results",  "Download Results",  class = "btn-sm w-100 mb-2"),
      downloadButton("tert_dl_controls", "Download Controls", class = "btn-sm w-100")
    )
  })

  make_filename <- function(df, suffix) {
    run <- gsub("[^A-Za-z0-9_]", "_", unique(df$`Run name`)[1])
    paste0(run, suffix, ".csv")
  }

  output$tert_dl_results <- downloadHandler(
    filename = function() make_filename(tert_results(), "_results"),
    content  = function(f) write_csv(tert_results(), f)
  )

  output$tert_dl_controls <- downloadHandler(
    filename = function() make_filename(tert_controls(), "_controls"),
    content  = function(f) write_csv(tert_controls(), f)
  )
}

shinyApp(ui, server)
