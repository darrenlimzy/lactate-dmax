# Lactate step test: LT1 (Log-log) and LT2 (Log-Poly-ModDmax) web app.
# The calculation lives in analysis.R (copied from the canonical script).
# Runs in the browser through Shinylive: entered data never leave the
# viewer's computer.

library(shiny)
library(bslib)
library(jsonlite)
source("analysis.R")

example_data <- paste(
  "step\tlength_min\tpower_w\tlactate\thr\tpartial",
  "0\t0\t0\t1.1\t72\tN",
  "1\t4\t100\t1.0\t112\tN",
  "2\t4\t125\t0.9\t119\tN",
  "3\t4\t150\t1.0\t127\tN",
  "4\t4\t175\t1.2\t136\tN",
  "5\t4\t200\t1.5\t145\tN",
  "6\t4\t225\t2.1\t155\tN",
  "7\t4\t250\t3.0\t165\tN",
  "8\t4\t275\t4.6\t175\tN",
  "9\t4\t300\t7.2\t185\tN",
  "10\t4\t\t9.1\t190\tY",
  sep = "\n"
)

# Reads pasted text (tabs from Excel, or commas) into the df_full shape
parse_steps <- function(txt) {
  txt <- trimws(txt)
  if (!nzchar(txt)) stop("Paste the stage data first.")
  sep <- if (grepl("\t", txt)) "\t" else ","
  raw <- read.table(text = txt, sep = sep, header = TRUE,
                    stringsAsFactors = FALSE, strip.white = TRUE,
                    na.strings = c("", "NA"), blank.lines.skip = TRUE)
  names(raw) <- tolower(trimws(names(raw)))
  need <- c("step", "length_min", "power_w", "lactate", "hr", "partial")
  missing <- setdiff(need, names(raw))
  if (length(missing))
    stop("Missing column(s): ", paste(missing, collapse = ", "),
         ". Use the header row from the template.")
  partial <- toupper(substr(trimws(as.character(raw$partial)), 1, 1)) %in% c("Y", "T")
  if (sum(partial) > 1) stop("Only one row can be marked partial = Y.")
  out <- tibble::tibble(
    step       = as.integer(raw$step),
    length     = as.integer(raw$length_min),
    intensity  = suppressWarnings(as.numeric(raw$power_w)),
    lactate    = as.numeric(raw$lactate),
    heart_rate = as.integer(round(as.numeric(raw$hr))),
    is_partial = partial
  )
  if (anyNA(out$lactate) || anyNA(out$heart_rate))
    stop("Every row needs a lactate and a heart rate value.")
  if (anyNA(out$intensity[!out$is_partial]))
    stop("Every row except the partial stage needs a power value.")
  if (!any(out$step == 0))
    stop("Add the resting baseline as step 0 with power 0.")
  if (sum(out$step > 0 & !out$is_partial) < 4)
    stop("At least 4 completed stages are needed for a 3rd-degree polynomial fit.")
  out
}

# SpexSG brand colours (SpexSG PPT template): red E12E3F, orange F79120
spex_red    <- "#E12E3F"
spex_orange <- "#F79120"

spex_theme <- bs_theme(
  version = 5,
  primary = spex_red,
  secondary = "#595959",
  success = "#2E7D32",
  "link-color" = spex_red,
  base_font = font_collection("Avenir Next", "Avenir", "Nunito Sans", "Helvetica Neue", "Arial", "sans-serif"),
  heading_font = font_collection("Avenir Next", "Avenir", "Nunito Sans", "Helvetica Neue", "Arial", "sans-serif"),
  "headings-font-weight" = 600
)

spex_css <- HTML(sprintf("
  .bslib-page-title.navbar, .navbar.bslib-page-title { background: #fff !important;
    border-bottom: 0; box-shadow: none; padding: 10px 16px 0; }
  .bslib-page-title .container-fluid, .bslib-page-title > div { display: block; width: 100%%; }
  .brand-header { width: 100%%; }
  .brand-bar { display: flex; align-items: center; gap: 18px; flex-wrap: wrap; width: 100%%; }
  .brand-bar img { height: 56px; width: auto; }
  .brand-bar .brand-divider { width: 1px; height: 40px; background: #D9D9D9; }
  .brand-bar .brand-title { color: #222; font-size: 1.25rem; font-weight: 600; line-height: 1.2; }
  .brand-bar .brand-sub { color: #595959; font-size: .85rem; font-weight: 400; }
  .brand-stripe { height: 6px; width: calc(100%% + 32px); margin: 10px -16px 0;
    background: linear-gradient(90deg, %1$s 0%%, %1$s 55%%, %2$s 100%%); }
  .nav-tabs .nav-link.active, .nav-tabs .nav-link.active:hover { color: %1$s !important;
    border-top: 3px solid %1$s; font-weight: 600; }
  .nav-tabs .nav-link { color: #595959; }
  .card-header { background: #F2F2F2; border-bottom: 1px solid #D9D9D9; font-weight: 600; }
  .bslib-sidebar-layout > .sidebar { background: #FAFAFA; border-right: 3px solid %1$s; }
  .sidebar h6 { color: %1$s; font-weight: 700; margin-top: .5rem; }
  .btn-default:not(.btn-primary) { border-color: %1$s; color: %1$s; }
  .btn-default:not(.btn-primary):hover { background: %1$s; color: #fff; }
  #run { color: #fff; font-weight: 600; letter-spacing: .02em; }
  table.table th { color: #fff; background: %1$s; }
  @media (max-width: 576px) { .brand-bar img { height: 32px; } .brand-bar .brand-divider { display: none; } }
", spex_red, spex_orange))

ui <- page_sidebar(
  title = div(
    class = "brand-header",
    div(
      class = "brand-bar",
      img(src = "spexsg.png", alt = "SpexSG"),
      img(src = "hpsi.png", alt = "HPSI"),
      div(class = "brand-divider"),
      div(
        div(class = "brand-title", "Lactate Step Test: LT1 and LT2"),
        div(class = "brand-sub", "Log-log LT1 and Log-Poly-ModDmax LT2 | HPSI Physiology")
      )
    ),
    div(class = "brand-stripe")
  ),
  window_title = "Lactate Step Test | HPSI Physiology",
  theme = spex_theme,
  tags$head(tags$style(spex_css)),
  sidebar = sidebar(
    width = 340,
    textInput("athlete", "Athlete name", "Example athlete"),
    numericInput("hr_max", "Max heart rate (bpm)", 195, min = 100, max = 230),
    h6("Final partial stage (for Wmax)"),
    numericInput("last_complete", "Last completed stage power (W)", 300),
    numericInput("next_stage", "Partial stage power (W)", 325),
    numericInput("stage_sec", "Full stage duration (s)", 240),
    numericInput("held_sec", "Time held in partial stage (s)", 90, min = 0),
    checkboxInput("include_partial", "Include the partial stage in the curve fit", TRUE),
    helpText("Untick if the steep top-end point makes the polynomial wiggle.",
             "The Wmax line still shows either way."),
    actionButton("run", "Calculate", class = "btn-primary")
  ),
  navset_card_tab(
    id = "results_tabs",
    nav_panel(
      "1. Data",
      p("Paste the stage data below, straight from Excel (tab-separated) or as CSV.",
        "Keep the header row. Step 0 is the resting baseline (power 0).",
        "Mark the final partial stage with partial = Y and leave its power blank:",
        "it is set to Wmax automatically."),
      fluidRow(
        column(6, fileInput("csv", NULL, accept = c(".csv", ".txt", ".tsv"),
                            buttonLabel = "Load CSV...", placeholder = "or paste below")),
        column(6, downloadButton("template", "Download CSV template", class = "btn-sm mt-1"))
      ),
      textAreaInput("steps", NULL, example_data, width = "100%", height = "330px",
                    resize = "vertical"),
      uiOutput("parse_msg")
    ),
    nav_panel(
      "2. Results",
      uiOutput("error_msg"),
      layout_columns(
        col_widths = c(5, 7),
        card(card_header("Thresholds"), tableOutput("lt_table"),
             uiOutput("wmax_line")),
        card(card_header("Training zones"), tableOutput("zone_table"))
      ),
      card(full_screen = TRUE, card_header("Lactate curve"),
           plotOutput("curve", height = "480px")),
      card(
        card_header("Downloads"),
        div(
          class = "d-flex flex-wrap gap-2",
          downloadButton("dl_curve", "Curve (PNG)"),
          downloadButton("dl_tables", "Tables (PNG)"),
          downloadButton("dl_csv", "Results (CSV)"),
          downloadButton("dl_json", "Results (JSON)")
        )
      )
    ),
    nav_panel(
      "About the method",
      markdown("
**LT1** is the Log-log breakpoint (`loglog_restrainer = 0.75`).
**LT2** is the Log-Poly-ModDmax point: the modified Dmax on a 3rd-degree
polynomial fit, starting from the point before the first rise found by the
log-log method. Both are computed with the R package
[lactater](https://cran.r-project.org/package=lactater) (`sport = \"cycling\"`,
baseline excluded from the fit).

**Wmax** = last completed power + (partial stage power - last completed power)
x time held / stage duration.

**Zones** (Burnley & Jones 2007 intensity domains):

| Zone | Power range |
|---|---|
| Z1 Light | 50% to <70% LT1 |
| Z2 Moderate | 70% LT1 to LT1 |
| Z3 Heavy | LT1 to 95% LT2 |
| Z4 LT2 | 95% to 105% LT2 |
| Z5 Severe | >105% LT2 |

Zone lactate comes from the fitted curve. Zone heart rate is linearly
interpolated between the measured stage values.

**Privacy.** The app runs entirely in your browser (R compiled to WebAssembly).
The data you enter are not sent to any server.
")
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$csv, {
    txt <- paste(readLines(input$csv$datapath, warn = FALSE), collapse = "\n")
    updateTextAreaInput(session, "steps", value = txt)
  })

  output$template <- downloadHandler(
    filename = "lactate_step_test_template.csv",
    content = function(file) writeLines(gsub("\t", ",", example_data), file)
  )

  wmax <- reactive({
    req(input$last_complete, input$next_stage, input$stage_sec)
    compute_wmax(input$last_complete, input$next_stage,
                 input$stage_sec, input$held_sec %||% 0)
  })

  parsed <- reactive(tryCatch(parse_steps(input$steps), error = function(e) e))

  output$parse_msg <- renderUI({
    d <- parsed()
    if (inherits(d, "error"))
      div(class = "alert alert-warning", conditionMessage(d))
    else
      div(class = "text-success small",
          sum(d$step > 0 & !d$is_partial), " completed stages",
          if (any(d$is_partial)) " + 1 partial stage" else "",
          " read. Press Calculate.")
  })

  result <- eventReactive(input$run, {
    d <- parsed()
    if (inherits(d, "error")) return(d)
    W <- wmax()
    d$intensity[d$is_partial] <- W
    out <- tryCatch(
      suppressMessages(analyse_step_test(d, input$athlete, input$hr_max, W,
                                         input$include_partial)),
      error = function(e) simpleError(paste("The analysis failed:", conditionMessage(e))))
    if (inherits(out, "error")) return(out)
    out$Wmax <- W
    out$data <- d
    out$meta <- list(athlete = input$athlete, hr_max = input$hr_max,
                     include_partial = input$include_partial)
    out
  })

  observeEvent(input$run, nav_select("results_tabs", "2. Results"))

  ok <- reactive({
    r <- result()
    req(!inherits(r, "error"))
    r
  })

  output$error_msg <- renderUI({
    r <- result()
    if (inherits(r, "error")) div(class = "alert alert-danger", conditionMessage(r))
  })

  output$lt_table   <- renderTable(
    dplyr::mutate(ok()$threshold_table, `Heart Rate (BPM)` = as.integer(`Heart Rate (BPM)`)),
    digits = 1, align = "c")
  output$zone_table <- renderTable(ok()$zones, align = "c")
  output$wmax_line  <- renderUI(p(strong("Wmax: "), round(ok()$Wmax, 1), " W"))
  output$curve      <- renderPlot(ok()$curve_plot, res = 96)

  stem <- function() gsub("[^A-Za-z0-9]+", "_", trimws(ok()$meta$athlete))

  output$dl_curve <- downloadHandler(
    filename = function() paste0(stem(), "_lactate_curve.png"),
    content = function(file)
      ggsave(file, ok()$curve_plot, width = 10, height = 6, dpi = 200, bg = "white")
  )
  output$dl_tables <- downloadHandler(
    filename = function() paste0(stem(), "_lactate_tables.png"),
    content = function(file)
      ggsave(file, ok()$tables_grob, width = 7, height = 5, dpi = 200, bg = "white")
  )
  output$dl_csv <- downloadHandler(
    filename = function() paste0(stem(), "_lactate_results.csv"),
    content = function(file) {
      r <- ok()
      lt <- r$threshold_table
      lines <- c(
        paste0("Athlete,", r$meta$athlete),
        paste0("HRmax (bpm),", r$meta$hr_max),
        paste0("Wmax (W),", round(r$Wmax, 1)),
        paste0("Partial stage in fit,", r$meta$include_partial),
        "",
        "Threshold,Power (W),Lactate (mmol/L),HR (bpm)",
        paste(lt$LT, lt$`Intensity (W)`, lt$`Lactate (mmol/L)`, lt$`Heart Rate (BPM)`, sep = ","),
        "",
        "Zone,Power (W),Lactate (mmol/L),HR (bpm),%HRmax",
        apply(r$zones, 1, function(z) paste0('"', z, '"', collapse = ",")),
        "",
        "Data entered",
        "step,length_min,power_w,lactate,hr,partial",
        with(r$data, paste(step, length, round(intensity, 1), lactate, heart_rate,
                           ifelse(is_partial, "Y", "N"), sep = ","))
      )
      writeLines(lines, file)
    }
  )
  output$dl_json <- downloadHandler(
    filename = function() paste0(stem(), "_lactate_results.json"),
    content = function(file) {
      r <- ok()
      writeLines(toJSON(results_list(r, r$meta$athlete, r$meta$hr_max, r$Wmax,
                                     r$meta$include_partial),
                        auto_unbox = TRUE, pretty = TRUE, digits = NA), file)
    }
  )
}

shinyApp(ui, server)
