# Lactate step test: LT1 (Log-log) and LT2 (Log-Poly-ModDmax) web app.
# The calculation lives in analysis.R (copied from the canonical script).
# Runs in the browser through Shinylive: entered data never leave the
# viewer's computer.

library(shiny)
library(bslib)
library(jsonlite)
library(rhandsontable)
source("analysis.R")

# Grid columns as shown to the user, and the names used in the code
grid_cols <- c("Step", "Stage length (s)", "Power (W)", "Lactate (mmol/L)",
               "HR (bpm)", "Partial *")
grid_vars <- c("step", "length_s", "power_w", "lactate", "hr", "partial")

example_grid <- data.frame(
  step     = 0:10,
  length_s = c(0, rep(240, 9), 90),
  power_w  = c(0, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325),
  lactate  = c(1.1, 1.0, 0.9, 1.0, 1.2, 1.5, 2.1, 3.0, 4.6, 7.2, 9.1),
  hr       = c(72, 112, 119, 127, 136, 145, 155, 165, 175, 185, 190),
  partial  = c(rep(FALSE, 10), TRUE)
)

blank_grid <- data.frame(
  step = c(0, 1:8), length_s = c(0, rep(NA, 8)), power_w = c(0, rep(NA, 8)),
  lactate = NA_real_, hr = NA_real_, partial = FALSE
)

# Reads a CSV (template layout, or tab-separated from Excel) into the grid
read_steps_file <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  sep <- if (grepl("\t", txt)) "\t" else ","
  raw <- read.table(text = txt, sep = sep, header = TRUE,
                    stringsAsFactors = FALSE, strip.white = TRUE,
                    na.strings = c("", "NA"), blank.lines.skip = TRUE)
  names(raw) <- tolower(trimws(names(raw)))
  if (!"length_s" %in% names(raw) && "length_min" %in% names(raw))
    raw$length_s <- as.numeric(raw$length_min) * 60
  missing <- setdiff(grid_vars, names(raw))
  if (length(missing))
    stop("Missing column(s): ", paste(missing, collapse = ", "),
         ". Use the header row from the CSV template.")
  data.frame(
    step     = as.numeric(raw$step),
    length_s = as.numeric(raw$length_s),
    power_w  = suppressWarnings(as.numeric(raw$power_w)),
    lactate  = as.numeric(raw$lactate),
    hr       = as.numeric(raw$hr),
    partial  = toupper(substr(trimws(as.character(raw$partial)), 1, 1)) %in% c("Y", "T", "1")
  )
}

# Drops the empty spare rows the grid keeps at the bottom
filled_rows <- function(g) {
  g[!(is.na(g$step) & is.na(g$length_s) & is.na(g$power_w) &
        is.na(g$lactate) & is.na(g$hr)), , drop = FALSE]
}

# Checks the grid and turns it into the df_full shape. The partial row
# holds the power the unfinished stage was set at and, as its stage
# length, the time held; Wmax comes from those two numbers.
steps_from_grid <- function(g) {
  g <- filled_rows(g)
  if (!nrow(g)) stop("Enter the stage data first.")
  partial <- as_tick(g$partial)
  if (sum(partial) > 1) stop("Tick Partial for one row only (the final, unfinished stage).")
  out <- tibble::tibble(
    step       = as.integer(g$step),
    length     = as.numeric(g$length_s) / 60,   # minutes, as in the canonical script
    intensity  = as.numeric(g$power_w),
    lactate    = as.numeric(g$lactate),
    heart_rate = as.integer(round(as.numeric(g$hr))),
    is_partial = partial
  )
  if (anyNA(out$step)) stop("Every row needs a step number.")
  if (anyDuplicated(out$step)) stop("Each step number can only appear once.")
  if (anyNA(out$intensity)) stop("Every row needs a power value (0 for the baseline).")
  if (anyNA(out$lactate) || anyNA(out$heart_rate))
    stop("Every row needs a lactate and a heart rate value.")
  if (!any(out$step == 0))
    stop("Add the resting baseline as step 0 with power 0.")
  out <- out[order(out$step), ]
  done <- out[out$step > 0 & !out$is_partial, ]
  if (nrow(done) < 4)
    stop("At least 4 completed stages are needed for a 3rd-degree polynomial fit.")

  last_complete <- done[nrow(done), ]
  if (any(out$is_partial)) {
    p <- out[out$is_partial, ]
    if (p$step != max(out$step)) stop("The partial stage must be the last step.")
    if (is.na(p$length) || p$length <= 0)
      stop("Enter the time held in the partial stage (its Stage length, in seconds).")
    if (is.na(last_complete$length) || last_complete$length <= 0)
      stop("Enter the Stage length of the completed stages (e.g. 240 s).")
    if (p$intensity <= last_complete$intensity)
      stop("The partial stage's power (", p$intensity, " W) must be higher than the last ",
           "completed stage (", last_complete$intensity, " W): enter the power it was set at.")
    if (p$length >= last_complete$length)
      stop("The partial stage's time held must be shorter than a full stage (",
           round(last_complete$length * 60), " s).")
    held_s <- p$length * 60
    full_s <- last_complete$length * 60
    W <- compute_wmax(last_complete$intensity, p$intensity, full_s, held_s)
    how <- sprintf("%g W + (%g - %g W) x %g / %g s", last_complete$intensity,
                   p$intensity, last_complete$intensity, held_s, full_s)
  } else {
    W <- last_complete$intensity
    how <- "no partial stage: last completed stage"
  }
  attr(out, "Wmax") <- W
  attr(out, "Wmax_how") <- how
  out
}

# Adds a blank stage after the last completed one (before the partial row,
# if there is one) with the next step number and the same stage length
add_stage <- function(g) {
  g <- filled_rows(g)
  if (!nrow(g)) return(blank_grid)
  part <- which(as_tick(g$partial))
  done <- if (length(part)) g[-part, , drop = FALSE] else g
  last <- done[nrow(done), ]
  new <- data.frame(
    step = if (nrow(done)) max(done$step, na.rm = TRUE) + 1 else 0,
    length_s = if (nrow(done) && !is.na(last$step) && last$step > 0) last$length_s else NA,
    power_w = NA_real_, lactate = NA_real_, hr = NA_real_, partial = FALSE
  )
  if (length(part)) {
    p <- g[part, , drop = FALSE]
    p$step <- new$step + 1
    rbind(done, new, p)
  } else rbind(g, new)
}

# Partial comes back from the grid as TRUE/FALSE, 0/1 or text; make it logical
as_tick <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "y", "yes")

render_grid <- function(df) {
  df <- data.frame(
    step = as.numeric(df$step), length_s = as.numeric(df$length_s),
    power_w = as.numeric(df$power_w), lactate = as.numeric(df$lactate),
    hr = as.numeric(df$hr), partial = as_tick(df$partial)
  )
  names(df) <- grid_cols
  rhandsontable(df, rowHeaders = NULL, stretchH = "all",
                height = 30 + 23 * (nrow(df) + 1),
                minSpareRows = 1, useTypes = TRUE) |>
    hot_col("Step", format = "0") |>
    hot_col("Stage length (s)", format = "0") |>
    hot_col("Power (W)", format = "0") |>
    hot_col("Lactate (mmol/L)", format = "0.0[0]") |>
    hot_col("HR (bpm)", format = "0") |>
    hot_col("Partial *", type = "checkbox", halign = "htCenter") |>
    hot_cols(colWidths = c(60, 120, 100, 120, 90, 80), halign = "htCenter") |>
    hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
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
  .data-card { max-width: 820px; border-top: 3px solid %1$s; }
  .handsontable th { background: #F2F2F2; color: #222; font-weight: 600; }
  .handsontable td { font-variant-numeric: tabular-nums; }
  .btn-add { background: %1$s !important; color: #fff !important; border-color: %1$s !important; }
  .wmax-box { margin: .25rem 0 .75rem; padding: .6rem .8rem; background: #fff;
    border: 1px solid #D9D9D9; border-left: 3px solid %1$s; border-radius: 4px; }
  .wmax-label { font-size: .8rem; color: #595959; font-weight: 600; text-transform: uppercase; }
  .wmax-value { font-size: 1.5rem; font-weight: 700; color: %1$s; line-height: 1.2; }
  .wmax-how { font-size: .78rem; color: #595959; }
  .partial-note { margin-top: .75rem; padding: .6rem .8rem; font-size: .85rem; color: #333;
    background: #FFF5F6; border-left: 3px solid %1$s; border-radius: 4px; }
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
    uiOutput("wmax_box"),
    checkboxInput("include_partial", "Include the partial stage in the curve fit", TRUE),
    helpText("Untick if the steep top-end point makes the polynomial wiggle.",
             "The Wmax line still shows either way."),
    actionButton("run", "Calculate", class = "btn-primary")
  ),
  navset_card_tab(
    id = "results_tabs",
    nav_panel(
      "1. Data",
      card(
        class = "data-card",
        card_header(
          div(
            class = "d-flex justify-content-between align-items-center flex-wrap gap-2 w-100",
            span("Stage data"),
            div(
              class = "d-flex flex-wrap gap-2",
              actionButton("load_example", "Example", class = "btn-sm"),
              actionButton("clear_grid", "Clear", class = "btn-sm"),
              downloadButton("template", "CSV template", class = "btn-sm")
            )
          )
        ),
        card_body(
          fillable = FALSE,
          p(class = "small text-muted mb-2",
            "One row per stage. Type into the grid, or paste straight from Excel: click the",
            "top-left cell, then Ctrl+V (Cmd+V on a Mac). Use Add stage for another row;",
            "right-click a row to insert or delete it. Step 0 is the resting baseline",
            "(length 0, power 0)."),
          rHandsontableOutput("grid"),
          actionButton("add_stage_below", "Add stage", icon = icon("plus"),
                       class = "btn-sm btn-add mt-2"),
          div(
            class = "partial-note",
            strong("* Partial: "),
            "tick the last stage if the athlete started it but could not finish it.",
            "For that row, enter the power the stage was set at and, as its Stage length,",
            "how long the athlete held it (in seconds). The app then works out Wmax:",
            "last completed power + (partial stage power - last completed power) x time held",
            "/ full stage length, and plots the partial point at Wmax. Tick one row only;",
            "leave it unticked if the test ended on a completed stage."
          ),
          div(
            class = "mt-3",
            fileInput("csv", NULL, accept = c(".csv", ".txt", ".tsv"), width = "100%",
                      buttonLabel = "Load CSV...",
                      placeholder = "Optional: a CSV saved from the template")
          ),
          uiOutput("parse_msg")
        )
      )
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
      card(full_screen = TRUE, card_header("Log-log visual inspection (LT1)"),
           plotOutput("loglog", height = "480px"),
           p(class = "small text-muted mb-0",
             "Use this to check that LT1 sits at a sensible break in the data. Grey line:",
             "the stage data as lactater interpolates them. Black lines: the two-segment fit.",
             "Shaded: the top 25% of the data, which the fit leaves out.")),
      card(
        card_header("Downloads"),
        div(
          class = "d-flex flex-wrap gap-2",
          downloadButton("dl_curve", "Curve (PNG)"),
          downloadButton("dl_loglog", "Log-log (PNG)"),
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
x time held / full stage length. The partial stage's power and time held come
from the row ticked Partial; the full stage length is that of the last completed
stage. The partial point is plotted (and fitted) at Wmax.

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

  # The grid is redrawn only when a new table is loaded; edits made in the
  # grid come back through input$grid. `current` holds whichever is newest.
  grid_seed <- reactiveVal(list(df = example_grid, n = 0))
  current   <- reactiveVal(example_grid)
  load_grid <- function(df) {
    grid_seed(list(df = df, n = grid_seed()$n + 1))
    current(df)
  }

  output$grid <- renderRHandsontable(render_grid(grid_seed()$df))

  observeEvent(input$grid, {
    g <- hot_to_r(input$grid)
    names(g) <- grid_vars
    current(g)
    # Typing or pasting past the last row adds rows; redraw so the grid
    # grows to fit them instead of scrolling inside its box
    if (nrow(filled_rows(g)) > nrow(filled_rows(grid_seed()$df))) load_grid(filled_rows(g))
  })

  observeEvent(input$load_example, load_grid(example_grid))
  observeEvent(input$add_stage_below, load_grid(add_stage(current())))
  observeEvent(input$clear_grid, load_grid(blank_grid))
  observeEvent(input$csv, {
    df <- tryCatch(read_steps_file(input$csv$datapath), error = function(e) e)
    if (inherits(df, "error"))
      showNotification(conditionMessage(df), type = "error", duration = 8)
    else load_grid(df)
  })

  output$template <- downloadHandler(
    filename = "lactate_step_test_template.csv",
    content = function(file) {
      t <- example_grid
      t$power_w <- ifelse(is.na(t$power_w), "", t$power_w)
      t$partial <- ifelse(t$partial, "Y", "N")
      write.csv(t, file, row.names = FALSE, quote = FALSE)
    }
  )

  parsed <- reactive(tryCatch(steps_from_grid(current()), error = function(e) e))

  output$wmax_box <- renderUI({
    d <- parsed()
    if (inherits(d, "error"))
      div(class = "wmax-box", div(class = "wmax-label", "Wmax"),
          div(class = "wmax-how", "Complete the stage data to see Wmax."))
    else
      div(class = "wmax-box",
          div(class = "wmax-label", "Wmax"),
          div(class = "wmax-value", round(attr(d, "Wmax"), 1), " W"),
          div(class = "wmax-how", attr(d, "Wmax_how")))
  })

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
    W <- attr(d, "Wmax")
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
  output$loglog     <- renderPlot(ok()$loglog_plot, res = 96)

  stem <- function() gsub("[^A-Za-z0-9]+", "_", trimws(ok()$meta$athlete))

  output$dl_curve <- downloadHandler(
    filename = function() paste0(stem(), "_lactate_curve.png"),
    content = function(file)
      ggsave(file, ok()$curve_plot, width = 10, height = 6, dpi = 200, bg = "white")
  )
  output$dl_loglog <- downloadHandler(
    filename = function() paste0(stem(), "_loglog.png"),
    content = function(file)
      ggsave(file, ok()$loglog_plot, width = 10, height = 6, dpi = 200, bg = "white")
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
        "step,length_s,power_w,lactate,hr,partial",
        with(r$data, paste(step, round(length * 60), round(intensity, 1), lactate, heart_rate,
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
