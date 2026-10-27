# =============================================================
# Lactate Step Test — LT1, LT2 and 5-Zone Derivation
# =============================================================
# The calculation below is copied from the canonical Lactate
# Analyzer V2 script (cycling-lab-test skill,
# assets/lactate_analyzer_v2.R) and wrapped in a function so the
# web app can call it. Any change to the method must be made in
# the canonical script first and copied here.
#
# Only difference: the plot's left x-limit is min(95, first stage
# power - 10) instead of a fixed 95 W, so tests that start below
# 105 W still show every point. Numbers are not affected.
# =============================================================

library(lactater)
library(dplyr)
library(ggplot2)
library(tibble)
library(gridExtra)
library(grid)

# Wmax from the final partial stage
compute_wmax <- function(last_complete_power, next_stage_power,
                         stage_duration_sec, time_completed_sec) {
  slope <- (next_stage_power - last_complete_power) / stage_duration_sec
  last_complete_power + slope * time_completed_sec
}

# df_full: step, length, intensity, lactate, heart_rate, is_partial
# (the partial row's intensity must already be set to Wmax)
analyse_step_test <- function(df_full, athlete_name, hr_max, Wmax,
                              include_partial = TRUE) {

  # df used for curve fitting and thresholds (respects the toggle)
  df <- if (include_partial) df_full else dplyr::filter(df_full, !is_partial)

  # ----------------------------
  # Fit the lactate curve
  # ----------------------------
  curve_data <- lactate_curve(
    .data = df,
    intensity_column = "intensity",
    lactate_column = "lactate",
    heart_rate_column = "heart_rate",
    fit = "3rd degree polynomial",
    include_baseline = FALSE,
    sport = "cycling"
  )

  # ----------------------------
  # Compute LT1 and LT2
  # ----------------------------
  thresholds <- lactate_threshold(
    .data = df,
    intensity_column = "intensity",
    lactate_column = "lactate",
    heart_rate_column = "heart_rate",
    method = c("Log-log", "Dmax"),
    fit = "3rd degree polynomial",
    include_baseline = FALSE,
    sport = "cycling",
    loglog_restrainer = 0.75,
    plot = FALSE
  ) %>%
    filter(method %in% c("Log-log", "Log-Poly-ModDmax")) %>%
    arrange(intensity) %>%
    mutate(LT = c("LT1", "LT2"))

  # ----------------------------
  # Interpolation helpers
  # ----------------------------
  lactate_at <- function(power) {
    approx(
      curve_data$lactate_curve$intensity,
      curve_data$lactate_curve$lactate,
      xout = power,
      rule = 2
    )$y
  }

  hr_at <- function(power) {
    approx(df$intensity, df$heart_rate, xout = power, rule = 2)$y
  }

  fmt_range <- function(lo, hi, dp = 0) {
    lo_r <- round(lo, dp)
    hi_r <- round(hi, dp)
    ifelse(
      lo_r == hi_r,
      format(lo_r, nsmall = dp),
      paste0(format(lo_r, nsmall = dp), " - ", format(hi_r, nsmall = dp))
    )
  }

  # ----------------------------
  # 5-Zone Training Zone Boundaries
  # ----------------------------
  # Burnley & Jones (2007) intensity-domain mapping:
  # T1 - Light:    50%  to <70% LT1   | Active recovery
  # T2 - Moderate: 70%  LT1 to LT1    | Moderate domain (endurance / Z2)
  # T3 - Heavy:    LT1  to 95% LT2    | Heavy domain (tempo / sweet-spot)
  # T4 - LT2:      95%  to 105% LT2   | Threshold band around MLSS
  # T5 - Severe:   > 105% LT2         | Severe domain (VO2max / anaerobic)

  LT1_power <- thresholds$intensity[thresholds$LT == "LT1"]
  LT2_power <- thresholds$intensity[thresholds$LT == "LT2"]

  zone_bounds <- tibble::tribble(
    ~Zone,           ~p_lo,                       ~p_hi,
    "Z1 - Light",    round(0.50 * LT1_power),     round(0.70 * LT1_power) - 1,
    "Z2 - Moderate", round(0.70 * LT1_power),     round(LT1_power),
    "Z3 - Heavy",    round(LT1_power) + 1,        round(0.95 * LT2_power),
    "Z4 - LT2",      round(0.95 * LT2_power) + 1, round(1.05 * LT2_power),
    "Z5 - Severe",   round(1.05 * LT2_power) + 1, NA_real_
  )

  # ----------------------------
  # Threshold summary table
  # ----------------------------
  threshold_table <- thresholds %>%
    transmute(
      LT,
      `Intensity (W)`    = round(intensity, 1),
      `Lactate (mmol/L)` = round(lactate, 1),
      `Heart Rate (BPM)` = round(heart_rate, 0)
    )

  # ----------------------------
  # Zones display table
  # ----------------------------
  zones <- zone_bounds %>%
    mutate(
      `Power (W)` = ifelse(
        is.na(p_hi),
        paste0(">", p_lo - 1),
        paste0(p_lo, " - ", p_hi)
      ),
      lac_lo = lactate_at(p_lo),
      lac_hi = ifelse(is.na(p_hi), NA_real_, lactate_at(p_hi)),
      `Lactate (mmol/L)` = ifelse(
        is.na(lac_hi),
        paste0(">", round(lac_lo, 1)),
        fmt_range(lac_lo, lac_hi, dp = 1)
      ),
      hr_lo = hr_at(p_lo),
      hr_hi = ifelse(is.na(p_hi), NA_real_, hr_at(p_hi)),
      `HR (bpm)` = ifelse(
        is.na(hr_hi),
        paste0(">", round(hr_lo)),
        fmt_range(hr_lo, hr_hi, dp = 0)
      ),
      `%HRmax` = ifelse(
        is.na(hr_hi),
        paste0(">", round(hr_lo / hr_max * 100)),
        fmt_range(hr_lo / hr_max * 100, hr_hi / hr_max * 100, dp = 0)
      )
    ) %>%
    select(Zone, `Power (W)`, `Lactate (mmol/L)`, `HR (bpm)`, `%HRmax`)

  # ----------------------------
  # Zone rectangles for the plot
  # ----------------------------
  x_right <- max(df_full$intensity)
  x_left  <- min(95, min(df_full$intensity[df_full$step > 0]) - 10)

  zone_rects <- zone_bounds %>%
    mutate(
      label     = sub(" .*", "", Zone),
      rect_xmin = p_lo,
      rect_xmax = ifelse(is.na(p_hi), Inf, p_hi),
      label_x   = ifelse(is.na(p_hi),
                         (p_lo + x_right) / 2,
                         (p_lo + p_hi) / 2),
      fill      = c("#9E9E9E", "#42A5F5", "#66BB6A", "#FFCA28", "#EF5350")
    )

  # ----------------------------
  # Lactate-power plot with full-height zone bands
  # ----------------------------
  y_axis_max   <- max(df_full$lactate) + 0.8
  band_top     <- max(df_full$lactate) + 0.5
  zone_label_y <- band_top * 0.93

  p <- ggplot() +
    geom_rect(
      data = zone_rects,
      aes(xmin = rect_xmin, xmax = rect_xmax,
          ymin = -Inf, ymax = band_top,
          fill = label),
      alpha = 0.35
    ) +
    geom_segment(
      data = thresholds,
      aes(x = intensity, xend = intensity, y = 0, yend = lactate),
      linetype = "dashed",
      linewidth = 0.8,
      color = "grey15"
    ) +
    geom_vline(
      xintercept = Wmax,
      linetype = "dotted",
      linewidth = 0.7,
      color = "grey25"
    ) +
    geom_line(
      data = curve_data$lactate_curve,
      aes(x = intensity, y = lactate),
      linewidth = 1
    ) +
    geom_point(
      data = df_full %>% filter(step > 0),
      aes(x = intensity, y = lactate, shape = is_partial),
      size = 2
    ) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1)) +
    geom_point(
      data = thresholds,
      aes(x = intensity, y = lactate),
      size = 3
    ) +
    geom_text(
      data = thresholds,
      aes(x = intensity, y = lactate, label = LT),
      size = 4.5,
      fontface = "bold",
      color = "grey10",
      vjust = -1
    ) +
    geom_text(
      data = zone_rects,
      aes(x = label_x, y = zone_label_y, label = label),
      size = 5,
      fontface = "bold",
      color = "grey15"
    ) +
    scale_fill_manual(values = setNames(zone_rects$fill, zone_rects$label)) +
    scale_y_continuous(limits = c(0, y_axis_max), expand = c(0, 0)) +
    coord_cartesian(xlim = c(x_left, x_right + 5)) +
    labs(
      x = "Power (W)",
      y = "Blood lactate (mmol/L)",
      title = paste0("Lactate curve with LT1, LT2 and training zones — ", athlete_name),
      subtitle = paste0("Wmax = ", round(Wmax, 1), " W (dotted line)")
    ) +
    theme_light() +
    theme(legend.position = "none")

  # ----------------------------
  # Table grobs
  # ----------------------------
  tbl_theme <- ttheme_minimal(
    base_size = 10,
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5),
      bg_params = list(fill = "white", col = "grey70")
    ),
    colhead = list(
      fg_params = list(fontface = "bold"),
      bg_params = list(fill = "grey95", col = "grey70")
    )
  )
  tables_combined <- arrangeGrob(
    tableGrob(threshold_table, rows = NULL, theme = tbl_theme),
    tableGrob(zones, rows = NULL, theme = tbl_theme),
    ncol = 1,
    heights = c(1, 1.4)
  )

  list(
    thresholds      = thresholds,
    threshold_table = threshold_table,
    zones           = zones,
    LT1_power       = LT1_power,
    LT2_power       = LT2_power,
    curve_plot      = p,
    tables_grob     = tables_combined
  )
}

# Same structure as the canonical script's lactate_results.json
results_list <- function(res, athlete_name, hr_max, Wmax, include_partial) {
  th <- res$thresholds
  list(
    athlete_name    = athlete_name,
    hr_max          = hr_max,
    Wmax            = round(Wmax, 1),
    include_partial = include_partial,
    thresholds = list(
      LT1 = list(
        power_w = round(res$LT1_power, 1),
        lactate = round(th$lactate[th$LT == "LT1"], 1),
        hr_bpm  = round(th$heart_rate[th$LT == "LT1"], 0)
      ),
      LT2 = list(
        power_w = round(res$LT2_power, 1),
        lactate = round(th$lactate[th$LT == "LT2"], 1),
        hr_bpm  = round(th$heart_rate[th$LT == "LT2"], 0)
      )
    ),
    zones = lapply(seq_len(nrow(res$zones)), function(i) {
      list(
        zone       = res$zones$Zone[i],
        power_w    = res$zones$`Power (W)`[i],
        lactate_mm = res$zones$`Lactate (mmol/L)`[i],
        hr_bpm     = res$zones$`HR (bpm)`[i],
        pct_hrmax  = res$zones$`%HRmax`[i]
      )
    })
  )
}
