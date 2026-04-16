############################################################
## 02_cross_sectional_analysis.R
##
## Purpose:
## Cross-sectional analyses for the ASAP manuscript.
##
## Expected input objects from 01_prepare_analysis_dataset.R:
## - analysis_data
## - cross_data_cog_master
##
## Main outputs:
## - table1
## - table2
## - results_cross
## - results_med_export
## - results_export_clean
##
## Notes:
## This script assumes that analysis-ready datasets have already
## been created in memory by running 01_prepare_analysis_dataset.R.
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260304)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(broom)
  library(purrr)
  library(stringr)
  library(scales)
  library(grid)
  library(ragg)
  library(writexl)
  library(openxlsx)
  library(mediation)
})

############################################################
## Input checks
############################################################

if (!exists("analysis_data")) {
  stop("Object `analysis_data` not found. Please run 01_prepare_analysis_dataset.R first.")
}

if (!exists("cross_data_cog_master")) {
  stop("Object `cross_data_cog_master` not found. Please run 01_prepare_analysis_dataset.R first.")
}

############################################################
## Helper functions
############################################################

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

fmt_mean_sd <- function(x) {
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_median_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.2f (%.2f, %.2f)", q[2], q[1], q[3])
}

fmt_n_pct <- function(x, val) {
  denom <- sum(!is.na(x))
  n <- sum(x == val, na.rm = TRUE)
  sprintf("%d/%d (%.1f%%)", n, denom, n / denom * 100)
}

clean_missing <- function(x) {
  x[x %in% c(-9, -8, -7, 88, 99, 888, 999, 998, 9998, 9999)] <- NA
  x
}

fmt_ci <- function(est, low, high) {
  sprintf("%.2f (%.2f, %.2f)", est, low, high)
}

force_utf8 <- function(x) {
  if (is.character(x)) {
    x <- iconv(x, from = "", to = "UTF-8", sub = "")
  }
  x
}

############################################################
## Labels
############################################################

get_exposure_label <- function(x) {
  dplyr::case_when(
    x == "SpO2"   ~ "SpO2",
    x == "lnULK1" ~ "ln(ULK1)",
    TRUE ~ x
  )
}

get_outcome_label <- function(x) {
  dplyr::case_when(
    x == "RAVLT_tot"   ~ "RAVLT total recall",
    x == "RAVLT_learn" ~ "RAVLT learning",
    x == "Stroop_t"    ~ "Stroop interference time",
    x == "Stroop_r"    ~ "Stroop interference ratio",
    TRUE ~ x
  )
}

############################################################
## ASAP plotting system
############################################################

ASAP_colors <- list(
  SpO2_fill = "#B7E1CD",
  SpO2_line = "#2E8B57",
  ULK1_fill = "#F2B3AA",
  ULK1_line = "#C54A3F",
  points = "#000000",
  ci_alpha = 0.25
)

ASAP_linewidth <- 1
ASAP_hist_linewidth <- 0.3
ASAP_tag_size <- 10
ASAP_pointsize <- 1.4

theme_ASAP <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "#333333"),
    axis.line = element_line(color = "#333333"),
    legend.position = "none"
  )

theme_ASAP_diag <- theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 8, color = "#333333"),
    axis.line = element_line(color = "#333333"),
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4)
  )

############################################################
## 1. Table 1 and Table 2
############################################################

analysis_data_tbl <- analysis_data %>%
  mutate(
    sex = clean_missing(sex),
    education = clean_missing(education),
    smoke_exposure = clean_missing(smoke_exposure),
    hypertension = clean_missing(hypertension),
    diabetes = clean_missing(diabetes),
    CVD = clean_missing(CVD),
    cpap_use = clean_missing(cpap_use),
    obesity = clean_missing(obesity),
    high_LDL = clean_missing(high_LDL),
    AHI_bin5 = clean_missing(AHI_bin5)
  )

table1 <- bind_rows(
  data.frame(
    Characteristic = "Age, years",
    Missing = sum(is.na(analysis_data_tbl$age)),
    ASAP1 = fmt_mean_sd(analysis_data_tbl$age)
  ),
  data.frame(
    Characteristic = "Male sex",
    Missing = sum(is.na(analysis_data_tbl$sex)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$sex, 1)
  ),
  data.frame(
    Characteristic = "Higher education, yes",
    Missing = sum(is.na(analysis_data_tbl$education)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$education, 1)
  ),
  data.frame(
    Characteristic = "BMI, kg/m2",
    Missing = sum(is.na(analysis_data_tbl$BMI)),
    ASAP1 = fmt_mean_sd(analysis_data_tbl$BMI)
  ),
  data.frame(
    Characteristic = "Obesity (BMI ≥30), yes",
    Missing = sum(is.na(analysis_data_tbl$obesity)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$obesity, 1)
  ),
  data.frame(
    Characteristic = "High LDL (≥3 mmol/L), yes",
    Missing = sum(is.na(analysis_data_tbl$high_LDL)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$high_LDL, 1)
  ),
  data.frame(
    Characteristic = "Smoking exposure (ever), yes",
    Missing = sum(is.na(analysis_data_tbl$smoke_exposure)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$smoke_exposure, 1)
  ),
  data.frame(
    Characteristic = "Hypertension, yes",
    Missing = sum(is.na(analysis_data_tbl$hypertension)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$hypertension, 1)
  ),
  data.frame(
    Characteristic = "Type 2 diabetes mellitus, yes",
    Missing = sum(is.na(analysis_data_tbl$diabetes)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$diabetes, 1)
  ),
  data.frame(
    Characteristic = "Cardiovascular disease, yes",
    Missing = sum(is.na(analysis_data_tbl$CVD)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$CVD, 1)
  ),
  data.frame(
    Characteristic = "CPAP use, yes",
    Missing = sum(is.na(analysis_data_tbl$cpap_use)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$cpap_use, 1)
  ),
  data.frame(
    Characteristic = "AHI ≥5, yes",
    Missing = sum(is.na(analysis_data_tbl$AHI_bin5)),
    ASAP1 = fmt_n_pct(analysis_data_tbl$AHI_bin5, 1)
  ),
  data.frame(
    Characteristic = "Total sleep time, hours",
    Missing = sum(is.na(analysis_data_tbl$TST)),
    ASAP1 = fmt_mean_sd(analysis_data_tbl$TST)
  ),
  data.frame(
    Characteristic = "Sleep efficiency, %",
    Missing = sum(is.na(analysis_data_tbl$sleep_efficiency)),
    ASAP1 = fmt_mean_sd(analysis_data_tbl$sleep_efficiency)
  )
)

table2 <- bind_rows(
  data.frame(
    Characteristic = "AHI, mean (SD)",
    Missing = sum(is.na(analysis_data_tbl$AHI_raw)),
    ASAP1 = fmt_mean_sd(analysis_data_tbl$AHI_raw)
  ),
  data.frame(
    Characteristic = "AHI, median (IQR)",
    Missing = sum(is.na(analysis_data_tbl$AHI_raw)),
    ASAP1 = fmt_median_iqr(analysis_data_tbl$AHI_raw)
  ),
  data.frame(
    Characteristic = "SpO2, mean (SD)",
    Missing = sum(is.na(analysis_data_tbl$SpO2)),
    ASAP1 = fmt_mean_sd(analysis_data_tbl$SpO2)
  ),
  data.frame(
    Characteristic = "SpO2, median (IQR)",
    Missing = sum(is.na(analysis_data_tbl$SpO2)),
    ASAP1 = fmt_median_iqr(analysis_data_tbl$SpO2)
  ),
  data.frame(
    Characteristic = "Serum ULK1, mean (SD)",
    Missing = sum(is.na(analysis_data_tbl$ULK1)),
    ASAP1 = fmt_mean_sd(analysis_data_tbl$ULK1)
  ),
  data.frame(
    Characteristic = "Serum ULK1, median (IQR)",
    Missing = sum(is.na(analysis_data_tbl$ULK1)),
    ASAP1 = fmt_median_iqr(analysis_data_tbl$ULK1)
  ),
  data.frame(
    Characteristic = "Ln ULK1, median (IQR)",
    Missing = sum(is.na(analysis_data_tbl$lnULK1)),
    ASAP1 = fmt_median_iqr(analysis_data_tbl$lnULK1)
  ),
  data.frame(
    Characteristic = "RAVLT total recall",
    Missing = sum(is.na(analysis_data$RAVLT_tot)),
    ASAP1 = fmt_mean_sd(analysis_data$RAVLT_tot)
  ),
  data.frame(
    Characteristic = "RAVLT learning",
    Missing = sum(is.na(analysis_data$RAVLT_learn)),
    ASAP1 = fmt_mean_sd(analysis_data$RAVLT_learn)
  ),
  data.frame(
    Characteristic = "Stroop interference time",
    Missing = sum(is.na(analysis_data$Stroop_t)),
    ASAP1 = fmt_mean_sd(analysis_data$Stroop_t)
  ),
  data.frame(
    Characteristic = "Stroop interference ratio",
    Missing = sum(is.na(analysis_data$Stroop_r)),
    ASAP1 = fmt_mean_sd(analysis_data$Stroop_r)
  )
)

############################################################
## 2. Supplementary Fig. S1
############################################################

plot_data_s1 <- analysis_data %>%
  dplyr::select(ULK1, lnULK1)

p1 <- ggplot(plot_data_s1, aes(x = ULK1)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = ASAP_colors$ULK1_fill,
    color = "grey40",
    linewidth = ASAP_hist_linewidth,
    alpha = 0.9
  ) +
  geom_density(
    color = ASAP_colors$ULK1_line,
    linewidth = ASAP_linewidth
  ) +
  labs(
    title = "ULK1 (raw)",
    x = "Serum ULK1",
    y = "Density"
  ) +
  theme_ASAP

p2 <- ggplot(plot_data_s1, aes(x = lnULK1)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = ASAP_colors$ULK1_fill,
    color = "grey40",
    linewidth = ASAP_hist_linewidth,
    alpha = 0.9
  ) +
  geom_density(
    color = ASAP_colors$ULK1_line,
    linewidth = ASAP_linewidth
  ) +
  labs(
    title = "ln(ULK1)",
    x = "ln(Serum ULK1)",
    y = "Density"
  ) +
  theme_ASAP

Figure_S1 <- (p1 | p2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = ASAP_tag_size, face = "bold"))

shapiro_ulk1 <- shapiro.test(na.omit(analysis_data$ULK1))
shapiro_lnulk1 <- shapiro.test(na.omit(analysis_data$lnULK1))

############################################################
## 3. Cross-sectional main analyses
############################################################

cross_data <- cross_data_cog_master %>%
  filter(!is.na(SpO2) | !is.na(lnULK1))

required_vars <- c(
  "SpO2", "lnULK1",
  "RAVLT_tot", "RAVLT_learn", "Stroop_t", "Stroop_r",
  "age", "sex", "education", "obesity", "smoke_exposure",
  "diabetes", "CVD", "hypertension",
  "TST", "sleep_efficiency", "cpap_use", "AHI_bin5"
)

missing_vars <- setdiff(required_vars, names(cross_data))
if (length(missing_vars) > 0) {
  stop(
    paste0(
      "The following required variables are missing in `cross_data`: ",
      paste(missing_vars, collapse = ", ")
    )
  )
}

cross_data <- cross_data %>%
  mutate(
    age = to_num(age),
    sex = to_num(sex),
    education = to_num(education),
    obesity = to_num(obesity),
    smoke_exposure = to_num(smoke_exposure),
    diabetes = to_num(diabetes),
    CVD = to_num(CVD),
    hypertension = to_num(hypertension),
    TST = to_num(TST),
    sleep_efficiency = to_num(sleep_efficiency),
    cpap_use = to_num(cpap_use),
    AHI_bin5 = to_num(AHI_bin5),
    SpO2 = to_num(SpO2),
    lnULK1 = to_num(lnULK1),
    RAVLT_tot = to_num(RAVLT_tot),
    RAVLT_learn = to_num(RAVLT_learn),
    Stroop_t = to_num(Stroop_t),
    Stroop_r = to_num(Stroop_r)
  )

exposures <- c("SpO2", "lnULK1")
outcomes <- c("RAVLT_tot", "RAVLT_learn", "Stroop_t", "Stroop_r")

model_formulas <- list(
  Model1 = "~ EXPOSURE",
  Model2 = "~ EXPOSURE + age + sex",
  Model3 = "~ EXPOSURE + age + sex + education + obesity + smoke_exposure",
  Model4 = "~ EXPOSURE + age + sex + education + obesity + smoke_exposure + diabetes + CVD + hypertension",
  Model5 = "~ EXPOSURE + age + sex + education + obesity + smoke_exposure + diabetes + CVD + hypertension + TST + sleep_efficiency + cpap_use + AHI_bin5"
)

run_cross_model <- function(exposure, outcome, formula_string, data) {
  f <- as.formula(
    paste0(outcome, " ", gsub("EXPOSURE", exposure, formula_string))
  )
  
  model <- lm(f, data = data)
  model_n <- nobs(model)
  
  broom::tidy(model, conf.int = TRUE) %>%
    filter(term == exposure) %>%
    mutate(
      exposure = exposure,
      outcome = outcome,
      N = model_n
    )
}

results_cross_list <- list()

for (outcome in outcomes) {
  for (exposure in exposures) {
    for (model_name in names(model_formulas)) {
      res <- try(
        run_cross_model(
          exposure = exposure,
          outcome = outcome,
          formula_string = model_formulas[[model_name]],
          data = cross_data
        ),
        silent = TRUE
      )
      
      if (!inherits(res, "try-error") && !is.null(res) && nrow(res) > 0) {
        res$model <- model_name
        results_cross_list[[length(results_cross_list) + 1]] <- res
      }
    }
  }
}

results_cross <- bind_rows(results_cross_list)

if (nrow(results_cross) == 0) {
  stop("No cross-sectional results were generated.")
}

results_cross <- results_cross %>%
  transmute(
    exposure = exposure,
    outcome = outcome,
    model = model,
    N = N,
    Beta = estimate,
    SE = std.error,
    CI_low = conf.low,
    CI_high = conf.high,
    t = statistic,
    p = p.value
  ) %>%
  arrange(exposure, outcome, model) %>%
  mutate(
    exposure_label = get_exposure_label(exposure),
    outcome_label = get_outcome_label(outcome),
    AHI_adjustment = case_when(
      model %in% c("Model1", "Model2", "Model3", "Model4") ~ "No AHI adjustment",
      model == "Model5" ~ "Model 5 + AHI binary (<5 vs >=5)",
      TRUE ~ NA_character_
    ),
    Beta_CI = sprintf("%.2f (%.2f, %.2f)", Beta, CI_low, CI_high),
    p_text = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

results_cross_main <- results_cross %>%
  filter(model %in% c("Model1", "Model2", "Model3", "Model4"))

results_cross_model5 <- results_cross %>%
  filter(model == "Model5")

############################################################
## 4. Fig. 1 forest plot
############################################################

forest_data <- results_cross %>%
  mutate(
    model = factor(model, levels = c("Model1", "Model2", "Model3", "Model4", "Model5")),
    outcome_label = factor(
      outcome_label,
      levels = c(
        "RAVLT total recall",
        "RAVLT learning",
        "Stroop interference time",
        "Stroop interference ratio"
      )
    ),
    exposure = factor(exposure, levels = c("SpO2", "lnULK1"))
  ) %>%
  arrange(exposure, outcome_label, model) %>%
  mutate(
    beta_ci = sprintf("%.2f (%.2f, %.2f)", Beta, CI_low, CI_high),
    p_text  = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
    sig = p < 0.05
  )

make_y_positions_block <- function(df) {
  model_slots <- c(17, 13, 9, 5, 1)
  block_height <- 24
  out_levels <- levels(df$outcome_label)
  y_vals <- numeric(nrow(df))
  
  for (i in seq_along(out_levels)) {
    idx <- which(df$outcome_label == out_levels[i])
    block_base <- (length(out_levels) - i) * block_height
    y_vals[idx] <- block_base + model_slots
  }
  
  y_vals
}

plot_forest_row_centered <- function(data, exposure_name, color, tag_letter) {
  df <- data %>%
    filter(exposure == exposure_name) %>%
    arrange(outcome_label, model)
  
  df$y <- make_y_positions_block(df)
  
  outcome_df <- df %>%
    group_by(outcome_label) %>%
    summarise(
      y_center = mean(y),
      outcome_display = str_wrap(first(as.character(outcome_label)), width = 16),
      .groups = "drop"
    )
  
  y_top <- max(df$y) + 6
  
  max_abs <- max(abs(c(df$CI_low, df$CI_high)), na.rm = TRUE)
  x_lim <- ceiling(max_abs * 10) / 10 + 0.2
  x_lim <- max(x_lim, 1.0)
  
  tag_panel <- wrap_elements(
    full = textGrob(
      tag_letter,
      x = 0.5, y = 0.98,
      gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "Arial")
    )
  )
  
  outcome_panel <- ggplot(outcome_df, aes(y = y_center)) +
    geom_text(
      aes(x = 0, label = outcome_display),
      hjust = 0.5,
      family = "Arial",
      size = 3.1
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Outcome",
      hjust = 0.5,
      family = "Arial",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(0, y_top + 1), breaks = NULL) +
    coord_cartesian(xlim = c(-0.9, 0.9), clip = "off") +
    theme_void()
  
  model_panel <- ggplot(df, aes(y = y)) +
    geom_text(
      aes(x = 0, label = model),
      hjust = 0.5,
      family = "Arial",
      size = 3
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Model",
      hjust = 0.5,
      family = "Arial",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(0, y_top + 1), breaks = NULL) +
    coord_cartesian(xlim = c(-0.7, 0.7), clip = "off") +
    theme_void()
  
  forest_panel <- ggplot(
    df,
    aes(y = y, x = Beta, xmin = CI_low, xmax = CI_high)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.4
    ) +
    geom_point(
      data = df %>% filter(sig),
      aes(x = Beta, y = y),
      inherit.aes = FALSE,
      size = 4.2,
      color = scales::alpha(color, 0.18)
    ) +
    geom_errorbarh(
      height = 0.12,
      linewidth = 0.5,
      color = color
    ) +
    geom_point(
      size = 1.6,
      color = color
    ) +
    scale_x_continuous(
      limits = c(-x_lim, x_lim),
      breaks = pretty(c(-x_lim, x_lim), n = 5)
    ) +
    scale_y_continuous(
      limits = c(0, y_top + 1),
      breaks = NULL
    ) +
    labs(
      title = unique(df$exposure_label),
      x = "Beta (95% CI)",
      y = NULL
    ) +
    theme_ASAP +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  beta_panel <- ggplot(df, aes(y = y)) +
    geom_text(
      aes(x = 0, label = beta_ci, fontface = ifelse(sig, "bold", "plain")),
      hjust = 0.5,
      family = "Courier",
      size = 2.8
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Beta (95% CI)",
      hjust = 0.5,
      family = "Arial",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(0, y_top + 1), breaks = NULL) +
    coord_cartesian(xlim = c(-1.1, 1.1), clip = "off") +
    theme_void()
  
  p_panel <- ggplot(df, aes(y = y)) +
    geom_text(
      aes(x = 0, label = p_text, fontface = ifelse(sig, "bold", "plain")),
      hjust = 0.5,
      family = "Courier",
      size = 2.8
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "p value",
      hjust = 0.5,
      family = "Arial",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(0, y_top + 1), breaks = NULL) +
    coord_cartesian(xlim = c(-0.45, 0.45), clip = "off") +
    theme_void()
  
  tag_panel +
    outcome_panel +
    model_panel +
    forest_panel +
    beta_panel +
    p_panel +
    plot_layout(widths = c(0.35, 1.8, 1.2, 4.8, 2.1, 0.9))
}

row_A <- plot_forest_row_centered(
  forest_data, "SpO2", ASAP_colors$SpO2_line, "A"
)

row_B <- plot_forest_row_centered(
  forest_data, "lnULK1", ASAP_colors$ULK1_line, "B"
)

Figure_cross_forest <- row_A / row_B

############################################################
## 5. Supplementary Fig. S3 diagnostics
############################################################

model5_formula <- "~ EXPOSURE + age + sex +
                    education + obesity + smoke_exposure +
                    diabetes + CVD + hypertension +
                    TST + sleep_efficiency + cpap_use +
                    AHI_bin5"

get_exposure_colors <- function(exposure) {
  if (exposure == "SpO2") {
    return(list(fill = ASAP_colors$SpO2_fill, line = ASAP_colors$SpO2_line))
  }
  if (exposure == "lnULK1") {
    return(list(fill = ASAP_colors$ULK1_fill, line = ASAP_colors$ULK1_line))
  }
  list(fill = "grey80", line = "grey40")
}

build_lm_diag_subplots <- function(data, exposure, outcome) {
  f <- as.formula(
    paste0(outcome, " ", gsub("EXPOSURE", exposure, model5_formula))
  )
  
  model <- lm(f, data = data)
  df_diag <- broom::augment(model)
  cols <- get_exposure_colors(exposure)
  
  p1 <- ggplot(df_diag, aes(x = .fitted, y = .resid)) +
    geom_point(size = ASAP_pointsize, alpha = 0.7, color = cols$fill) +
    geom_smooth(method = "loess", se = FALSE, color = cols$line, linewidth = 0.8) +
    labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
    theme_ASAP_diag
  
  p2 <- ggplot(df_diag, aes(sample = .std.resid)) +
    stat_qq(size = ASAP_pointsize, alpha = 0.7, color = cols$fill) +
    stat_qq_line(color = cols$line, linewidth = 0.8) +
    labs(title = "Normal Q-Q", x = "Theoretical quantiles", y = "Standardized residuals") +
    theme_ASAP_diag
  
  p3 <- ggplot(df_diag, aes(x = .fitted, y = sqrt(abs(.std.resid)))) +
    geom_point(size = ASAP_pointsize, alpha = 0.7, color = cols$fill) +
    geom_smooth(method = "loess", se = FALSE, color = cols$line, linewidth = 0.8) +
    labs(title = "Scale-Location", x = "Fitted values", y = "Sqrt(|Standardized residuals|)") +
    theme_ASAP_diag
  
  p4 <- ggplot(df_diag, aes(x = .hat, y = .std.resid)) +
    geom_point(size = ASAP_pointsize, alpha = 0.7, color = cols$fill) +
    geom_smooth(method = "loess", se = FALSE, color = cols$line, linewidth = 0.8) +
    labs(title = "Residuals vs Leverage", x = "Leverage", y = "Standardized residuals") +
    theme_ASAP_diag
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

build_refined_model_panel <- function(data, exposure, outcome, panel_id) {
  subs <- build_lm_diag_subplots(data, exposure, outcome)
  
  panel_title <- paste0(
    get_exposure_label(exposure),
    " x ",
    get_outcome_label(outcome)
  )
  
  tag_grob <- wrap_elements(
    full = textGrob(
      label = panel_id,
      x = 0, y = 0.5,
      just = c("left", "center"),
      gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "sans")
    )
  )
  
  title_grob <- wrap_elements(
    full = textGrob(
      label = panel_title,
      x = 0.5, y = 0.5,
      gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "sans")
    )
  )
  
  top_strip <- tag_grob + title_grob +
    plot_layout(widths = c(0.08, 0.92))
  
  inner_patch <- (subs$p1 | subs$p2) / (subs$p3 | subs$p4)
  inner_grob  <- patchwork::patchworkGrob(inner_patch)
  
  border_grob <- rectGrob(
    x = 0.5, y = 0.5,
    width = 0.985, height = 0.985,
    just = "center",
    gp = gpar(fill = NA, col = "grey45", lwd = 1.2, lty = "dashed")
  )
  
  boxed_diag <- wrap_elements(
    full = grobTree(inner_grob, border_grob)
  )
  
  top_strip / boxed_diag +
    plot_layout(heights = c(0.09, 0.91)) &
    theme(plot.margin = margin(8, 8, 8, 8))
}

panel_map <- tibble::tribble(
  ~panel_id, ~exposure, ~outcome,      ~caption_text,
  "A",       "SpO2",    "RAVLT_tot",   "SpO2 x RAVLT total recall",
  "B",       "SpO2",    "RAVLT_learn", "SpO2 x RAVLT learning",
  "C",       "SpO2",    "Stroop_t",    "SpO2 x Stroop interference time",
  "D",       "SpO2",    "Stroop_r",    "SpO2 x Stroop interference ratio",
  "E",       "lnULK1",  "RAVLT_tot",   "ln(ULK1) x RAVLT total recall",
  "F",       "lnULK1",  "RAVLT_learn", "ln(ULK1) x RAVLT learning",
  "G",       "lnULK1",  "Stroop_t",    "ln(ULK1) x Stroop interference time",
  "H",       "lnULK1",  "Stroop_r",    "ln(ULK1) x Stroop interference ratio"
)

diag_panels <- vector("list", length = nrow(panel_map))

for (i in seq_len(nrow(panel_map))) {
  diag_panels[[i]] <- build_refined_model_panel(
    data = cross_data,
    exposure = panel_map$exposure[i],
    outcome = panel_map$outcome[i],
    panel_id = panel_map$panel_id[i]
  )
}

plot_list <- list(
  diag_panels[[1]], plot_spacer(), diag_panels[[2]],
  diag_panels[[3]], plot_spacer(), diag_panels[[4]],
  diag_panels[[5]], plot_spacer(), diag_panels[[6]],
  diag_panels[[7]], plot_spacer(), diag_panels[[8]]
)

Figure_S3_diag <- wrap_plots(
  plotlist = plot_list,
  ncol = 3,
  byrow = TRUE,
  widths = c(1, 0.01, 1),
  heights = c(1, 1, 1, 1)
)

############################################################
## 6. Mediation analyses
############################################################

exposures_main <- c("SpO2")
mediator_main  <- "lnULK1"

exposure_reverse  <- "lnULK1"
mediators_reverse <- c("SpO2")

covariates_main <- c(
  "age", "sex", "education",
  "obesity", "smoke_exposure",
  "diabetes", "CVD", "hypertension",
  "TST", "sleep_efficiency", "cpap_use",
  "AHI_bin5"
)

required_vars_med <- unique(c(
  exposures_main,
  mediator_main,
  exposure_reverse,
  mediators_reverse,
  outcomes,
  covariates_main
))

missing_vars_med <- setdiff(required_vars_med, names(cross_data))

if (length(missing_vars_med) > 0) {
  stop(
    paste0(
      "Missing required mediation variables in `cross_data`: ",
      paste(missing_vars_med, collapse = ", ")
    )
  )
}

run_mediation <- function(data, exposure, mediator, outcome, covariates) {
  vars_needed <- c(exposure, mediator, outcome, covariates)
  df <- data[, vars_needed]
  df <- df[complete.cases(df), ]
  
  if (nrow(df) < 30) {
    return(NULL)
  }
  
  colnames(df)[colnames(df) == exposure] <- "X"
  colnames(df)[colnames(df) == mediator] <- "M"
  colnames(df)[colnames(df) == outcome]  <- "Y"
  
  med_model <- lm(M ~ X + ., data = df[, c("M", "X", covariates)])
  out_model <- lm(Y ~ X + M + ., data = df[, c("Y", "X", "M", covariates)])
  
  environment(med_model$terms) <- environment()
  environment(out_model$terms) <- environment()
  
  med_res <- mediate(
    med_model,
    out_model,
    treat = "X",
    mediator = "M",
    boot = TRUE,
    sims = 1000
  )
  
  tibble(
    exposure = exposure,
    mediator = mediator,
    outcome  = outcome,
    N = nrow(df),
    
    ACME  = med_res$d0,
    ACME_low = med_res$d0.ci[1],
    ACME_high = med_res$d0.ci[2],
    ACME_p = med_res$d0.p,
    
    ADE  = med_res$z0,
    ADE_low = med_res$z0.ci[1],
    ADE_high = med_res$z0.ci[2],
    ADE_p = med_res$z0.p,
    
    Total  = med_res$tau.coef,
    Total_low = med_res$tau.ci[1],
    Total_high = med_res$tau.ci[2],
    Total_p = med_res$tau.p
  )
}

results_forward <- map_dfr(exposures_main, function(exp) {
  map_dfr(outcomes, function(out) {
    run_mediation(
      data = cross_data,
      exposure = exp,
      mediator = mediator_main,
      outcome  = out,
      covariates = covariates_main
    )
  })
}) %>%
  mutate(direction = "SpO2 -> ULK1 -> Cognition")

results_reverse <- map_dfr(mediators_reverse, function(med) {
  map_dfr(outcomes, function(out) {
    run_mediation(
      data = cross_data,
      exposure = exposure_reverse,
      mediator = med,
      outcome  = out,
      covariates = covariates_main
    )
  })
}) %>%
  mutate(direction = "ULK1 -> SpO2 -> Cognition")

results_med_all <- bind_rows(results_forward, results_reverse)

if (nrow(results_med_all) == 0) {
  stop("No mediation results generated.")
}

results_med_all <- results_med_all %>%
  mutate(
    exposure_label = get_exposure_label(exposure),
    mediator_label = get_exposure_label(mediator),
    outcome_label = get_outcome_label(outcome)
  )

results_med_export <- results_med_all %>%
  mutate(
    ACME_ci  = fmt_ci(ACME, ACME_low, ACME_high),
    ADE_ci   = fmt_ci(ADE, ADE_low, ADE_high),
    Total_ci = fmt_ci(Total, Total_low, Total_high),
    p_ACME  = ifelse(ACME_p  < 0.001, "<0.001", sprintf("%.3f", ACME_p)),
    p_ADE   = ifelse(ADE_p   < 0.001, "<0.001", sprintf("%.3f", ADE_p)),
    p_Total = ifelse(Total_p < 0.001, "<0.001", sprintf("%.3f", Total_p))
  )

results_med_export <- results_med_export[, c(
  "direction",
  "exposure", "exposure_label",
  "mediator", "mediator_label",
  "outcome", "outcome_label",
  "N",
  "ACME_ci", "p_ACME",
  "ADE_ci", "p_ADE",
  "Total_ci", "p_Total"
)]

############################################################
## 7. Interaction analyses
############################################################

covariates_int <- c(
  "age", "sex", "education",
  "obesity", "smoke_exposure",
  "diabetes", "CVD", "hypertension",
  "TST", "sleep_efficiency", "cpap_use",
  "AHI_bin5"
)

required_vars_int <- unique(c(exposures, outcomes, covariates_int))
missing_vars_int <- setdiff(required_vars_int, names(cross_data))

if (length(missing_vars_int) > 0) {
  stop(
    paste0(
      "Missing required interaction variables in `cross_data`: ",
      paste(missing_vars_int, collapse = ", ")
    )
  )
}

exp_pairs <- combn(exposures, 2, simplify = FALSE)

run_interaction <- function(data, exposure1, exposure2, outcome, covariates) {
  vars_needed <- c(exposure1, exposure2, outcome, covariates)
  df <- data[, vars_needed]
  df <- df[complete.cases(df), ]
  
  if (nrow(df) < 30) {
    return(NULL)
  }
  
  formula_str <- paste0(
    outcome, " ~ ",
    exposure1, " * ", exposure2, " + ",
    paste(covariates, collapse = " + ")
  )
  
  model <- lm(as.formula(formula_str), data = df)
  coef_name <- paste0(exposure1, ":", exposure2)
  coef_table <- summary(model)$coefficients
  
  if (!(coef_name %in% rownames(coef_table))) {
    coef_name <- paste0(exposure2, ":", exposure1)
  }
  
  if (!(coef_name %in% rownames(coef_table))) {
    return(NULL)
  }
  
  tibble(
    exposure1 = exposure1,
    exposure2 = exposure2,
    outcome   = outcome,
    N         = nrow(df),
    beta      = coef_table[coef_name, "Estimate"],
    se        = coef_table[coef_name, "Std. Error"],
    p         = coef_table[coef_name, "Pr(>|t|)"]
  )
}

results_interaction <- map_dfr(exp_pairs, function(pair) {
  map_dfr(outcomes, function(out) {
    run_interaction(
      data = cross_data,
      exposure1 = pair[1],
      exposure2 = pair[2],
      outcome   = out,
      covariates = covariates_int
    )
  })
})

if (nrow(results_interaction) == 0) {
  stop("No interaction results generated.")
}

results_interaction <- results_interaction %>%
  mutate(
    exposure1_label = get_exposure_label(exposure1),
    exposure2_label = get_exposure_label(exposure2),
    outcome_label = get_outcome_label(outcome),
    interaction = paste0(exposure1_label, " x ", exposure2_label),
    beta_ci = sprintf("%.2f (%.2f)", beta, se),
    p_text = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

results_export <- results_interaction[, c(
  "interaction",
  "outcome",
  "outcome_label",
  "N",
  "beta",
  "se",
  "p",
  "beta_ci",
  "p_text"
)]

results_export_clean <- results_export %>%
  mutate(across(everything(), force_utf8)) %>%
  mutate(across(where(is.character), ~ gsub("[^ -~]", "", .)))

############################################################
## 8. Export selected outputs
############################################################

write_xlsx(
  list(
    Table1 = table1,
    Table2 = table2,
    Cross_sectional_all_results = results_cross,
    Model1_to_Model4 = results_cross_main,
    Model5_AHIbin5 = results_cross_model5,
    Mediation_results = results_med_export
  ),
  "ASAP_cross_sectional_outputs.xlsx"
)

wb <- createWorkbook()
addWorksheet(wb, "interaction_results")
writeData(wb, sheet = 1, x = results_export_clean)
setColWidths(wb, sheet = 1, cols = 1:ncol(results_export_clean), widths = "auto")
saveWorkbook(
  wb,
  file = "ASAP_interaction_results.xlsx",
  overwrite = TRUE
)

############################################################
## 9. Export figures
############################################################

ggsave(
  filename = "Supplementary_Fig_S1_ULK1_distribution.png",
  plot = Figure_S1,
  width = 180,
  height = 90,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = "Supplementary_Fig_S1_ULK1_distribution.tiff",
  plot = Figure_S1,
  width = 180,
  height = 90,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Fig_1_cross_sectional_forest.png",
  plot = Figure_cross_forest,
  width = 250,
  height = 210,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = "Fig_1_cross_sectional_forest.tiff",
  plot = Figure_cross_forest,
  width = 250,
  height = 210,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Supplementary_Fig_S3_cross_sectional_diagnostics.png",
  plot = Figure_S3_diag,
  width = 280,
  height = 370,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Supplementary_Fig_S3_cross_sectional_diagnostics.tiff",
  plot = Figure_S3_diag,
  width = 280,
  height = 370,
  units = "mm",
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw",
  bg = "white"
)

message("Cross-sectional analyses completed successfully.")