############################################################
## 03_survival_analysis.R
##
## Purpose:
## Survival analyses for the ASAP manuscript.
##
## Expected input objects from 01_prepare_analysis_dataset.R:
## - surv_data_all
## - surv_data_ulk1
##
## Main outputs:
## - continuous_results
## - joint_results
## - cox_ph_results
##
## Notes:
## This script assumes that analysis-ready datasets have already
## been created in memory by running 01_prepare_analysis_dataset.R.
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260304)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(writexl)
  library(stringr)
  library(scales)
  library(grid)
  library(ragg)
})

############################################################
## Input checks
############################################################

if (!exists("surv_data_all")) {
  stop("Object `surv_data_all` not found. Please run 01_prepare_analysis_dataset.R first.")
}

if (!exists("surv_data_ulk1")) {
  stop("Object `surv_data_ulk1` not found. Please run 01_prepare_analysis_dataset.R first.")
}

############################################################
## Helper functions
############################################################

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

to_bin01 <- function(x) {
  xn <- suppressWarnings(as.numeric(as.character(x)))
  dplyr::case_when(
    is.na(xn) ~ NA_integer_,
    xn %in% c(1) ~ 1L,
    xn %in% c(0, 2) ~ 0L,
    TRUE ~ NA_integer_
  )
}

safe_drop_na <- function(df, vars) {
  vars <- intersect(vars, names(df))
  tidyr::drop_na(df, dplyr::all_of(vars))
}

sanitize_for_xlsx <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  for (nm in names(df)) {
    v <- df[[nm]]
    
    if (is.list(v)) {
      df[[nm]] <- vapply(v, function(x) paste(x, collapse = ";"), character(1))
      next
    }
    
    if (is.numeric(v)) {
      v[!is.finite(v)] <- NA_real_
      df[[nm]] <- v
      next
    }
    
    if (is.integer(v)) {
      df[[nm]] <- v
      next
    }
    
    vv <- as.character(v)
    vv <- enc2utf8(vv)
    vv <- gsub("[\r\n\t]", " ", vv)
    df[[nm]] <- vv
  }
  
  df
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

############################################################
## ASAP plotting system
############################################################

ASAP_colors <- list(
  SpO2_fill = "#B7E1CD",
  SpO2_line = "#2E8B57",
  ULK1_fill = "#F2B3AA",
  ULK1_line = "#C54A3F"
)

theme_ASAP_surv <- theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "#333333"),
    axis.line = element_line(color = "#333333"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, color = "#333333"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_blank()
  )

theme_ASAP_forest_surv <- theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "#333333"),
    axis.line = element_line(color = "#333333"),
    legend.position = "none"
  )

theme_ASAP_diag <- ggplot2::theme_classic(base_family = "sans", base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = ggplot2::element_text(size = 9),
    axis.text = ggplot2::element_text(size = 9, color = "#333333"),
    axis.line = ggplot2::element_line(color = "#333333"),
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4)
  )

############################################################
## 1. Rebuild survival analysis datasets
############################################################

surv_all <- surv_data_all
surv_ulk1 <- surv_data_ulk1

surv_all <- surv_all %>%
  mutate(
    stime_years    = to_num(stime_years),
    MACEPLUS       = to_bin01(MACEPLUS),
    SpO2           = to_num(SpO2),
    age            = to_num(age),
    sex            = to_num(sex),
    obesity        = to_bin01(obesity),
    smoke_exposure = to_bin01(smoke_exposure),
    diabetes       = to_bin01(diabetes),
    CVD            = to_bin01(CVD),
    hypertension   = to_bin01(hypertension),
    TST            = to_num(TST),
    high_LDL       = to_bin01(high_LDL),
    cpap_use       = to_bin01(cpap_use),
    AHI_bin5       = as.factor(AHI_bin5)
  )

surv_ulk1 <- surv_ulk1 %>%
  mutate(
    stime_years    = to_num(stime_years),
    MACEPLUS       = to_bin01(MACEPLUS),
    SpO2           = to_num(SpO2),
    lnULK1         = to_num(lnULK1),
    age            = to_num(age),
    sex            = to_num(sex),
    obesity        = to_bin01(obesity),
    smoke_exposure = to_bin01(smoke_exposure),
    diabetes       = to_bin01(diabetes),
    CVD            = to_bin01(CVD),
    hypertension   = to_bin01(hypertension),
    TST            = to_num(TST),
    high_LDL       = to_bin01(high_LDL),
    cpap_use       = to_bin01(cpap_use),
    AHI_bin5       = as.factor(AHI_bin5)
  )

############################################################
## 2. Grouped exposure variables for KM and joint Cox
############################################################

cut_spo2  <- 94.7
cut_lnulk <- 3.98

surv_all <- surv_all %>%
  mutate(
    SpO2_grp = case_when(
      !is.na(SpO2) & SpO2 >= cut_spo2 ~ "High SpO2",
      !is.na(SpO2) & SpO2 <  cut_spo2 ~ "Low SpO2",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    SpO2_grp = factor(SpO2_grp, levels = c("High SpO2", "Low SpO2"))
  )

surv_ulk1 <- surv_ulk1 %>%
  mutate(
    SpO2_grp = case_when(
      !is.na(SpO2) & SpO2 >= cut_spo2 ~ "High SpO2",
      !is.na(SpO2) & SpO2 <  cut_spo2 ~ "Low SpO2",
      TRUE ~ NA_character_
    ),
    ULK1_grp = case_when(
      !is.na(lnULK1) & lnULK1 >= cut_lnulk ~ "High ULK1",
      !is.na(lnULK1) & lnULK1 <  cut_lnulk ~ "Low ULK1",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    SpO2_grp = factor(SpO2_grp, levels = c("High SpO2", "Low SpO2")),
    ULK1_grp = factor(ULK1_grp, levels = c("High ULK1", "Low ULK1"))
  ) %>%
  mutate(
    SpO2_ULK1_grp = case_when(
      SpO2_grp == "High SpO2" & ULK1_grp == "High ULK1" ~ "High SpO2 + High ULK1",
      SpO2_grp == "Low SpO2"  & ULK1_grp == "High ULK1" ~ "Low SpO2 + High ULK1",
      SpO2_grp == "High SpO2" & ULK1_grp == "Low ULK1"  ~ "High SpO2 + Low ULK1",
      SpO2_grp == "Low SpO2"  & ULK1_grp == "Low ULK1"  ~ "Low SpO2 + Low ULK1",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    SpO2_ULK1_grp = relevel(
      factor(
        SpO2_ULK1_grp,
        levels = c(
          "High SpO2 + High ULK1",
          "Low SpO2 + High ULK1",
          "High SpO2 + Low ULK1",
          "Low SpO2 + Low ULK1"
        )
      ),
      ref = "High SpO2 + High ULK1"
    )
  )

############################################################
## 3. Kaplan–Meier curves (Fig. 2)
############################################################

km_spo2_data <- surv_all %>%
  filter(!is.na(stime_years), !is.na(MACEPLUS), !is.na(SpO2_grp))

km_ulk1_data <- surv_ulk1 %>%
  filter(!is.na(stime_years), !is.na(MACEPLUS), !is.na(ULK1_grp))

km_spo2_ulk1_data <- surv_ulk1 %>%
  filter(!is.na(stime_years), !is.na(MACEPLUS), !is.na(SpO2_ULK1_grp))

fit_spo2 <- survfit(Surv(stime_years, MACEPLUS) ~ SpO2_grp, data = km_spo2_data)
fit_ulk1 <- survfit(Surv(stime_years, MACEPLUS) ~ ULK1_grp, data = km_ulk1_data)
fit_spo2_ulk1 <- survfit(Surv(stime_years, MACEPLUS) ~ SpO2_ULK1_grp, data = km_spo2_ulk1_data)

x_breaks_surv <- c(0, 5, 10, 15)
y_breaks_surv <- c(0.0, 0.1, 0.2, 0.3)
y_labels_surv <- c("0.0", "0.1", "0.2", "0.3")
y_upper_surv  <- 0.33

combo_colors <- list(
  SpO2_ULK1_high_high = "#ecb06f",
  SpO2_ULK1_low_high  = ASAP_colors$ULK1_line,
  SpO2_ULK1_high_low  = ASAP_colors$SpO2_line,
  SpO2_ULK1_low_low   = "#A6A6A6"
)

plot_A_core <- ggsurvplot(
  fit = fit_spo2,
  data = km_spo2_data,
  fun = "event",
  conf.int = FALSE,
  risk.table = FALSE,
  censor = TRUE,
  size = 0.7,
  palette = c(ASAP_colors$SpO2_line, ASAP_colors$SpO2_fill),
  legend.title = "",
  legend.labs = c("Group 1: High SpO2", "Group 2: Low SpO2"),
  ggtheme = theme_ASAP_surv
)$plot +
  scale_x_continuous(limits = c(0, 15), breaks = x_breaks_surv, expand = expansion(mult = c(0.04, 0.02))) +
  scale_y_continuous(limits = c(0, y_upper_surv), breaks = y_breaks_surv, labels = y_labels_surv, expand = expansion(mult = c(0.04, 0.02))) +
  labs(title = NULL, x = "Years", y = "Probability of MACE") +
  theme_ASAP_surv +
  theme(legend.position = c(0.04, 0.96), legend.justification = c(0, 1))

plot_B_core <- ggsurvplot(
  fit = fit_ulk1,
  data = km_ulk1_data,
  fun = "event",
  conf.int = FALSE,
  risk.table = FALSE,
  censor = TRUE,
  size = 0.7,
  palette = c(ASAP_colors$ULK1_line, ASAP_colors$ULK1_fill),
  legend.title = "",
  legend.labs = c("Group 1: High ULK1", "Group 2: Low ULK1"),
  ggtheme = theme_ASAP_surv
)$plot +
  scale_x_continuous(limits = c(0, 15), breaks = x_breaks_surv, expand = expansion(mult = c(0.04, 0.02))) +
  scale_y_continuous(limits = c(0, y_upper_surv), breaks = y_breaks_surv, labels = y_labels_surv, expand = expansion(mult = c(0.04, 0.02))) +
  labs(title = NULL, x = "Years", y = "Probability of MACE") +
  theme_ASAP_surv +
  theme(legend.position = c(0.04, 0.96), legend.justification = c(0, 1))

plot_C_core <- ggsurvplot(
  fit = fit_spo2_ulk1,
  data = km_spo2_ulk1_data,
  fun = "event",
  conf.int = FALSE,
  risk.table = FALSE,
  censor = TRUE,
  size = 0.7,
  palette = c(
    combo_colors$SpO2_ULK1_high_high,
    combo_colors$SpO2_ULK1_low_high,
    combo_colors$SpO2_ULK1_high_low,
    combo_colors$SpO2_ULK1_low_low
  ),
  legend.title = "",
  legend.labs = c(
    "Group 1: High SpO2 + High ULK1",
    "Group 2: Low SpO2 + High ULK1",
    "Group 3: High SpO2 + Low ULK1",
    "Group 4: Low SpO2 + Low ULK1"
  ),
  ggtheme = theme_ASAP_surv
)$plot +
  scale_x_continuous(limits = c(0, 15), breaks = x_breaks_surv, expand = expansion(mult = c(0.04, 0.02))) +
  scale_y_continuous(limits = c(0, y_upper_surv), breaks = y_breaks_surv, labels = y_labels_surv, expand = expansion(mult = c(0.04, 0.02))) +
  labs(title = NULL, x = "Years", y = "Probability of MACE") +
  theme_ASAP_surv +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.direction = "vertical")

add_panel_header <- function(plot_obj, panel_id, panel_title) {
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
  
  top_strip <- tag_grob + title_grob + plot_layout(widths = c(0.08, 0.92))
  
  top_strip / plot_obj +
    plot_layout(heights = c(0.10, 0.90)) &
    theme(plot.margin = margin(6, 6, 6, 6))
}

plot_A <- add_panel_header(plot_A_core, "A", "SpO2 (median split)")
plot_B <- add_panel_header(plot_B_core, "B", "ln(ULK1) (median split)")
plot_C <- add_panel_header(plot_C_core, "C", "SpO2 + ULK1")

row_top <- wrap_plots(plot_A, plot_B, ncol = 2)
row_bottom <- wrap_plots(plot_spacer(), plot_C, plot_spacer(), ncol = 3, widths = c(0.12, 1, 0.12))

Figure_2_KM <- wrap_plots(row_top, row_bottom, ncol = 1, heights = c(1, 1.05))

############################################################
## 4. Cox analyses for MACE
############################################################

required_vars_all <- c(
  "stime_years", "MACEPLUS", "SpO2",
  "age", "sex", "obesity", "smoke_exposure",
  "diabetes", "CVD", "hypertension",
  "TST", "high_LDL", "cpap_use", "AHI_bin5"
)

required_vars_ulk1 <- c(
  "stime_years", "MACEPLUS", "SpO2", "lnULK1", "SpO2_ULK1_grp",
  "age", "sex", "obesity", "smoke_exposure",
  "diabetes", "CVD", "hypertension",
  "TST", "high_LDL", "cpap_use", "AHI_bin5"
)

missing_all <- setdiff(required_vars_all, names(surv_all))
missing_ulk1 <- setdiff(required_vars_ulk1, names(surv_ulk1))

if (length(missing_all) > 0) {
  stop(paste0("Missing required variables in `surv_all`: ", paste(missing_all, collapse = ", ")))
}

if (length(missing_ulk1) > 0) {
  stop(paste0("Missing required variables in `surv_ulk1`: ", paste(missing_ulk1, collapse = ", ")))
}

model_covs <- list(
  Model1 = c(),
  Model2 = c("age", "sex"),
  Model3 = c("age", "sex", "obesity", "smoke_exposure"),
  Model4 = c("age", "sex", "obesity", "smoke_exposure", "diabetes", "CVD", "hypertension"),
  Model5 = c("age", "sex", "obesity", "smoke_exposure", "diabetes", "CVD", "hypertension", "TST", "high_LDL", "cpap_use", "AHI_bin5")
)

continuous_exposures <- tibble::tribble(
  ~dataset,    ~exposure, ~exposure_label,
  "surv_all",  "SpO2",    "SpO2",
  "surv_ulk1", "lnULK1",  "ln(ULK1)"
)

joint_exposures <- tibble::tribble(
  ~dataset,     ~joint_var,       ~joint_label,  ~reference_group,
  "surv_ulk1",  "SpO2_ULK1_grp",  "SpO2 + ULK1", "High SpO2 + High ULK1"
)

run_cox_continuous <- function(df, exposure, exposure_label, covars, model_name, dataset_name, min_n = 30) {
  covars_use <- intersect(covars, names(df))
  needed <- unique(c("stime_years", "MACEPLUS", exposure, covars_use))
  df_m <- safe_drop_na(df, needed)
  
  if (nrow(df_m) < min_n) return(NULL)
  if (length(unique(df_m[[exposure]])) < 2) return(NULL)
  
  f <- as.formula(
    paste0("Surv(stime_years, MACEPLUS) ~ ", paste(c(exposure, covars_use), collapse = " + "))
  )
  
  fit <- try(coxph(f, data = df_m), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure) %>%
    mutate(
      analysis_type = "Continuous",
      dataset = dataset_name,
      exposure = exposure,
      exposure_label = exposure_label,
      model = model_name,
      n = fit$n,
      events = fit$nevent,
      covariates = ifelse(length(covars_use) == 0, "Unadjusted", paste(covars_use, collapse = ";"))
    ) %>%
    transmute(
      analysis_type = analysis_type,
      dataset = dataset,
      exposure = exposure,
      exposure_label = exposure_label,
      model = model,
      term = term,
      n = n,
      events = events,
      covariates = covariates,
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      SE = std.error,
      z = statistic,
      p = p.value
    )
}

run_cox_joint <- function(df, joint_var, joint_label, reference_group, covars, model_name, dataset_name, min_n = 30) {
  covars_use <- intersect(covars, names(df))
  needed <- unique(c("stime_years", "MACEPLUS", joint_var, covars_use))
  df_m <- safe_drop_na(df, needed)
  
  if (nrow(df_m) < min_n) return(NULL)
  if (length(unique(df_m[[joint_var]])) < 2) return(NULL)
  
  f <- as.formula(
    paste0("Surv(stime_years, MACEPLUS) ~ ", paste(c(joint_var, covars_use), collapse = " + "))
  )
  
  fit <- try(coxph(f, data = df_m), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(startsWith(term, joint_var)) %>%
    mutate(
      analysis_type = "Joint_group",
      dataset = dataset_name,
      joint_var = joint_var,
      joint_label = joint_label,
      reference_group = reference_group,
      model = model_name,
      n = fit$n,
      events = fit$nevent,
      covariates = ifelse(length(covars_use) == 0, "Unadjusted", paste(covars_use, collapse = ";"))
    ) %>%
    mutate(comparison_group = gsub(paste0("^", joint_var), "", term)) %>%
    transmute(
      analysis_type = analysis_type,
      dataset = dataset,
      joint_var = joint_var,
      joint_label = joint_label,
      reference_group = reference_group,
      comparison_group = comparison_group,
      model = model,
      term = term,
      n = n,
      events = events,
      covariates = covariates,
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      SE = std.error,
      z = statistic,
      p = p.value
    )
}

continuous_results_list <- list()

for (i in seq_len(nrow(continuous_exposures))) {
  dataset_name <- continuous_exposures$dataset[i]
  exposure_i <- continuous_exposures$exposure[i]
  exposure_label_i <- continuous_exposures$exposure_label[i]
  
  df_i <- if (dataset_name == "surv_all") surv_all else surv_ulk1
  
  for (mname in names(model_covs)) {
    res <- run_cox_continuous(
      df = df_i,
      exposure = exposure_i,
      exposure_label = exposure_label_i,
      covars = model_covs[[mname]],
      model_name = mname,
      dataset_name = dataset_name,
      min_n = 30
    )
    
    if (!is.null(res) && nrow(res) > 0) {
      continuous_results_list[[length(continuous_results_list) + 1]] <- res
    }
  }
}

continuous_results <- bind_rows(continuous_results_list) %>%
  arrange(exposure, model)

joint_results_list <- list()

for (i in seq_len(nrow(joint_exposures))) {
  dataset_name <- joint_exposures$dataset[i]
  joint_var_i <- joint_exposures$joint_var[i]
  joint_label_i <- joint_exposures$joint_label[i]
  reference_i <- joint_exposures$reference_group[i]
  
  df_i <- if (dataset_name == "surv_all") surv_all else surv_ulk1
  
  for (mname in names(model_covs)) {
    res <- run_cox_joint(
      df = df_i,
      joint_var = joint_var_i,
      joint_label = joint_label_i,
      reference_group = reference_i,
      covars = model_covs[[mname]],
      model_name = mname,
      dataset_name = dataset_name,
      min_n = 30
    )
    
    if (!is.null(res) && nrow(res) > 0) {
      joint_results_list[[length(joint_results_list) + 1]] <- res
    }
  }
}

joint_results <- bind_rows(joint_results_list) %>%
  arrange(joint_var, model, comparison_group)

continuous_results <- continuous_results %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    p_text = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

joint_results <- joint_results %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    p_text = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

############################################################
## 5. Fig. 3 and Fig. 4 survival forest plots
############################################################

forest_surv_cont <- continuous_results %>%
  mutate(
    model = factor(model, levels = c("Model1", "Model2", "Model3", "Model4", "Model5")),
    exposure_label = factor(exposure_label, levels = c("SpO2", "ln(ULK1)"))
  ) %>%
  arrange(exposure_label, model) %>%
  mutate(
    hr_ci = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    p_text = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
    sig = p < 0.05
  )

forest_surv_joint <- joint_results %>%
  mutate(
    model = factor(model, levels = c("Model1", "Model2", "Model3", "Model4", "Model5")),
    joint_label = factor(joint_label, levels = c("SpO2 + ULK1")),
    comparison_group = case_when(
      comparison_group == "Low SpO2 + High ULK1"  ~ "Group 2: Low SpO2 + High ULK1",
      comparison_group == "High SpO2 + Low ULK1"  ~ "Group 3: High SpO2 + Low ULK1",
      comparison_group == "Low SpO2 + Low ULK1"   ~ "Group 4: Low SpO2 + Low ULK1",
      TRUE ~ comparison_group
    )
  ) %>%
  arrange(joint_label, comparison_group, model) %>%
  mutate(
    hr_ci = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    p_text = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
    sig = p < 0.05
  )

make_y_positions_equal <- function(n, start = NULL, step = 1.2) {
  if (is.null(start)) {
    start <- (n - 1) * step
  }
  seq(from = start, by = -step, length.out = n)
}

plot_surv_cont_row <- function(data, exposure_name, color, tag_letter) {
  df <- data %>%
    filter(exposure_label == exposure_name) %>%
    arrange(model)
  
  df$y <- make_y_positions_equal(nrow(df), step = 1.15)
  
  y_top <- max(df$y) + 0.9
  y_bottom <- min(df$y) - 0.7
  
  max_hr <- max(df$CI_high, na.rm = TRUE)
  min_hr <- min(df$CI_low, na.rm = TRUE)
  
  x_max <- max(2.1, ceiling(max_hr * 10) / 10 + 0.15)
  x_min <- min(0.5, floor(min_hr * 10) / 10 - 0.10)
  x_min <- max(0.35, x_min)
  
  tag_panel <- wrap_elements(
    full = textGrob(
      tag_letter,
      x = 0.5, y = 0.92,
      gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "sans")
    )
  )
  
  exposure_panel <- ggplot(data.frame(y = mean(df$y)), aes(y = y)) +
    geom_text(
      aes(x = 0, label = exposure_name),
      hjust = 0.5,
      family = "sans",
      size = 3.3
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Exposure",
      hjust = 0.5,
      family = "sans",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.8, 0.8), clip = "off") +
    theme_void()
  
  model_panel <- ggplot(df, aes(y = y)) +
    geom_text(
      aes(x = 0, label = model),
      hjust = 0.5,
      family = "sans",
      size = 3
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Model",
      hjust = 0.5,
      family = "sans",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.7, 0.7), clip = "off") +
    theme_void()
  
  forest_panel <- ggplot(df, aes(y = y, x = HR, xmin = CI_low, xmax = CI_high)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_point(
      data = df %>% filter(sig),
      aes(x = HR, y = y),
      inherit.aes = FALSE,
      size = 4.2,
      color = alpha(color, 0.18)
    ) +
    geom_errorbarh(height = 0.10, linewidth = 0.5, color = color) +
    geom_point(size = 1.6, color = color) +
    scale_x_continuous(limits = c(x_min, x_max), breaks = pretty(c(x_min, x_max), n = 5)) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    labs(title = exposure_name, x = "Hazard ratio (95% CI)", y = NULL) +
    theme_ASAP_forest_surv +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  hr_panel <- ggplot(df, aes(y = y)) +
    geom_text(
      aes(x = 0, label = hr_ci, fontface = ifelse(sig, "bold", "plain")),
      hjust = 0.5,
      family = "Courier",
      size = 2.8
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "HR (95% CI)",
      hjust = 0.5,
      family = "sans",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
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
      family = "sans",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.45, 0.45), clip = "off") +
    theme_void()
  
  tag_panel + exposure_panel + model_panel + forest_panel + hr_panel + p_panel +
    plot_layout(widths = c(0.35, 1.10, 1.05, 4.35, 2.00, 0.85))
}

plot_surv_joint_row <- function(data, joint_name, color_map, tag_letter) {
  df <- data %>%
    filter(joint_label == joint_name) %>%
    arrange(comparison_group, model)
  
  df$comparison_group <- factor(df$comparison_group, levels = unique(as.character(df$comparison_group)))
  comp_levels <- levels(df$comparison_group)
  model_step <- 1.05
  block_gap <- 1.15
  
  y_vals <- numeric(nrow(df))
  current_top <- (length(comp_levels) - 1) * (5 * model_step + block_gap)
  
  for (i in seq_along(comp_levels)) {
    idx <- which(df$comparison_group == comp_levels[i])
    y_vals[idx] <- seq(from = current_top + 4 * model_step, by = -model_step, length.out = 5)
    current_top <- current_top - (5 * model_step + block_gap)
  }
  
  df$y <- y_vals
  
  comp_df <- df %>%
    group_by(comparison_group) %>%
    summarise(
      y_center = mean(y),
      comp_display = str_wrap(first(as.character(comparison_group)), width = 22),
      .groups = "drop"
    )
  
  y_top <- max(df$y) + 1.2
  y_bottom <- min(df$y) - 0.7
  
  max_hr <- max(df$CI_high, na.rm = TRUE)
  min_hr <- min(df$CI_low, na.rm = TRUE)
  
  x_max <- max(3.5, ceiling(max_hr * 5) / 5 + 0.2)
  x_min <- min(0.2, floor(min_hr * 10) / 10 - 0.1)
  x_min <- max(0, x_min)
  
  df <- df %>%
    mutate(plot_color = recode(as.character(comparison_group), !!!color_map))
  
  tag_panel <- wrap_elements(
    full = textGrob(
      tag_letter,
      x = 0.5, y = 0.98,
      gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "sans")
    )
  )
  
  comp_panel <- ggplot(comp_df, aes(y = y_center)) +
    geom_text(
      aes(x = 0, label = comp_display),
      hjust = 0.5,
      family = "sans",
      size = 3.0
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Comparison",
      hjust = 0.5,
      family = "sans",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-1.2, 1.2), clip = "off") +
    theme_void()
  
  model_panel <- ggplot(df, aes(y = y)) +
    geom_text(
      aes(x = 0, label = model),
      hjust = 0.5,
      family = "sans",
      size = 3
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Model",
      hjust = 0.5,
      family = "sans",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.75, 0.75), clip = "off") +
    theme_void()
  
  forest_panel <- ggplot(df, aes(y = y, x = HR, xmin = CI_low, xmax = CI_high)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_point(
      data = df %>% filter(sig),
      aes(x = HR, y = y, color = plot_color),
      inherit.aes = FALSE,
      size = 4.2,
      alpha = 0.18,
      show.legend = FALSE
    ) +
    geom_errorbarh(aes(color = plot_color), height = 0.10, linewidth = 0.5, show.legend = FALSE) +
    geom_point(aes(color = plot_color), size = 1.6, show.legend = FALSE) +
    scale_color_identity() +
    scale_x_continuous(limits = c(x_min, x_max), breaks = pretty(c(x_min, x_max), n = 5)) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    labs(title = NULL, x = "Hazard ratio (95% CI)", y = NULL) +
    theme_ASAP_forest_surv +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  hr_panel <- ggplot(df, aes(y = y)) +
    geom_text(
      aes(x = 0, label = hr_ci, fontface = ifelse(sig, "bold", "plain")),
      hjust = 0.5,
      family = "Courier",
      size = 2.7
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "HR (95% CI)",
      hjust = 0.5,
      family = "sans",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-1.2, 1.2), clip = "off") +
    theme_void()
  
  p_panel <- ggplot(df, aes(y = y)) +
    geom_text(
      aes(x = 0, label = p_text, fontface = ifelse(sig, "bold", "plain")),
      hjust = 0.5,
      family = "Courier",
      size = 2.7
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "p value",
      hjust = 0.5,
      family = "sans",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.45, 0.45), clip = "off") +
    theme_void()
  
  tag_panel + comp_panel + model_panel + forest_panel + hr_panel + p_panel +
    plot_layout(widths = c(0.35, 2.10, 1.05, 4.60, 2.20, 0.85))
}

row_A <- plot_surv_cont_row(forest_surv_cont, "SpO2", ASAP_colors$SpO2_line, "A")
row_B <- plot_surv_cont_row(forest_surv_cont, "ln(ULK1)", ASAP_colors$ULK1_line, "B")
Figure_survival_continuous_forest <- row_A / row_B

spo2_ulk1_colors <- c(
  "Group 2: Low SpO2 + High ULK1" = ASAP_colors$ULK1_line,
  "Group 3: High SpO2 + Low ULK1" = ASAP_colors$SpO2_line,
  "Group 4: Low SpO2 + Low ULK1"  = "#A6A6A6"
)

Figure_survival_joint_forest <- plot_surv_joint_row(
  forest_surv_joint,
  "SpO2 + ULK1",
  spo2_ulk1_colors,
  ""
)

############################################################
## 6. Cox PH diagnostics (Supplementary Fig. S4)
############################################################

model5_covs <- c(
  "age", "sex", "obesity", "smoke_exposure",
  "diabetes", "CVD", "hypertension",
  "TST", "high_LDL", "cpap_use", "AHI_bin5"
)

diag_models <- tibble::tribble(
  ~analysis_type, ~dataset,    ~exposure_var,    ~exposure_label, ~plot_color,
  "continuous",   "surv_all",  "SpO2",           "SpO2",          ASAP_colors$SpO2_line,
  "continuous",   "surv_ulk1", "lnULK1",         "ln(ULK1)",      ASAP_colors$ULK1_line,
  "joint",        "surv_ulk1", "SpO2_ULK1_grp",  "SpO2 + ULK1",   ASAP_colors$ULK1_line
)

fit_cox_model5 <- function(df, exposure_var, covars, min_n = 30) {
  covars_use <- intersect(covars, names(df))
  needed <- unique(c("stime_years", "MACEPLUS", exposure_var, covars_use))
  df_m <- safe_drop_na(df, needed)
  
  if (nrow(df_m) < min_n) return(NULL)
  if (length(unique(df_m[[exposure_var]])) < 2) return(NULL)
  
  f <- as.formula(
    paste0("Surv(stime_years, MACEPLUS) ~ ", paste(c(exposure_var, covars_use), collapse = " + "))
  )
  
  fit <- try(survival::coxph(f, data = df_m), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

extract_zph_table <- function(fit, model_label, analysis_type, exposure_var, exposure_label) {
  zph <- survival::cox.zph(fit, transform = "km")
  ztab <- as.data.frame(zph$table)
  ztab$term <- rownames(ztab)
  rownames(ztab) <- NULL
  names(ztab)[1:3] <- c("chisq", "df", "p")
  
  ztab %>%
    dplyr::mutate(
      model = model_label,
      analysis_type = analysis_type,
      exposure_var = exposure_var,
      exposure_label = exposure_label,
      n = fit$n,
      events = fit$nevent
    ) %>%
    dplyr::select(
      analysis_type, exposure_var, exposure_label, model,
      term, chisq, df, p, n, events
    )
}

extract_zph_plot_data <- function(zph_obj, terms_keep = NULL) {
  ymat <- as.data.frame(zph_obj$y)
  ymat$time_transformed <- zph_obj$x
  
  long_df <- ymat %>%
    tidyr::pivot_longer(
      cols = -time_transformed,
      names_to = "term",
      values_to = "scaled_schoenfeld"
    )
  
  if (!is.null(terms_keep)) {
    long_df <- long_df %>% dplyr::filter(term %in% terms_keep)
  }
  
  long_df
}

plot_cox_zph_core <- function(plot_df, line_color) {
  term_levels <- unique(plot_df$term)
  p_list <- vector("list", length(term_levels))
  
  for (i in seq_along(term_levels)) {
    term_i <- term_levels[i]
    
    df_i <- plot_df %>%
      dplyr::filter(term == term_i)
    
    p_list[[i]] <- ggplot2::ggplot(
      df_i,
      ggplot2::aes(x = time_transformed, y = scaled_schoenfeld)
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey60",
        linewidth = 0.4
      ) +
      ggplot2::geom_point(
        size = 1.4,
        alpha = 0.7,
        color = line_color
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        se = FALSE,
        color = line_color,
        linewidth = 0.8
      ) +
      ggplot2::labs(
        title = NULL,
        x = "Transformed time",
        y = "Scaled Schoenfeld residuals"
      ) +
      theme_ASAP_diag
  }
  
  patchwork::wrap_plots(p_list, ncol = min(2, length(p_list)))
}

build_cox_diag_panel <- function(plot_df, panel_id, panel_title, line_color) {
  core_plot <- plot_cox_zph_core(plot_df = plot_df, line_color = line_color)
  
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
  
  top_strip <- tag_grob + title_grob + plot_layout(widths = c(0.08, 0.92))
  core_grob <- patchwork::patchworkGrob(core_plot)
  
  border_grob <- rectGrob(
    x = 0.5, y = 0.5,
    width = 0.985, height = 0.985,
    just = "center",
    gp = gpar(fill = NA, col = "grey45", lwd = 1.1, lty = "dashed")
  )
  
  boxed_diag <- wrap_elements(full = grobTree(core_grob, border_grob))
  
  top_strip / boxed_diag +
    plot_layout(heights = c(0.10, 0.90)) &
    theme(plot.margin = margin(8, 8, 8, 8))
}

zph_results_list <- list()
plot_df_list <- list()

for (i in seq_len(nrow(diag_models))) {
  dataset_name   <- diag_models$dataset[i]
  analysis_type  <- diag_models$analysis_type[i]
  exposure_var_i <- diag_models$exposure_var[i]
  exposure_lab_i <- diag_models$exposure_label[i]
  color_i        <- diag_models$plot_color[i]
  
  df_i <- if (dataset_name == "surv_all") surv_all else surv_ulk1
  
  fit_i <- fit_cox_model5(
    df = df_i,
    exposure_var = exposure_var_i,
    covars = model5_covs,
    min_n = 30
  )
  
  if (is.null(fit_i)) next
  
  zph_tab_i <- extract_zph_table(
    fit = fit_i,
    model_label = "Model5",
    analysis_type = analysis_type,
    exposure_var = exposure_var_i,
    exposure_label = exposure_lab_i
  )
  
  zph_results_list[[length(zph_results_list) + 1]] <- zph_tab_i
  
  zph_i <- survival::cox.zph(fit_i, transform = "km")
  
  if (analysis_type == "continuous") {
    terms_keep_i <- exposure_var_i
  } else {
    terms_keep_i <- rownames(zph_i$table)
    terms_keep_i <- terms_keep_i[startsWith(terms_keep_i, exposure_var_i)]
  }
  
  plot_df_i <- extract_zph_plot_data(zph_obj = zph_i, terms_keep = terms_keep_i)
  
  if (nrow(plot_df_i) == 0) next
  
  plot_df_list[[exposure_var_i]] <- list(
    plot_df = plot_df_i,
    exposure_label = exposure_lab_i,
    analysis_type = analysis_type,
    plot_color = color_i
  )
}

cox_ph_results <- dplyr::bind_rows(zph_results_list) %>%
  dplyr::arrange(analysis_type, exposure_var, term)

required_plot_keys <- c("SpO2", "lnULK1", "SpO2_ULK1_grp")
missing_plot_keys <- setdiff(required_plot_keys, names(plot_df_list))

if (length(missing_plot_keys) > 0) {
  stop(
    paste0(
      "Missing combined diagnostic plot data for: ",
      paste(missing_plot_keys, collapse = ", ")
    )
  )
}

panel_A <- build_cox_diag_panel(
  plot_df = plot_df_list[["SpO2"]]$plot_df,
  panel_id = "A",
  panel_title = "SpO2",
  line_color = plot_df_list[["SpO2"]]$plot_color
)

panel_B <- build_cox_diag_panel(
  plot_df = plot_df_list[["lnULK1"]]$plot_df,
  panel_id = "B",
  panel_title = "ln(ULK1)",
  line_color = plot_df_list[["lnULK1"]]$plot_color
)

panel_C <- build_cox_diag_panel(
  plot_df = plot_df_list[["SpO2_ULK1_grp"]]$plot_df,
  panel_id = "C",
  panel_title = "SpO2 + ULK1",
  line_color = plot_df_list[["SpO2_ULK1_grp"]]$plot_color
)

row_top_diag <- wrap_plots(panel_A, panel_B, ncol = 2)
row_bottom_diag <- wrap_plots(plot_spacer(), panel_C, plot_spacer(), ncol = 3, widths = c(0.10, 1, 0.10))

Figure_Supp_Cox_PH <- wrap_plots(row_top_diag, row_bottom_diag, ncol = 1, heights = c(1, 1.10))

############################################################
## 7. Export selected outputs
############################################################

continuous_out <- sanitize_for_xlsx(continuous_results)
joint_out <- sanitize_for_xlsx(joint_results)
cox_ph_out <- sanitize_for_xlsx(cox_ph_results)

write_xlsx(
  list(
    Continuous_Cox = continuous_out,
    Joint_Group_Cox = joint_out,
    Cox_PH_Assumption_Test = cox_ph_out
  ),
  path = "ASAP_survival_outputs.xlsx"
)

############################################################
## 8. Export figures
############################################################

ggsave(
  filename = "Fig_2_Kaplan_Meier_survival_curves.png",
  plot = Figure_2_KM,
  width = 170,
  height = 180,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = "Fig_2_Kaplan_Meier_survival_curves.tiff",
  plot = Figure_2_KM,
  width = 170,
  height = 180,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Fig_3_continuous_survival_forest.png",
  plot = Figure_survival_continuous_forest,
  width = 200,
  height = 160,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = "Fig_3_continuous_survival_forest.tiff",
  plot = Figure_survival_continuous_forest,
  width = 200,
  height = 160,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Fig_4_joint_survival_forest.png",
  plot = Figure_survival_joint_forest,
  width = 200,
  height = 160,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = "Fig_4_joint_survival_forest.tiff",
  plot = Figure_survival_joint_forest,
  width = 200,
  height = 160,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Supplementary_Fig_S4_Cox_PH_diagnostics.png",
  plot = Figure_Supp_Cox_PH,
  width = 250,
  height = 245,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Supplementary_Fig_S4_Cox_PH_diagnostics.tiff",
  plot = Figure_Supp_Cox_PH,
  width = 250,
  height = 245,
  units = "mm",
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw",
  bg = "white"
)

message("Survival analyses completed successfully.")