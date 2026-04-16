############################################################
## 04_multiple_imputation_sensitivity.R
##
## Purpose:
## Multiple imputation sensitivity analyses for the ASAP manuscript.
##
## Expected input objects from 01_prepare_analysis_dataset.R:
## - cross_data_cog_master
## - surv_data_all
## - surv_data_ulk1
##
## Main outputs:
## - imp_cross
## - results_cross_mi
## - imp_surv_all
## - imp_surv_ulk1
## - continuous_results_mi
## - joint_results_mi
##
## Notes:
## - Only missing covariates are imputed.
## - Exposures and outcomes are included as predictors but are not imputed.
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260304)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(mice)
  library(survival)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(scales)
  library(grid)
  library(writexl)
})

############################################################
## Input checks
############################################################

required_objects <- c("cross_data_cog_master", "surv_data_all", "surv_data_ulk1")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "The following required objects are missing: ",
      paste(missing_objects, collapse = ", "),
      ". Please run 01_prepare_analysis_dataset.R first."
    )
  )
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
  ULK1_line = "#C54A3F"
)

theme_ASAP <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "#333333"),
    axis.line = element_line(color = "#333333"),
    legend.position = "none"
  )

theme_ASAP_forest_surv <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "#333333"),
    axis.line = element_line(color = "#333333"),
    legend.position = "none"
  )

############################################################
## 1. Cross-sectional MI sensitivity analysis
############################################################

cross_data <- cross_data_cog_master %>%
  filter(!is.na(SpO2) | !is.na(lnULK1))

exposures <- c("SpO2", "lnULK1")
outcomes <- c("RAVLT_tot", "RAVLT_learn", "Stroop_t", "Stroop_r")

required_vars_cross <- c(
  exposures,
  outcomes,
  "age",
  "sex",
  "education",
  "obesity",
  "smoke_exposure",
  "diabetes",
  "CVD",
  "hypertension",
  "TST",
  "sleep_efficiency",
  "cpap_use",
  "AHI_bin5"
)

missing_vars_cross <- setdiff(required_vars_cross, names(cross_data))
if (length(missing_vars_cross) > 0) {
  stop(
    paste0(
      "Missing required cross-sectional MI variables: ",
      paste(missing_vars_cross, collapse = ", ")
    )
  )
}

mi_cross_data <- cross_data %>%
  dplyr::select(all_of(required_vars_cross))

mi_missing_summary_cross <- data.frame(
  Variable = names(mi_cross_data),
  Missing_n = sapply(mi_cross_data, function(x) sum(is.na(x))),
  Missing_pct = round(sapply(mi_cross_data, function(x) mean(is.na(x)) * 100), 2)
) %>%
  arrange(desc(Missing_pct))

ini_cross <- mice(
  mi_cross_data,
  maxit = 0,
  printFlag = FALSE
)

meth_cross <- ini_cross$method
pred_cross <- ini_cross$predictorMatrix

meth_cross[] <- ""

## Continuous covariates
meth_cross["age"]              <- "pmm"
meth_cross["TST"]              <- "pmm"
meth_cross["sleep_efficiency"] <- "pmm"

## Binary covariates
binary_covs_cross <- c(
  "sex",
  "education",
  "obesity",
  "smoke_exposure",
  "diabetes",
  "CVD",
  "hypertension",
  "cpap_use",
  "AHI_bin5"
)

meth_cross[binary_covs_cross] <- "logreg"

## All variables can be predictors; diagonal excluded
pred_cross[,] <- 1
diag(pred_cross) <- 0

imp_cross <- mice(
  data = mi_cross_data,
  m = 20,
  maxit = 20,
  method = meth_cross,
  predictorMatrix = pred_cross,
  printFlag = TRUE,
  seed = 20260304
)

model_formulas_cross <- list(
  Model1 = "~ EXPOSURE",
  Model2 = "~ EXPOSURE + age + sex",
  Model3 = "~ EXPOSURE + age + sex + education + obesity + smoke_exposure",
  Model4 = "~ EXPOSURE + age + sex + education + obesity + smoke_exposure + diabetes + CVD + hypertension",
  Model5 = "~ EXPOSURE + age + sex + education + obesity + smoke_exposure + diabetes + CVD + hypertension + TST + sleep_efficiency + cpap_use + AHI_bin5"
)

run_cross_model_mi <- function(exposure, outcome, formula_string, imp_obj) {
  formula_text <- paste0(
    outcome,
    " ",
    gsub("EXPOSURE", exposure, formula_string)
  )
  
  f <- as.formula(formula_text)
  completed_list <- mice::complete(imp_obj, action = "all")
  
  fit_list <- lapply(completed_list, function(dat) {
    stats::lm(formula = f, data = dat)
  })
  
  fit_mira <- mice::as.mira(fit_list)
  pooled_fit <- mice::pool(fit_mira)
  
  pooled_res <- summary(
    pooled_fit,
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    dplyr::filter(term == exposure) %>%
    dplyr::mutate(
      exposure = exposure,
      outcome = outcome
    )
  
  dat1 <- completed_list[[1]]
  model_n <- dat1 %>%
    dplyr::select(all_of(all.vars(f))) %>%
    stats::na.omit() %>%
    nrow()
  
  pooled_res %>%
    dplyr::mutate(N = model_n)
}

results_cross_mi_list <- list()

for (outcome in outcomes) {
  for (exposure in exposures) {
    for (model_name in names(model_formulas_cross)) {
      res <- run_cross_model_mi(
        exposure = exposure,
        outcome = outcome,
        formula_string = model_formulas_cross[[model_name]],
        imp_obj = imp_cross
      )
      
      res$model <- model_name
      results_cross_mi_list[[length(results_cross_mi_list) + 1]] <- res
    }
  }
}

results_cross_mi <- bind_rows(results_cross_mi_list) %>%
  transmute(
    exposure = exposure,
    outcome = outcome,
    model = model,
    N = N,
    Beta = estimate,
    SE = std.error,
    CI_low = `2.5 %`,
    CI_high = `97.5 %`,
    t = statistic,
    df = df,
    p = p.value
  ) %>%
  arrange(exposure, outcome, model) %>%
  mutate(
    exposure_label = get_exposure_label(exposure),
    outcome_label = get_outcome_label(outcome)
  )

############################################################
## 2. Cross-sectional MI forest plot
############################################################

forest_data_mi <- results_cross_mi %>%
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
    )
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

plot_forest_row_centered_mi <- function(data, exposure_name, color, tag_letter) {
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

row_A_mi <- plot_forest_row_centered_mi(
  forest_data_mi, "SpO2", ASAP_colors$SpO2_line, "A"
)

row_B_mi <- plot_forest_row_centered_mi(
  forest_data_mi, "lnULK1", ASAP_colors$ULK1_line, "B"
)

Figure_cross_forest_mi <- row_A_mi / row_B_mi

############################################################
## 3. Survival MI sensitivity analysis
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

mi_surv_all_vars <- c(
  "stime_years", "MACEPLUS",
  "SpO2",
  "age", "sex",
  "obesity", "smoke_exposure",
  "diabetes", "CVD", "hypertension",
  "TST", "high_LDL", "cpap_use",
  "AHI_bin5"
)

mi_surv_ulk1_vars <- c(
  "stime_years", "MACEPLUS",
  "SpO2", "lnULK1",
  "age", "sex",
  "obesity", "smoke_exposure",
  "diabetes", "CVD", "hypertension",
  "TST", "high_LDL", "cpap_use",
  "AHI_bin5"
)

mi_surv_all_data <- surv_all %>%
  dplyr::select(all_of(mi_surv_all_vars))

mi_surv_ulk1_data <- surv_ulk1 %>%
  dplyr::select(all_of(mi_surv_ulk1_vars))

mi_surv_all_missing <- data.frame(
  Variable = names(mi_surv_all_data),
  Missing_n = sapply(mi_surv_all_data, function(x) sum(is.na(x))),
  Missing_pct = round(sapply(mi_surv_all_data, function(x) mean(is.na(x)) * 100), 2)
) %>%
  arrange(desc(Missing_pct))

mi_surv_ulk1_missing <- data.frame(
  Variable = names(mi_surv_ulk1_data),
  Missing_n = sapply(mi_surv_ulk1_data, function(x) sum(is.na(x))),
  Missing_pct = round(sapply(mi_surv_ulk1_data, function(x) mean(is.na(x)) * 100), 2)
) %>%
  arrange(desc(Missing_pct))

ini_all <- mice(
  mi_surv_all_data,
  maxit = 0,
  printFlag = FALSE
)

ini_ulk1 <- mice(
  mi_surv_ulk1_data,
  maxit = 0,
  printFlag = FALSE
)

meth_all <- ini_all$method
pred_all <- ini_all$predictorMatrix

meth_ulk1 <- ini_ulk1$method
pred_ulk1 <- ini_ulk1$predictorMatrix

meth_all[] <- ""
meth_ulk1[] <- ""

## Continuous covariates
meth_all["age"]  <- "pmm"
meth_all["TST"]  <- "pmm"

meth_ulk1["age"] <- "pmm"
meth_ulk1["TST"] <- "pmm"

## Binary covariates
binary_covs_surv <- c(
  "sex",
  "obesity",
  "smoke_exposure",
  "diabetes",
  "CVD",
  "hypertension",
  "high_LDL",
  "cpap_use",
  "AHI_bin5"
)

meth_all[binary_covs_surv]  <- "logreg"
meth_ulk1[binary_covs_surv] <- "logreg"

pred_all[,] <- 1
diag(pred_all) <- 0

pred_ulk1[,] <- 1
diag(pred_ulk1) <- 0

imp_surv_all <- mice(
  data = mi_surv_all_data,
  m = 20,
  maxit = 20,
  method = meth_all,
  predictorMatrix = pred_all,
  printFlag = TRUE,
  seed = 20260304
)

imp_surv_ulk1 <- mice(
  data = mi_surv_ulk1_data,
  m = 20,
  maxit = 20,
  method = meth_ulk1,
  predictorMatrix = pred_ulk1,
  printFlag = TRUE,
  seed = 20260304
)

create_joint_groups <- function(df) {
  cut_spo2  <- 94.7
  cut_lnulk <- 3.98
  
  df %>%
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
}

model_covs_surv <- list(
  Model1 = c(),
  Model2 = c("age", "sex"),
  Model3 = c("age", "sex", "obesity", "smoke_exposure"),
  Model4 = c("age", "sex", "obesity", "smoke_exposure", "diabetes", "CVD", "hypertension"),
  Model5 = c("age", "sex", "obesity", "smoke_exposure", "diabetes", "CVD", "hypertension", "TST", "high_LDL", "cpap_use", "AHI_bin5")
)

continuous_exposures_surv <- tibble::tribble(
  ~imp_obj,         ~dataset_name, ~exposure, ~exposure_label,
  "imp_surv_all",   "surv_all",    "SpO2",    "SpO2",
  "imp_surv_ulk1",  "surv_ulk1",   "lnULK1",  "ln(ULK1)"
)

joint_exposures_surv <- tibble::tribble(
  ~imp_obj,         ~dataset_name, ~joint_var,       ~joint_label,  ~reference_group,
  "imp_surv_ulk1",  "surv_ulk1",   "SpO2_ULK1_grp",  "SpO2 + ULK1", "High SpO2 + High ULK1"
)

run_cox_continuous_mi <- function(imp_obj, exposure, exposure_label, covars, model_name, dataset_name, min_n = 30) {
  covars_use <- covars
  completed_list <- mice::complete(imp_obj, action = "all")
  
  formula_text <- paste0(
    "Surv(stime_years, MACEPLUS) ~ ",
    paste(c(exposure, covars_use), collapse = " + ")
  )
  f <- as.formula(formula_text)
  
  fit_list <- list()
  n_vec <- c()
  event_vec <- c()
  
  for (i in seq_along(completed_list)) {
    dat_i <- completed_list[[i]]
    needed <- unique(c("stime_years", "MACEPLUS", exposure, covars_use))
    dat_i <- safe_drop_na(dat_i, needed)
    
    if (nrow(dat_i) < min_n) next
    if (length(unique(dat_i[[exposure]])) < 2) next
    
    fit_i <- try(coxph(f, data = dat_i), silent = TRUE)
    if (inherits(fit_i, "try-error")) next
    
    fit_list[[length(fit_list) + 1]] <- fit_i
    n_vec <- c(n_vec, fit_i$n)
    event_vec <- c(event_vec, fit_i$nevent)
  }
  
  if (length(fit_list) == 0) return(NULL)
  
  fit_mira <- mice::as.mira(fit_list)
  pooled_fit <- mice::pool(fit_mira)
  
  summary(
    pooled_fit,
    conf.int = TRUE,
    conf.level = 0.95,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(term = as.character(term)) %>%
    dplyr::filter(term == exposure) %>%
    dplyr::mutate(
      analysis_type = "Continuous",
      dataset = dataset_name,
      exposure = exposure,
      exposure_label = exposure_label,
      model = model_name,
      n = round(mean(n_vec)),
      events = round(mean(event_vec)),
      covariates = ifelse(length(covars_use) == 0, "Unadjusted", paste(covars_use, collapse = ";"))
    ) %>%
    dplyr::transmute(
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
      CI_low = `2.5 %`,
      CI_high = `97.5 %`,
      SE = std.error,
      z = statistic,
      p = p.value
    )
}

run_cox_joint_mi <- function(imp_obj, joint_var, joint_label, reference_group, covars, model_name, dataset_name, min_n = 30) {
  covars_use <- covars
  completed_list <- mice::complete(imp_obj, action = "all")
  
  formula_text <- paste0(
    "Surv(stime_years, MACEPLUS) ~ ",
    paste(c(joint_var, covars_use), collapse = " + ")
  )
  f <- as.formula(formula_text)
  
  fit_list <- list()
  n_vec <- c()
  event_vec <- c()
  
  for (i in seq_along(completed_list)) {
    dat_i <- completed_list[[i]] %>%
      create_joint_groups()
    
    needed <- unique(c("stime_years", "MACEPLUS", joint_var, covars_use))
    dat_i <- safe_drop_na(dat_i, needed)
    
    if (nrow(dat_i) < min_n) next
    if (length(unique(dat_i[[joint_var]])) < 2) next
    
    fit_i <- try(coxph(f, data = dat_i), silent = TRUE)
    if (inherits(fit_i, "try-error")) next
    
    fit_list[[length(fit_list) + 1]] <- fit_i
    n_vec <- c(n_vec, fit_i$n)
    event_vec <- c(event_vec, fit_i$nevent)
  }
  
  if (length(fit_list) == 0) return(NULL)
  
  fit_mira <- mice::as.mira(fit_list)
  pooled_fit <- mice::pool(fit_mira)
  
  pooled_summary <- summary(
    pooled_fit,
    conf.int = TRUE,
    conf.level = 0.95,
    exponentiate = TRUE
  )
  
  if (!"term" %in% names(pooled_summary)) return(NULL)
  
  pooled_summary$term <- as.character(pooled_summary$term)
  
  pooled_summary %>%
    dplyr::filter(base::startsWith(term, joint_var)) %>%
    dplyr::mutate(
      analysis_type = "Joint_group",
      dataset = dataset_name,
      joint_var = joint_var,
      joint_label = joint_label,
      reference_group = reference_group,
      model = model_name,
      n = round(mean(n_vec)),
      events = round(mean(event_vec)),
      covariates = ifelse(length(covars_use) == 0, "Unadjusted", paste(covars_use, collapse = ";"))
    ) %>%
    dplyr::mutate(
      comparison_group = gsub(paste0("^", joint_var), "", term)
    ) %>%
    dplyr::transmute(
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
      CI_low = `2.5 %`,
      CI_high = `97.5 %`,
      SE = std.error,
      z = statistic,
      p = p.value
    )
}

continuous_results_mi_list <- list()

for (i in seq_len(nrow(continuous_exposures_surv))) {
  imp_name <- continuous_exposures_surv$imp_obj[i]
  dataset_name <- continuous_exposures_surv$dataset_name[i]
  exposure_i <- continuous_exposures_surv$exposure[i]
  exposure_label_i <- continuous_exposures_surv$exposure_label[i]
  
  imp_i <- get(imp_name)
  
  for (mname in names(model_covs_surv)) {
    res <- run_cox_continuous_mi(
      imp_obj = imp_i,
      exposure = exposure_i,
      exposure_label = exposure_label_i,
      covars = model_covs_surv[[mname]],
      model_name = mname,
      dataset_name = dataset_name,
      min_n = 30
    )
    
    if (!is.null(res) && nrow(res) > 0) {
      continuous_results_mi_list[[length(continuous_results_mi_list) + 1]] <- res
    }
  }
}

continuous_results_mi <- bind_rows(continuous_results_mi_list) %>%
  arrange(exposure, model)

joint_results_mi_list <- list()

for (i in seq_len(nrow(joint_exposures_surv))) {
  imp_name <- joint_exposures_surv$imp_obj[i]
  dataset_name <- joint_exposures_surv$dataset_name[i]
  joint_var_i <- joint_exposures_surv$joint_var[i]
  joint_label_i <- joint_exposures_surv$joint_label[i]
  reference_i <- joint_exposures_surv$reference_group[i]
  
  imp_i <- get(imp_name)
  
  for (mname in names(model_covs_surv)) {
    res <- run_cox_joint_mi(
      imp_obj = imp_i,
      joint_var = joint_var_i,
      joint_label = joint_label_i,
      reference_group = reference_i,
      covars = model_covs_surv[[mname]],
      model_name = mname,
      dataset_name = dataset_name,
      min_n = 30
    )
    
    if (!is.null(res) && nrow(res) > 0) {
      joint_results_mi_list[[length(joint_results_mi_list) + 1]] <- res
    }
  }
}

joint_results_mi <- bind_rows(joint_results_mi_list) %>%
  arrange(joint_var, model, comparison_group)

continuous_results_mi <- continuous_results_mi %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    p_text = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

joint_results_mi <- joint_results_mi %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    p_text = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

############################################################
## 4. Survival MI forest plots
############################################################

forest_surv_cont_mi <- continuous_results_mi %>%
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

forest_surv_joint_mi <- joint_results_mi %>%
  mutate(
    model = factor(model, levels = c("Model1", "Model2", "Model3", "Model4", "Model5")),
    joint_label = factor(joint_label, levels = c("SpO2 + ULK1")),
    comparison_group = case_when(
      joint_label == "SpO2 + ULK1" & comparison_group == "Low SpO2 + High ULK1"  ~ "Group 2: Low SpO2 + High ULK1",
      joint_label == "SpO2 + ULK1" & comparison_group == "High SpO2 + Low ULK1"  ~ "Group 3: High SpO2 + Low ULK1",
      joint_label == "SpO2 + ULK1" & comparison_group == "Low SpO2 + Low ULK1"   ~ "Group 4: Low SpO2 + Low ULK1",
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

plot_surv_cont_row_mi <- function(data, exposure_name, color, tag_letter) {
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
      x = 0.5, y = 0.98,
      gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "Arial")
    )
  )
  
  exposure_panel <- ggplot(data.frame(y = mean(df$y)), aes(y = y)) +
    geom_text(
      aes(x = 0, label = exposure_name),
      hjust = 0.5,
      family = "Arial",
      size = 3.3
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Exposure",
      hjust = 0.5,
      family = "Arial",
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
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.7, 0.7), clip = "off") +
    theme_void()
  
  forest_panel <- ggplot(df, aes(y = y, x = HR, xmin = CI_low, xmax = CI_high)) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.4
    ) +
    geom_point(
      data = df %>% filter(sig),
      aes(x = HR, y = y),
      inherit.aes = FALSE,
      size = 4.2,
      color = alpha(color, 0.18)
    ) +
    geom_errorbarh(
      height = 0.10,
      linewidth = 0.5,
      color = color
    ) +
    geom_point(
      size = 1.6,
      color = color
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = pretty(c(x_min, x_max), n = 5)
    ) +
    scale_y_continuous(
      limits = c(y_bottom, y_top + 0.5),
      breaks = NULL
    ) +
    labs(
      title = exposure_name,
      x = "Hazard ratio (95% CI)",
      y = NULL
    ) +
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
      family = "Arial",
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
      family = "Arial",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.45, 0.45), clip = "off") +
    theme_void()
  
  tag_panel +
    exposure_panel +
    model_panel +
    forest_panel +
    hr_panel +
    p_panel +
    plot_layout(widths = c(0.35, 1.10, 1.05, 4.35, 2.00, 0.85))
}

plot_surv_joint_row_mi <- function(data, joint_name, color_map, tag_letter) {
  df <- data %>%
    filter(joint_label == joint_name) %>%
    arrange(comparison_group, model)
  
  df$comparison_group <- factor(
    df$comparison_group,
    levels = unique(as.character(df$comparison_group))
  )
  
  comp_levels <- levels(df$comparison_group)
  model_step <- 1.05
  block_gap <- 1.15
  
  y_vals <- numeric(nrow(df))
  current_top <- (length(comp_levels) - 1) * (5 * model_step + block_gap)
  
  for (i in seq_along(comp_levels)) {
    idx <- which(df$comparison_group == comp_levels[i])
    y_vals[idx] <- seq(
      from = current_top + 4 * model_step,
      by = -model_step,
      length.out = 5
    )
    current_top <- current_top - (5 * model_step + block_gap)
  }
  
  df$y <- y_vals
  
  comp_df <- df %>%
    group_by(comparison_group) %>%
    summarise(
      y_center = mean(y),
      comp_display = str_wrap(first(as.character(comparison_group)), width = 20),
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
    mutate(
      plot_color = recode(as.character(comparison_group), !!!color_map)
    )
  
  tag_panel <- wrap_elements(
    full = textGrob(
      tag_letter,
      x = 0.5, y = 0.98,
      gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "Arial")
    )
  )
  
  comp_panel <- ggplot(comp_df, aes(y = y_center)) +
    geom_text(
      aes(x = 0, label = comp_display),
      hjust = 0.5,
      family = "Arial",
      size = 3.0
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "Comparison",
      hjust = 0.5,
      family = "Arial",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-1.1, 1.1), clip = "off") +
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
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.75, 0.75), clip = "off") +
    theme_void()
  
  forest_panel <- ggplot(df, aes(y = y, x = HR, xmin = CI_low, xmax = CI_high)) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.4
    ) +
    geom_point(
      data = df %>% filter(sig),
      aes(x = HR, y = y, color = plot_color),
      inherit.aes = FALSE,
      size = 4.2,
      alpha = 0.18,
      show.legend = FALSE
    ) +
    geom_errorbarh(
      aes(color = plot_color),
      height = 0.10,
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    geom_point(
      aes(color = plot_color),
      size = 1.6,
      show.legend = FALSE
    ) +
    scale_color_identity() +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = pretty(c(x_min, x_max), n = 5)
    ) +
    scale_y_continuous(
      limits = c(y_bottom, y_top + 0.5),
      breaks = NULL
    ) +
    labs(
      title = joint_name,
      x = "Hazard ratio (95% CI)",
      y = NULL
    ) +
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
      size = 2.7
    ) +
    annotate(
      "text",
      x = 0, y = y_top,
      label = "HR (95% CI)",
      hjust = 0.5,
      family = "Arial",
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
      family = "Arial",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(limits = c(y_bottom, y_top + 0.5), breaks = NULL) +
    coord_cartesian(xlim = c(-0.45, 0.45), clip = "off") +
    theme_void()
  
  tag_panel +
    comp_panel +
    model_panel +
    forest_panel +
    hr_panel +
    p_panel +
    plot_layout(widths = c(0.35, 1.95, 1.05, 4.60, 2.20, 0.85))
}

row_A_surv_mi <- plot_surv_cont_row_mi(
  forest_surv_cont_mi,
  "SpO2",
  ASAP_colors$SpO2_line,
  "A"
)

row_B_surv_mi <- plot_surv_cont_row_mi(
  forest_surv_cont_mi,
  "ln(ULK1)",
  ASAP_colors$ULK1_line,
  "B"
)

Figure_survival_continuous_forest_mi <- row_A_surv_mi / row_B_surv_mi

spo2_ulk1_colors <- c(
  "Group 2: Low SpO2 + High ULK1" = ASAP_colors$ULK1_line,
  "Group 3: High SpO2 + Low ULK1" = ASAP_colors$SpO2_line,
  "Group 4: Low SpO2 + Low ULK1"  = "#A6A6A6"
)

Figure_survival_joint_forest_mi <- plot_surv_joint_row_mi(
  forest_surv_joint_mi,
  "SpO2 + ULK1",
  spo2_ulk1_colors,
  ""
)

############################################################
## 5. Export selected outputs
############################################################

cross_out_mi <- sanitize_for_xlsx(results_cross_mi)
continuous_out_mi <- sanitize_for_xlsx(continuous_results_mi)
joint_out_mi <- sanitize_for_xlsx(joint_results_mi)
miss_cross_out <- sanitize_for_xlsx(mi_missing_summary_cross)
miss_all_out <- sanitize_for_xlsx(mi_surv_all_missing)
miss_ulk1_out <- sanitize_for_xlsx(mi_surv_ulk1_missing)

write_xlsx(
  list(
    MI_missing_cross = miss_cross_out,
    Cross_sectional_MI_results = cross_out_mi,
    MI_missing_surv_all = miss_all_out,
    MI_missing_surv_ulk1 = miss_ulk1_out,
    Continuous_Cox_MI = continuous_out_mi,
    Joint_Group_Cox_MI = joint_out_mi
  ),
  path = "ASAP_MI_sensitivity_outputs.xlsx"
)

############################################################
## 6. Export figures
############################################################

ggsave(
  filename = "Figure_cross_sectional_forest_MI.png",
  plot = Figure_cross_forest_mi,
  width = 250,
  height = 210,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = "Figure_cross_sectional_forest_MI.tiff",
  plot = Figure_cross_forest_mi,
  width = 250,
  height = 210,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Figure_survival_continuous_forest_MI.png",
  plot = Figure_survival_continuous_forest_mi,
  width = 200,
  height = 160,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = "Figure_survival_continuous_forest_MI.tiff",
  plot = Figure_survival_continuous_forest_mi,
  width = 200,
  height = 160,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Figure_survival_joint_forest_MI.png",
  plot = Figure_survival_joint_forest_mi,
  width = 200,
  height = 160,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = "Figure_survival_joint_forest_MI.tiff",
  plot = Figure_survival_joint_forest_mi,
  width = 200,
  height = 160,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

message("Multiple imputation sensitivity analyses completed successfully.")