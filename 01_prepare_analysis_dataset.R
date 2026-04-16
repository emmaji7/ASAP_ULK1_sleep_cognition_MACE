############################################################
## 01_prepare_analysis_dataset.R
##
## Purpose:
## Prepare analysis-ready datasets for the ASAP manuscript.
##
## Note:
## The underlying participant-level data are restricted and are
## not distributed with this repository. Users must first obtain
## the required approvals and load the source dataset into R.
##
## Expected input:
## A data frame named `asap_raw` containing the required variables.
##
## Output objects created in memory:
## - analysis_data
## - surv_data_all
## - surv_data_ulk1
## - cross_data_cog_master
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260304)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(haven)
})

############################################################
## Helper functions
############################################################

to_num <- function(x) {
  if (inherits(x, "labelled")) x <- haven::zap_labels(x)
  suppressWarnings(as.numeric(as.character(x)))
}

############################################################
## Input check
############################################################

if (!exists("asap_raw")) {
  stop(
    paste(
      "Input object `asap_raw` was not found.",
      "Please load the restricted source dataset into R before running this script."
    )
  )
}

############################################################
## Participant-level master dataset
############################################################

analysis_master <- asap_raw %>%
  mutate(
    ID_std = str_pad(as.character(PatientID), width = 6, pad = "0")
  ) %>%
  arrange(ID_std) %>%
  distinct(ID_std, .keep_all = TRUE)

############################################################
## Define ULK1 availability
############################################################

analysis_master <- analysis_master %>%
  mutate(
    ULK1_nonmissing = !is.na(a1_lnSerumULK1_1)
  )

############################################################
## Construct analysis variables
############################################################

END_DATE_STATA <- 22400

analysis_data <- analysis_master %>%
  mutate(
    ########################################################
    ## Examination date
    ########################################################
    sXusdag = as.integer(sXusdag),
    sXusmnd = as.integer(sXusmnd),
    sXusar  = as.integer(sXusar),
    
    sXusar_year = case_when(
      sXusar %in% 6:9   ~ 2000L + sXusar,
      sXusar %in% 10:99 ~ 1900L + sXusar,
      TRUE ~ NA_integer_
    ),
    
    exam_date = as.Date(
      sprintf("%04d-%02d-%02d", sXusar_year, sXusmnd, sXusdag)
    ),
    
    exam_date_stata = as.integer(exam_date - as.Date("1960-01-01")),
    
    ########################################################
    ## Event indicators
    ########################################################
    I21 = if_else(I21 == 1, 1L, 0L, missing = 0L),
    I50 = if_else(I50 == 1, 1L, 0L, missing = 0L),
    I64 = if_else(I64 == 1, 1L, 0L, missing = 0L),
    Dod = if_else(Dod == 1, 1L, 0L, missing = 0L),
    
    MACEPLUS = if_else(I21 == 1 | I50 == 1 | I64 == 1 | Dod == 1, 1L, 0L),
    
    t_I21   = if_else(I21 == 1, as.numeric(time_to_I21), NA_real_),
    t_I50   = if_else(I50 == 1, as.numeric(time_to_I50), NA_real_),
    t_I64   = if_else(I64 == 1, as.numeric(time_to_I64), NA_real_),
    t_death = if_else(Dod == 1, as.numeric(time_to_death), NA_real_),
    
    t_event_raw = pmin(t_I21, t_I50, t_I64, t_death, na.rm = TRUE),
    t_event = if_else(is.infinite(t_event_raw), NA_real_, t_event_raw),
    
    t_censor = END_DATE_STATA - exam_date_stata,
    stime_days = if_else(MACEPLUS == 1, t_event, t_censor),
    stime_years = stime_days / 365.24,
    
    ########################################################
    ## Core variables
    ########################################################
    age = to_num(a1age),
    BMI = to_num(a1BMI),
    sex = to_num(sex),
    education = to_num(highereducation),
    
    ########################################################
    ## Sleep / exposure
    ########################################################
    AHI_raw = to_num(a1AHI),
    AHI_bin5 = case_when(
      AHI_raw < 5  ~ 0L,
      AHI_raw >= 5 ~ 1L,
      TRUE ~ NA_integer_
    ),
    
    SpO2 = to_num(a1Average_SpO2),
    
    ########################################################
    ## Biomarker
    ########################################################
    ULK1 = to_num(a1_SerumULK1_1),
    lnULK1 = to_num(a1_lnSerumULK1_1),
    
    ########################################################
    ## Smoking (ONLY keep smoke_exposure)
    ########################################################
    smoking_status = to_num(sXroykestatus),
    
    smoke_exposure = case_when(
      smoking_status %in% c(1, 2) ~ 1L,
      smoking_status == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    ########################################################
    ## Clinical covariates
    ########################################################
    hypertension = to_num(a1BloodP),
    diabetes = to_num(a1diabet),
    CVD = to_num(a1CVD),
    
    TST = to_num(a1TSTh),
    sleep_efficiency = to_num(a1sleep_efficiency),
    cpap_use = to_num(CPAP),
    
    LDL = to_num(LDLkoles),
    high_LDL = case_when(
      LDL >= 3 ~ 1L,
      LDL < 3  ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    obesity = case_when(
      BMI >= 30 ~ 1L,
      BMI < 30  ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    ########################################################
    ## Cognitive outcomes
    ########################################################
    RAVLT_tot = to_num(a1RAVL_sumA),
    RAVLT_learn = to_num(a1RAVL_learn),
    Stroop_t = to_num(a1S3tid),
    Stroop_r = to_num(a1Sratio_10)
  )

############################################################
## Define analysis populations
############################################################

surv_data_all <- analysis_data

surv_data_ulk1 <- analysis_data %>%
  filter(ULK1_nonmissing)

cross_data_cog_master <- analysis_data %>%
  filter(
    !is.na(RAVLT_tot) |
      !is.na(RAVLT_learn) |
      !is.na(Stroop_t) |
      !is.na(Stroop_r)
  )

############################################################
## Final check
############################################################

required_vars <- c(
  "smoke_exposure", "SpO2", "lnULK1",
  "MACEPLUS", "stime_years"
)

missing_vars <- setdiff(required_vars, names(analysis_data))

if (length(missing_vars) > 0) {
  stop(paste("Missing required variables:", paste(missing_vars, collapse = ", ")))
}

message("Analysis-ready datasets created successfully.")