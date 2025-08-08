# Set up ----
library(readr)
library(dplyr)
library(tidyr)
library(purrr)   
library(janitor)
library(tidyverse)
library(skimr)
library(naniar)
library(caret)
library(ordinal)
library(tidymodels)
library(tictoc)


get_script_path <- function() {
  # RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    sp <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(sp) && sp != "") {
      return(normalizePath(dirname(sp)))
    }
  }
  # Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args)
  if (length(file_arg) > 0) {
    sp <- gsub("--file=", "", args[file_arg])
    return(normalizePath(dirname(sp)))
  }
  
  stop("Could not determine the script path. Please consider setting it manually!")
}

script_path <- get_script_path()

# my function library ----
setwd(script_path)
source(file = file.path("..", "Old R Stuff", "Function_Library.R"))

# Loading in data ----
raw_train <- readRDS("SpotSepsis Data/Derived/RAW_train_df.rds")
folds     <- readRDS("SpotSepsis Data/Derived/folds5x5.rds")

# Parallel compute ----
library(future)
plan(multisession, workers = parallel::detectCores() - 1)

# Recipe ----
baseline_recipe <- 
  recipe(SFI_5cat ~ CRP + PROC + age.months + infection + weight, data = raw_train) %>%
  
  update_role(weight, new_role = "case_weight") %>%
  
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(all_numeric_predictors())

baseline_spec <- multinom_reg(penalty = 0) %>%
  set_engine("glmnet")

baseline_wf <- workflow() %>%
  add_recipe(baseline_recipe) %>%
  add_model(baseline_spec)

# Cross validation ----
set.seed(55378008)
tic()
baseline_cv_res <- fit_resamples(
  baseline_wf,
  resamples = folds,
  metrics   = metric_set(roc_auc, pr_auc, accuracy, sens, spec),
  control   = control_grid(
    save_pred = TRUE,
    verbose   = TRUE,
    allow_par = TRUE
  )
)
toc()

# Results ----
all_baseline_preds <- collect_predictions(baseline_cv_res)

per_fold_per_class_metrics <- all_baseline_preds %>%
  group_by(id) %>%
  summarise(
    pr_auc_Non_Severe = pr_auc_vec(truth = factor(ifelse(SFI_5cat == "Non-Severe", "event", "other")), estimate = .data[[".pred_Non-Severe"]]),
    pr_auc_Probable_Non_Severe = pr_auc_vec(truth = factor(ifelse(SFI_5cat == "Probable Non-Severe", "event", "other")), estimate = .data[[".pred_Probable Non-Severe"]]),
    pr_auc_Probable_Severe = pr_auc_vec(truth = factor(ifelse(SFI_5cat == "Probable Severe", "event", "other")), estimate = .data[[".pred_Probable Severe"]]),
    pr_auc_Onset_gt24h = pr_auc_vec(truth = factor(ifelse(SFI_5cat == "Onset greater than 24 hours", "event", "other")), estimate = .data[[".pred_Onset greater than 24 hours"]]),
    pr_auc_Onset_lt24h = pr_auc_vec(truth = factor(ifelse(SFI_5cat == "Onset within 24 hours", "event", "other")), estimate = .data[[".pred_Onset within 24 hours"]]),
    .groups = "drop"
  )

manual_pr_auc_summary <- per_fold_per_class_metrics %>%
  summarise(
    across(starts_with("pr_auc_"), ~mean(.x, na.rm = TRUE), .names = "{.col}")
  ) %>%
  rowwise() %>%
  mutate(
    manual_macro_avg_pr_auc = mean(c_across(starts_with("pr_auc_")), na.rm = TRUE)
  ) %>%
  ungroup()



other_metrics <- collect_metrics(baseline_cv_res) %>%
  select(.metric, mean) %>%
  pivot_wider(names_from = .metric, values_from = mean) %>%
  rename(
    roc_auc_macro = roc_auc,
    sens_macro = sens,
    spec_macro = spec
  )


final_baseline_results <- bind_cols(manual_pr_auc_summary, other_metrics) %>%
  mutate(
    model_name = "Logistic Regression",
    feature_set = "Tier 1 Literature"
  ) %>%
  # Reorder columns for a clean final output
  select(
    model_name,
    feature_set,
    manual_macro_avg_pr_auc,
    roc_auc_macro,
    sens_macro,
    spec_macro,
    # Select the per-class PR AUCs, using any_of() for safety
    any_of(c("pr_auc_Non_Severe", "pr_auc_Probable_Non_Severe", "pr_auc_Probable_Severe",
             "pr_auc_Onset_gt24h", "pr_auc_Onset_lt24h", "accuracy"))
  )

# View the final, correct summary table for your baseline model
print(final_baseline_results)

# Save the final results to a CSV file
write_csv(final_baseline_results, "Trained Models/DAG_informed/baseline_model_final_metrics.csv")
