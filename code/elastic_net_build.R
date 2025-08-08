# Set up ----
library(dplyr)
library(ggplot2)
library(rpart)
library(caret)
library(rpart.plot)
library(vip)
library(pdp)
library(tidymodels)
library(tictoc)
library(parsnip)
library(dials)
library(xgboost)

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

# Loading in data/recipe/folds ----
raw_train <- readRDS("SpotSepsis Data/Derived/RAW_train_df.rds")
recipe_all <- readRDS("SpotSepsis Data/Derived/elastic_net_ml_recipe_blueprint.rds")
recipe_bio <- readRDS("SpotSepsis Data/Derived/elastic_net_ml_recipe_bio.rds")  
recipe_clin <- readRDS("SpotSepsis Data/Derived/elastic_net_ml_recipe_clin.rds")  
folds     <- readRDS("SpotSepsis Data/Derived/folds5x5.rds")

# elastic net shared specs ----
penalised_reg_spec <- multinom_reg(
  mode    = "classification",
  penalty = tune(),     
  mixture = tune() # proportion of lasso and ridge penalty
) %>% 
  set_engine("glmnet")


# --- Penalised Regression tuning grid ----
penalised_reg_grid <- grid_regular(
  penalty(), 
  mixture(), 
  levels = 8 # 64 
)

print(penalised_reg_grid)

# parallel compute ----
library(future)
plan(multisession, workers = parallel::detectCores() - 1)
plan()

# workflows for elastic net ----
en_workflow_all <- workflow() %>% 
  add_recipe(recipe_all) %>% 
  add_model(penalised_reg_spec)

en_workflow_clin <- workflow() %>% 
  add_recipe(recipe_clin) %>% 
  add_model(penalised_reg_spec)

en_workflow_bio <- workflow() %>% 
  add_recipe(recipe_bio) %>% 
  add_model(penalised_reg_spec)

# All vars ----
set.seed(55378008)

tic()
en_res <- tune_grid(
  en_workflow_all,
  resamples = folds,
  grid      = penalised_reg_grid, 
  metrics   = metric_set(roc_auc, pr_auc, accuracy, sens, yardstick::spec), 
  control   = control_grid(
    save_pred = TRUE,
    verbose   = TRUE,
    allow_par = TRUE
  )
)
toc()

# Results Processing
# predictions from the tuning results for every model
all_preds <- collect_predictions(en_res)

write.csv(all_preds, "Trained Models/Penalised Regression/all_feature_elastic_net_predictions_all.csv")

# Per-Class PR AUC for Each Model and Each Fold ---
per_fold_per_class_metrics <- all_preds %>%
  group_by(.config, id) %>% 
  summarise(
    pr_auc_Non_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Non-Severe", "event", "other")),
      estimate = .data[[".pred_Non-Severe"]]
    ),
    pr_auc_Probable_Non_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Probable Non-Severe", "event", "other")),
      estimate = .data[[".pred_Probable Non-Severe"]]
    ),
    pr_auc_Probable_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Probable Severe", "event", "other")),
      estimate = .data[[".pred_Probable Severe"]]
    ),
    pr_auc_Onset_gt24h = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Onset greater than 24 hours", "event", "other")),
      estimate = .data[[".pred_Onset greater than 24 hours"]]
    ),
    pr_auc_Onset_lt24h = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Onset within 24 hours", "event", "other")),
      estimate = .data[[".pred_Onset within 24 hours"]]
    ),
    .groups = "drop"
  )


# per-fold scores to get the final result for each model
final_manual_summary <- per_fold_per_class_metrics %>%
  group_by(.config) %>%
  summarise(
    across(
      starts_with("pr_auc_"), 
      ~mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  )

hyperparam_values <- collect_metrics(en_res) %>%
  distinct(.config, penalty, mixture)

final_results_table <- final_manual_summary %>%
  rowwise() %>%
  mutate(
    manual_macro_avg = mean(c_across(ends_with("_mean")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(hyperparam_values, by = ".config") %>%
  arrange(desc(manual_macro_avg))

best_params_from_manual_results <- final_results_table %>%
  
  # Tier 1: Find the model(s) with the highest manual macro-average PR AUC
  slice_max(order_by = manual_macro_avg, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 2 (Elastic Net Simplicity): Of those, find the simplest (highest penalty)
  slice_max(order_by = penalty, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 3 (Elastic Net Simplicity): Of those, find the simplest (highest mixture)
  slice_max(order_by = mixture, n = 1, with_ties = TRUE, na_rm = TRUE)

if (nrow(best_params_from_manual_results) > 1) {
  warning("A perfect tie still exists after all rules. Selecting the first model from the final tied group.")
  
  final_winner_row <- best_params_from_manual_results %>% 
    slice(1)
} else {
  final_winner_row <- best_params_from_manual_results
}

en_champion_params <- final_winner_row %>%
  select(penalty, mixture)

print(en_champion_params)

best_en_preds <- collect_predictions(en_res, parameters = en_champion_params)

# plotting PR curve
best_en_preds %>%
  pr_curve(truth = SFI_5cat, `.pred_Non-Severe`, `.pred_Probable Non-Severe`, 
           `.pred_Probable Severe`, `.pred_Onset greater than 24 hours`, 
           `.pred_Onset within 24 hours`) %>%
  autoplot() +
  labs(title = "Precision-Recall Curve for the Best (all feature) Penalised Regression Model")

write.csv(best_en_preds, "Trained Models/Penalised Regression/all_feature_elastic_net_predictions_best.csv")

# workflow for best performing model
final_en_wf <- finalize_workflow(en_workflow_all, en_champion_params)

# fitting finalised workflow to entire training set
EN_FULL_MODEL <- fit(final_en_wf, data = raw_train)

# feature importance scores
## uses the 'permutation' method specified in the model engine
en_importance_scores <- vi(EN_FULL_MODEL)

print("Top Feature Importance Scores:");print(en_importance_scores)

vip(EN_FULL_MODEL, num_features = 20) +
  labs(title = "Top 20 Most Important Features in the Elastic Net Model")

# Save the importance scores to a new file
write.csv(en_importance_scores, "Trained Models/Penalised Regression/all_feature_elastic_net_feature_importance.csv")

# Saving full model
saveRDS(EN_FULL_MODEL, "Trained Models/Penalised Regression/EN_FULL_MODEL.rds")


# Clin vars only ----
set.seed(55378008)

tic()
en_res <- tune_grid(
  en_workflow_clin,
  resamples = folds,
  grid      = penalised_reg_grid, 
  metrics   = metric_set(roc_auc, pr_auc, accuracy, sens, yardstick::spec), 
  control   = control_grid(
    save_pred = TRUE,
    verbose   = TRUE,
    allow_par = TRUE
  )
)
toc()

# Results Processing
# predictions from the tuning results for every model
all_preds <- collect_predictions(en_res)

write.csv(all_preds, "Trained Models/Penalised Regression/clin_elastic_net_predictions_all.csv")

# Per-Class PR AUC for Each Model and Each Fold ---
per_fold_per_class_metrics <- all_preds %>%
  group_by(.config, id) %>% 
  summarise(
    pr_auc_Non_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Non-Severe", "event", "other")),
      estimate = .data[[".pred_Non-Severe"]]
    ),
    pr_auc_Probable_Non_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Probable Non-Severe", "event", "other")),
      estimate = .data[[".pred_Probable Non-Severe"]]
    ),
    pr_auc_Probable_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Probable Severe", "event", "other")),
      estimate = .data[[".pred_Probable Severe"]]
    ),
    pr_auc_Onset_gt24h = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Onset greater than 24 hours", "event", "other")),
      estimate = .data[[".pred_Onset greater than 24 hours"]]
    ),
    pr_auc_Onset_lt24h = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Onset within 24 hours", "event", "other")),
      estimate = .data[[".pred_Onset within 24 hours"]]
    ),
    .groups = "drop"
  )


# per-fold scores to get the final result for each model
final_manual_summary <- per_fold_per_class_metrics %>%
  group_by(.config) %>%
  summarise(
    across(
      starts_with("pr_auc_"), 
      ~mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  )

hyperparam_values <- collect_metrics(en_res) %>%
  distinct(.config, penalty, mixture)

final_results_table <- final_manual_summary %>%
  rowwise() %>%
  mutate(
    manual_macro_avg = mean(c_across(ends_with("_mean")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(hyperparam_values, by = ".config") %>%
  arrange(desc(manual_macro_avg))

best_params_from_manual_results <- final_results_table %>%
  
  # Tier 1: Find the model(s) with the highest manual macro-average PR AUC
  slice_max(order_by = manual_macro_avg, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 2 (Elastic Net Simplicity): Of those, find the simplest (highest penalty)
  slice_max(order_by = penalty, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 3 (Elastic Net Simplicity): Of those, find the simplest (highest mixture)
  slice_max(order_by = mixture, n = 1, with_ties = TRUE, na_rm = TRUE)

if (nrow(best_params_from_manual_results) > 1) {
  warning("A perfect tie still exists after all rules. Selecting the first model from the final tied group.")
  
  final_winner_row <- best_params_from_manual_results %>% 
    slice(1)
} else {
  final_winner_row <- best_params_from_manual_results
}

en_champion_params <- final_winner_row %>%
  select(penalty, mixture)

print(en_champion_params)

best_en_preds <- collect_predictions(en_res, parameters = en_champion_params)

# plotting PR curve
best_en_preds %>%
  pr_curve(truth = SFI_5cat, `.pred_Non-Severe`, `.pred_Probable Non-Severe`, 
           `.pred_Probable Severe`, `.pred_Onset greater than 24 hours`, 
           `.pred_Onset within 24 hours`) %>%
  autoplot() +
  labs(title = "Precision-Recall Curve for the Best (all feature) Penalised Regression Model")

write.csv(best_en_preds, "Trained Models/Penalised Regression/clin_elastic_net_predictions_best.csv")

# workflow for best performing model
final_en_wf <- finalize_workflow(en_workflow_clin, en_champion_params)

# fitting finalised workflow to entire training set
EN_CLIN_MODEL <- fit(final_en_wf, data = raw_train)

# feature importance scores
## uses the 'permutation' method specified in the model engine
en_importance_scores <- vi(EN_CLIN_MODEL)

print("Top Feature Importance Scores:");print(en_importance_scores)

vip(EN_CLIN_MODEL, num_features = 20) +
  labs(title = "Top 20 Most Important Features in the Elastic Net (Clinical Feature) Model")

# Save the importance scores to a new file
write.csv(en_importance_scores, "Trained Models/Penalised Regression/clin_elastic_net_feature_importance.csv")

# Saving full model
saveRDS(EN_CLIN_MODEL, "Trained Models/Penalised Regression/EN_CLIN_MODEL.rds")


# Lab vars only ----
set.seed(55378008)

tic()
en_res <- tune_grid(
  en_workflow_bio,
  resamples = folds,
  grid      = penalised_reg_grid, 
  metrics   = metric_set(roc_auc, pr_auc, accuracy, sens, yardstick::spec), 
  control   = control_grid(
    save_pred = TRUE,
    verbose   = TRUE,
    allow_par = TRUE
  )
)
toc()

# Results Processing
# predictions from the tuning results for every model
all_preds <- collect_predictions(en_res)

write.csv(all_preds, "Trained Models/Penalised Regression/bio_elastic_net_predictions_all.csv")

# Per-Class PR AUC for Each Model and Each Fold ---
per_fold_per_class_metrics <- all_preds %>%
  group_by(.config, id) %>% 
  summarise(
    pr_auc_Non_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Non-Severe", "event", "other")),
      estimate = .data[[".pred_Non-Severe"]]
    ),
    pr_auc_Probable_Non_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Probable Non-Severe", "event", "other")),
      estimate = .data[[".pred_Probable Non-Severe"]]
    ),
    pr_auc_Probable_Severe = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Probable Severe", "event", "other")),
      estimate = .data[[".pred_Probable Severe"]]
    ),
    pr_auc_Onset_gt24h = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Onset greater than 24 hours", "event", "other")),
      estimate = .data[[".pred_Onset greater than 24 hours"]]
    ),
    pr_auc_Onset_lt24h = pr_auc_vec(
      truth = factor(ifelse(SFI_5cat == "Onset within 24 hours", "event", "other")),
      estimate = .data[[".pred_Onset within 24 hours"]]
    ),
    .groups = "drop"
  )


# per-fold scores to get the final result for each model
final_manual_summary <- per_fold_per_class_metrics %>%
  group_by(.config) %>%
  summarise(
    across(
      starts_with("pr_auc_"), 
      ~mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  )

hyperparam_values <- collect_metrics(en_res) %>%
  distinct(.config, penalty, mixture)

final_results_table <- final_manual_summary %>%
  rowwise() %>%
  mutate(
    manual_macro_avg = mean(c_across(ends_with("_mean")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(hyperparam_values, by = ".config") %>%
  arrange(desc(manual_macro_avg))

best_params_from_manual_results <- final_results_table %>%
  
  # Tier 1: Find the model(s) with the highest manual macro-average PR AUC
  slice_max(order_by = manual_macro_avg, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 2 (Elastic Net Simplicity): Of those, find the simplest (highest penalty)
  slice_max(order_by = penalty, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 3 (Elastic Net Simplicity): Of those, find the simplest (highest mixture)
  slice_max(order_by = mixture, n = 1, with_ties = TRUE, na_rm = TRUE)

if (nrow(best_params_from_manual_results) > 1) {
  warning("A perfect tie still exists after all rules. Selecting the first model from the final tied group.")
  
  final_winner_row <- best_params_from_manual_results %>% 
    slice(1)
} else {
  final_winner_row <- best_params_from_manual_results
}

en_champion_params <- final_winner_row %>%
  select(penalty, mixture)

print(en_champion_params)

best_en_preds <- collect_predictions(en_res, parameters = en_champion_params)

# plotting PR curve
best_en_preds %>%
  pr_curve(truth = SFI_5cat, `.pred_Non-Severe`, `.pred_Probable Non-Severe`, 
           `.pred_Probable Severe`, `.pred_Onset greater than 24 hours`, 
           `.pred_Onset within 24 hours`) %>%
  autoplot() +
  labs(title = "Precision-Recall Curve for the Best (bio feature) Penalised Regression Model")

write.csv(best_en_preds, "Trained Models/Penalised Regression/bio_elastic_net_predictions_best.csv")

# workflow for best performing model
final_en_wf <- finalize_workflow(en_workflow_bio, en_champion_params)

# fitting finalised workflow to entire training set
EN_BIO_MODEL <- fit(final_en_wf, data = raw_train)

# feature importance scores
## uses the 'permutation' method specified in the model engine
en_importance_scores <- vi(EN_BIO_MODEL)

print("Top Feature Importance Scores:");print(en_importance_scores)

vip(EN_BIO_MODEL, num_features = 20) +
  labs(title = "Top 20 Most Important Features in the Elastic Net (Lab Feature) Model")

# Save the importance scores to a new file
write.csv(en_importance_scores, "Trained Models/Penalised Regression/bio_elastic_net_feature_importance.csv")

# Saving full model
saveRDS(EN_BIO_MODEL, "Trained Models/Penalised Regression/EN_BIO_MODEL.rds")

