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
recipe_all <- readRDS("SpotSepsis Data/Derived/svm_ml_recipe_blueprint.rds")
recipe_bio <- readRDS("SpotSepsis Data/Derived/svm_ml_recipe_bio.rds")  
recipe_clin <- readRDS("SpotSepsis Data/Derived/svm_ml_recipe_clin.rds") 
folds     <- readRDS("SpotSepsis Data/Derived/folds5x5.rds")

# SVM Shared specs ----
# Support Vector Machine with a Radial Basis Function (RBF) kernel,
svm_spec <- svm_rbf(
  mode      = "classification",
  cost      = tune(),      # The cost of misclassifications.
  rbf_sigma = tune()       # The RBF kernel's scale parameter.
) %>% 
  set_engine("kernlab") # Using the 'kernlab' package as the engine.


# tuning grid ranges (log)
cost_range      <- cost(range = c(-5, 3), trans = log10_trans())    
rbf_sigma_range <- rbf_sigma(range = c(-5, 0), trans = log10_trans())  

# random 64‑point grid
set.seed(123)
svm_grid <- grid_random(
  cost_range,
  rbf_sigma_range,
  size = 64
)

print(svm_grid, n = 64)

# parallel compute ----
library(future)
plan(multisession, workers = parallel::detectCores() - 1)
plan()

# workflows for SVM ----
svm_workflow_all <- workflow() %>% 
  add_recipe(recipe_all) %>% 
  add_model(svm_spec)          

svm_workflow_clin <- workflow() %>% 
  add_recipe(recipe_clin) %>% 
  add_model(svm_spec)

svm_workflow_bio <- workflow() %>% 
  add_recipe(recipe_bio) %>% 
  add_model(svm_spec)

# Clin vars ----
set.seed(55378008)

tic()
svm_res <- tune_grid(
  svm_workflow_clin, 
  resamples = folds,
  grid      = svm_grid, 
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
all_preds <- collect_predictions(svm_res)

write.csv(all_preds, "Trained Models/SVM/clin_svm_all_preds.csv")

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

hyperparam_values <- collect_metrics(svm_res) %>%
  distinct(.config, cost, rbf_sigma)


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
  
  # Tier 2 (SVM Simplicity): Of those, find the simplest (lowest cost)
  slice_min(order_by = cost, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 3 (SVM Simplicity): Of those, find the simplest (highest rbf_sigma)
  slice_max(order_by = rbf_sigma, n = 1, with_ties = TRUE, na_rm = TRUE)

if (nrow(best_params_from_manual_results) > 1) {
  warning("A perfect tie still exists after all rules. Selecting the first model from the final tied group.")
  
  final_winner_row <- best_params_from_manual_results %>% 
    slice(1)
} else {
  final_winner_row <- best_params_from_manual_results
}

svm_champion_params <- final_winner_row %>%
  select(cost, rbf_sigma)

print(svm_champion_params)

best_svm_preds <- collect_predictions(svm_res, parameters = svm_champion_params)

# plotting PR curve
best_svm_preds %>%
  pr_curve(truth = SFI_5cat, `.pred_Non-Severe`, `.pred_Probable Non-Severe`, 
           `.pred_Probable Severe`, `.pred_Onset greater than 24 hours`, 
           `.pred_Onset within 24 hours`) %>%
  autoplot() +
  labs(title = "Precision-Recall Curve for the Best (clinical feature) SVM Model")

write.csv(best_svm_preds, "Trained Models/SVM/clin_svm_predictions_best.csv")

# workflow for best performing model
final_svm_wf <- finalize_workflow(svm_workflow_clin, svm_champion_params)

# fitting finalised workflow to entire training set
svm_CLIN_MODEL <- fit(final_svm_wf, data = raw_train)

# Saving full model
saveRDS(svm_CLIN_MODEL, "Trained Models/SVM/svm_CLIN_MODEL.rds")


# [CLIN ONLY] feature importance extraction using bespoke multi-class permutation method ----
svm_fit <- svm_CLIN_MODEL %>%
  pull_workflow_fit()

final_recipe <- svm_CLIN_MODEL %>% 
  extract_recipe()

baked_train_data <- bake(final_recipe, new_data = raw_train)

col_name_map <- list(
  "Non-Severe" = ".pred_Non-Severe",
  "Probable Non-Severe" = ".pred_Probable Non-Severe",
  "Probable Severe" = ".pred_Probable Severe",
  "Onset greater than 24 hours" = ".pred_Onset greater than 24 hours",
  "Onset within 24 hours" = ".pred_Onset within 24 hours"
)

# outcome classes 
outcome_levels <- names(col_name_map)

# empty list to store the importance table for each class
per_class_importance_list <- list()

message("Starting final permutation importance analysis for all 5 classes...")
tic() 

for (outcome_level in outcome_levels) {
  
  message(paste("--- Calculating importance for class:", outcome_level, "---"))
  
  prob_col_name <- col_name_map[[outcome_level]]
  
  pred_wrapper <- function(object, newdata) {
    predict(object, new_data = newdata, type = "prob")[[prob_col_name]]
  }
  
  baseline_score <- pr_auc_vec(
    truth = factor(ifelse(baked_train_data$SFI_5cat == outcome_level, "event", "other")),
    estimate = pred_wrapper(svm_fit, baked_train_data)
  )
  
  importance_results <- list()
  predictor_names <- baked_train_data %>% select(-SFI_5cat) %>% colnames()
  for (feature in predictor_names) {
    set.seed(345) 
    shuffled_data <- baked_train_data %>% mutate("{feature}" := sample(.data[[feature]]))
    
    shuffled_score <- pr_auc_vec(
      truth = factor(ifelse(shuffled_data$SFI_5cat == outcome_level, "event", "other")),
      estimate = pred_wrapper(svm_fit, shuffled_data)
    )
    
    importance_results[[feature]] <- baseline_score - shuffled_score
  }
  
  list_name_clean <- gsub("[^[:alnum:]]", "_", outcome_level) 
  class_importance_df <- importance_results %>%
    enframe(name = "Variable", value = paste0("Importance_", list_name_clean))
  
  per_class_importance_list[[outcome_level]] <- class_importance_df
}

toc() 

# combining results
final_importance_table <- reduce(
  per_class_importance_list,
  full_join,
  by = "Variable"
) %>%
  rowwise() %>%
  mutate(Combined_Importance = mean(c_across(starts_with("Importance_")), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(Combined_Importance))

svm_final_importance <- final_importance_table %>%
  select(Variable, Importance = Combined_Importance) %>%
  arrange(desc(Importance))

write.csv(svm_final_importance, "Trained Models/SVM/clin_svm_feature_importance.csv")

# All vars ----
set.seed(55378008)

tic()
svm_res <- tune_grid(
  svm_workflow_all, 
  resamples = folds,
  grid      = svm_grid, 
  metrics   = metric_set(roc_auc, pr_auc, accuracy, sens, yardstick::spec), 
  control   = control_grid(
    save_pred = TRUE,
    verbose   = TRUE,
    allow_par = TRUE
  )
)
toc() # 5114.324 sec elapsed

# Results Processing
# predictions from the tuning results for every model
all_preds <- collect_predictions(svm_res)

write.csv(all_preds, "Trained Models/SVM/all_feature_svm_all_preds.csv")

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

hyperparam_values <- collect_metrics(svm_res) %>%
  distinct(.config, cost, rbf_sigma)


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
  
  # Tier 2 (SVM Simplicity): Of those, find the simplest (lowest cost)
  slice_min(order_by = cost, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 3 (SVM Simplicity): Of those, find the simplest (highest rbf_sigma)
  slice_max(order_by = rbf_sigma, n = 1, with_ties = TRUE, na_rm = TRUE)

if (nrow(best_params_from_manual_results) > 1) {
  warning("A perfect tie still exists after all rules. Selecting the first model from the final tied group.")
  
  final_winner_row <- best_params_from_manual_results %>% 
    slice(1)
} else {
  final_winner_row <- best_params_from_manual_results
}

svm_champion_params <- final_winner_row %>%
  select(cost, rbf_sigma)

print(svm_champion_params)

best_svm_preds <- collect_predictions(svm_res, parameters = svm_champion_params)

# plotting PR curve
best_svm_preds %>%
  pr_curve(truth = SFI_5cat, `.pred_Non-Severe`, `.pred_Probable Non-Severe`, 
           `.pred_Probable Severe`, `.pred_Onset greater than 24 hours`, 
           `.pred_Onset within 24 hours`) %>%
  autoplot() +
  labs(title = "Precision-Recall Curve for the Best (all feature) SVM Model")

write.csv(best_svm_preds, "Trained Models/SVM/all_feature_svm_predictions_best.csv")

# workflow for best performing model
final_svm_wf <- finalize_workflow(svm_workflow_all, svm_champion_params)

# fitting finalised workflow to entire training set
svm_FULL_MODEL <- fit(final_svm_wf, data = raw_train)

# Saving full model
saveRDS(svm_FULL_MODEL, "Trained Models/SVM/svm_FULL_MODEL.rds")

# [ALL VAR] feature importance extraction using bespoke multi-class permutation method ----
svm_fit <- svm_FULL_MODEL %>%
  pull_workflow_fit()

final_recipe <- svm_FULL_MODEL %>% 
  extract_recipe()

baked_train_data <- bake(final_recipe, new_data = raw_train)

col_name_map <- list(
  "Non-Severe" = ".pred_Non-Severe",
  "Probable Non-Severe" = ".pred_Probable Non-Severe",
  "Probable Severe" = ".pred_Probable Severe",
  "Onset greater than 24 hours" = ".pred_Onset greater than 24 hours",
  "Onset within 24 hours" = ".pred_Onset within 24 hours"
)

# outcome classes 
outcome_levels <- names(col_name_map)

# empty list to store the importance table for each class
per_class_importance_list <- list()

message("Starting final permutation importance analysis for all 5 classes...")
tic() 

for (outcome_level in outcome_levels) {
  
  message(paste("--- Calculating importance for class:", outcome_level, "---"))
  
  prob_col_name <- col_name_map[[outcome_level]]
  
  pred_wrapper <- function(object, newdata) {
    predict(object, new_data = newdata, type = "prob")[[prob_col_name]]
  }
  
  baseline_score <- pr_auc_vec(
    truth = factor(ifelse(baked_train_data$SFI_5cat == outcome_level, "event", "other")),
    estimate = pred_wrapper(svm_fit, baked_train_data)
  )
  
  importance_results <- list()
  predictor_names <- baked_train_data %>% select(-SFI_5cat) %>% colnames()
  for (feature in predictor_names) {
    set.seed(345) 
    shuffled_data <- baked_train_data %>% mutate("{feature}" := sample(.data[[feature]]))
    
    shuffled_score <- pr_auc_vec(
      truth = factor(ifelse(shuffled_data$SFI_5cat == outcome_level, "event", "other")),
      estimate = pred_wrapper(svm_fit, shuffled_data)
    )
    
    importance_results[[feature]] <- baseline_score - shuffled_score
  }
  
  list_name_clean <- gsub("[^[:alnum:]]", "_", outcome_level) 
  class_importance_df <- importance_results %>%
    enframe(name = "Variable", value = paste0("Importance_", list_name_clean))
  
  per_class_importance_list[[outcome_level]] <- class_importance_df
}

toc() 

# combining results
final_importance_table <- reduce(
  per_class_importance_list,
  full_join,
  by = "Variable"
) %>%
  rowwise() %>%
  mutate(Combined_Importance = mean(c_across(starts_with("Importance_")), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(Combined_Importance))

svm_final_importance <- final_importance_table %>%
  select(Variable, Importance = Combined_Importance) %>%
  arrange(desc(Importance))

write.csv(svm_final_importance, "Trained Models/SVM/all_feature_svm_feature_importance.csv")
# Bio vars ----
set.seed(55378008)

tic()
svm_res <- tune_grid(
  svm_workflow_bio, 
  resamples = folds,
  grid      = svm_grid, 
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
all_preds <- collect_predictions(svm_res)

write.csv(all_preds, "Trained Models/SVM/bio_svm_all_preds.csv")

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

hyperparam_values <- collect_metrics(svm_res) %>%
  distinct(.config, cost, rbf_sigma)


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
  
  # Tier 2 (SVM Simplicity): Of those, find the simplest (lowest cost)
  slice_min(order_by = cost, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
  
  # Tier 3 (SVM Simplicity): Of those, find the simplest (highest rbf_sigma)
  slice_max(order_by = rbf_sigma, n = 1, with_ties = TRUE, na_rm = TRUE)

if (nrow(best_params_from_manual_results) > 1) {
  warning("A perfect tie still exists after all rules. Selecting the first model from the final tied group.")
  
  final_winner_row <- best_params_from_manual_results %>% 
    slice(1)
} else {
  final_winner_row <- best_params_from_manual_results
}

svm_champion_params <- final_winner_row %>%
  select(cost, rbf_sigma)

print(svm_champion_params)

best_svm_preds <- collect_predictions(svm_res, parameters = svm_champion_params)

# plotting PR curve
best_svm_preds %>%
  pr_curve(truth = SFI_5cat, `.pred_Non-Severe`, `.pred_Probable Non-Severe`, 
           `.pred_Probable Severe`, `.pred_Onset greater than 24 hours`, 
           `.pred_Onset within 24 hours`) %>%
  autoplot() +
  labs(title = "Precision-Recall Curve for the Best (lab feature) SVM Model")

write.csv(best_svm_preds, "Trained Models/SVM/bio_svm_predictions_best.csv")

# workflow for best performing model
final_svm_wf <- finalize_workflow(svm_workflow_bio, svm_champion_params)

# fitting finalised workflow to entire training set
svm_BIO_MODEL <- fit(final_svm_wf, data = raw_train)

# Saving full model
saveRDS(svm_BIO_MODEL, "Trained Models/SVM/svm_BIO_MODEL.rds")



# [BIO ONLY] feature importance extraction using bespoke multi-class permutation method ----
svm_fit <- svm_BIO_MODEL %>%
  pull_workflow_fit()

final_recipe <- svm_BIO_MODEL %>% 
  extract_recipe()

baked_train_data <- bake(final_recipe, new_data = raw_train)

col_name_map <- list(
  "Non-Severe" = ".pred_Non-Severe",
  "Probable Non-Severe" = ".pred_Probable Non-Severe",
  "Probable Severe" = ".pred_Probable Severe",
  "Onset greater than 24 hours" = ".pred_Onset greater than 24 hours",
  "Onset within 24 hours" = ".pred_Onset within 24 hours"
)

# outcome classes 
outcome_levels <- names(col_name_map)

# empty list to store the importance table for each class
per_class_importance_list <- list()

message("Starting final permutation importance analysis for all 5 classes...")
tic() 

for (outcome_level in outcome_levels) {
  
  message(paste("--- Calculating importance for class:", outcome_level, "---"))
  
  prob_col_name <- col_name_map[[outcome_level]]
  
  pred_wrapper <- function(object, newdata) {
    predict(object, new_data = newdata, type = "prob")[[prob_col_name]]
  }
  
  baseline_score <- pr_auc_vec(
    truth = factor(ifelse(baked_train_data$SFI_5cat == outcome_level, "event", "other")),
    estimate = pred_wrapper(svm_fit, baked_train_data)
  )
  
  importance_results <- list()
  predictor_names <- baked_train_data %>% select(-SFI_5cat) %>% colnames()
  for (feature in predictor_names) {
    set.seed(345) 
    shuffled_data <- baked_train_data %>% mutate("{feature}" := sample(.data[[feature]]))
    
    shuffled_score <- pr_auc_vec(
      truth = factor(ifelse(shuffled_data$SFI_5cat == outcome_level, "event", "other")),
      estimate = pred_wrapper(svm_fit, shuffled_data)
    )
    
    importance_results[[feature]] <- baseline_score - shuffled_score
  }
  
  list_name_clean <- gsub("[^[:alnum:]]", "_", outcome_level) 
  class_importance_df <- importance_results %>%
    enframe(name = "Variable", value = paste0("Importance_", list_name_clean))
  
  per_class_importance_list[[outcome_level]] <- class_importance_df
}

toc() 

# combining results
final_importance_table <- reduce(
  per_class_importance_list,
  full_join,
  by = "Variable"
) %>%
  rowwise() %>%
  mutate(Combined_Importance = mean(c_across(starts_with("Importance_")), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(Combined_Importance))

svm_final_importance <- final_importance_table %>%
  select(Variable, Importance = Combined_Importance) %>%
  arrange(desc(Importance))

write.csv(svm_final_importance, "Trained Models/SVM/bio_svm_feature_importance.csv")

