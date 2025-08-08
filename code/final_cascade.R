# Set up ----
library(tidyverse)
library(tidymodels)
library(probably)
library(tictoc)

full_results_table <- read_csv("Trained Models/all_ml_model_results.csv", show_col_types = FALSE)

raw_train <- readRDS("SpotSepsis Data/Derived/RAW_train_df.rds")
raw_test <- readRDS("SpotSepsis Data/Derived/RAW_test_df.rds")

# Bespoke function ----
# class with max probability per row
max_prob_class <- function(df) {
  prob_cols <- df %>% dplyr::select(starts_with(".pred_"))
  levels   <- gsub("\\.pred_", "", colnames(prob_cols))
  factor(levels[max.col(prob_cols)], levels = levels)
}

# Predictor sets ----
clin_predictors <- c(
  "site", "age.months", "sex", "bgcombyn", "adm.recent", "wfaz", "waste", "stunt",
  "cidysymp", "prior.care", "travel.time.bin", "urti", "lrti", "diarrhoeal",
  "neuro", "auf", "syndrome.resp", "syndrome.nonresp", "pneumo", "sev.pneumo",
  "ensapro", "vomit.all", "seiz", "pfacleth", "not.alert", "danger.sign",
  "hr.all", "rr.all", "oxy.ra", "envhtemp", "crt.long", "LqSOFA",
  "parenteral_screen", "SIRS"
)

all_predictors_list <- c(
  clin_predictors, "weight", "ANG1", "ANG2", "CHI3L", "CRP", "CXCl10", "IL1ra",
  "IL6", "IL8", "IL10", "PROC", "TNFR1", "STREM1", "VEGFR1", "supar", "lblac",
  "lbglu", "enescbchb1"
)

# Stage one model hyperparameters (rand forest) ----
s1_champion_spec <- rand_forest(
  mode  = "classification",
  mtry  = 17,
  min_n = 2,
  trees = 1000 
) %>% 
  set_engine("ranger", importance = "permutation")

# Specialist meta-learner model hyperparamters ----
s2_specialist_specs <- list(
  "Probable_Severe" = boost_tree(mode = "classification", trees = 1000,
                                 tree_depth = 11, min_n = 11, learn_rate = 0.006878251, mtry = 7, sample_size = 0.6941596) %>%
    set_engine("xgboost"),
  "Onset_greater_than_24_hours" = boost_tree(mode = "classification", trees = 1000,
                                             tree_depth = 5, min_n = 13, learn_rate = 0.017956412,
                                             mtry = 42, sample_size = 0.9170272) %>%
    set_engine("xgboost"), 
  "Non_Severe" = boost_tree(mode = "classification", trees = 1000,
                            tree_depth = 14, min_n = 7, learn_rate = 0.007688208, mtry = 37, sample_size = 0.6773915) %>% set_engine("xgboost"), 
  "Onset_within_24_hours" = rand_forest(mode = "classification", trees = 1000, mtry = 17, min_n = 2) %>% set_engine("ranger")  
)


# Recipes for stage one and stage two model ----
# Stage 1 (clinical predictors only)
s1_recipe <- recipe(SFI_5cat ~ ., data = raw_train) %>%
  step_rm(weight_bin) %>%
  update_role(weight, new_role = "case_weight") %>%
  step_rm(-all_of(clin_predictors), -has_role("case_weight"), -all_outcomes()) %>%
  step_indicate_na(all_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_zv(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors())

# Stage 2 (all predictors)
recipe_all <- recipe(SFI_5cat ~ ., data = raw_train) %>%
  step_rm(weight_bin) %>%
  
  update_role(weight, new_role = "case_weight") %>%
  step_rm(-all_of(all_predictors_list), -has_role("case_weight"), -all_outcomes()) %>%
  step_indicate_na(all_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_log(ANG1, ANG2, CHI3L, CRP, CXCl10, IL1ra, IL6, IL8,
           IL10, PROC, TNFR1, STREM1, VEGFR1, supar, offset = 1)


# Calibrating stage one model and tuning handoff thresholds ----
message("Calibrating Stage 1 and tuning handoff thresholds...")

# subset of the training data for calibration;
set.seed(55378008)

cal_split <- initial_split(raw_train, prop = 0.8, strata = SFI_5cat)
train_sub <- training(cal_split)
cal_sub   <- testing(cal_split)

# Training stage one model on full training data set (clinical features only) ----
message("Training Stage 1 model...")
s1_workflow <- workflow(s1_recipe, s1_champion_spec)
s1_fit <- fit(s1_workflow, data = train_sub)

# raw probabilities from the fitted Stage 1 model
cal_raw <- predict(s1_fit, new_data = cal_sub, type = "prob")

# standardising column names
cal_raw <- setNames(cal_raw, gsub("[- ]", "_", colnames(cal_raw)))

# combining predictions with true outcomes
cal_data_for_fit <- bind_cols(
  cal_raw,
  cal_sub %>% dplyr::select(SFI_5cat)
)

# fitting the calibrator on the prepared data frame
calibrator <- cal_estimate_isotonic(
  cal_data_for_fit,
  truth = SFI_5cat
)

find_best_thresholds <- function(fit, cal_data, calibrator) {
  
  cal_probs <- predict(fit, cal_data, type = "prob") %>%
    setNames(gsub("[- ]", "_", colnames(.)))
  
  cal_probs <- cal_probs %>%
    cal_apply(calibrator) %>%
    bind_cols(truth = cal_data$SFI_5cat, .) %>%
    # --- FIX IS HERE ---
    mutate(
      # 1. Standardize the character values in the truth column
      truth = gsub("[- ]", "_", truth),
      # 2. Create the estimate column
      .pred_class = max_prob_class(.)
    ) %>%
    # 3. Ensure both are factors with the exact same levels
    mutate(
      truth = factor(truth, levels = levels(.pred_class))
    )
  
  # The rest of the function remains the same
  threshold_grid <- expand_grid(
    thresh_non_severe = seq(0.70, 0.95, by = 0.05),
    thresh_onset_24h  = seq(0.70, 0.95, by = 0.05)
  )
  
  tune_res <- threshold_grid %>%
    mutate(
      precision = map2_dbl(thresh_non_severe, thresh_onset_24h, ~ {
        df <- cal_probs %>%
          filter(.pred_Non_Severe > .x | .pred_Onset_within_24_hours > .y)
        if (nrow(df) == 0) return(0)
        # This will now work correctly
        yardstick::accuracy(df, truth = truth, estimate = .pred_class)$.estimate
      })
    )
  
  tune_res %>%
    slice_max(order_by = precision, n = 1, with_ties = FALSE) %>%
    dplyr::select(starts_with("thresh"))
}

optimal_thresholds <- find_best_thresholds(s1_fit, cal_sub, calibrator)

message(paste("Optimal thresholds found: Non-Severe >", round(optimal_thresholds$thresh_non_severe, 2),
              "| Onset < 24h >", round(optimal_thresholds$thresh_onset_24h, 2)))



# Training the stage two meta learner ----
# Training the specialist models
message("Training Stage 2 specialist models...")
s2_specialist_fits <- map(s2_specialist_specs, ~{
  fit(workflow(recipe_all, .x), data = train_sub)
})

# Data prep for the meta learner
message("Preparing data for meta-learner...")

# predictions from specialists on the calibration set
base_preds_for_stacking <- map_dfc(s2_specialist_fits, ~{
  predict(.x, new_data = cal_sub, type = "prob")
})

meta_df <- bind_cols(base_preds_for_stacking, SFI_5cat = cal_sub$SFI_5cat)



# Defining + training meta learner
message("Training the meta-learner...")
meta_wf <- workflow(
  recipe(SFI_5cat ~ ., data = meta_df) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_pca(all_numeric_predictors(), threshold = 0.90),
  multinom_reg() %>% set_engine("nnet", trace = FALSE)
)

s2_meta_learner_fit <- fit(meta_wf, data = meta_df)

# Saving meta learner data for later reference ----
saveRDS(cal_sub, "Final Cascade Meta Learner Datasets/cal_sub.rds")
saveRDS(train_sub, "Final Cascade Meta Learner Datasets/train_sub.rds")
saveRDS(cal_split, "Final Cascade Meta Learner Datasets/cal_split.rds")
saveRDS(meta_df, "Final Cascade Meta Learner Datasets/meta_df.rds")

# Evaluating on completely held out test data ----
message("Evaluating final cascade on the test set...")

severe_classes_std <- c("Onset_greater_than_24_hours", "Onset_within_24_hours")

tic("Total Evaluation Time")

# Calibrated S1 Predictions on Test Set
s1_test_preds <- predict(s1_fit, raw_test, type = "prob") %>%
  setNames(gsub("[- ]", "_", colnames(.))) %>%
  cal_apply(calibrator) %>%
  mutate(.pred_class_s1 = max_prob_class(.))

# Uncertain Cases in Test Set
uncertain_idx_test <- which(
  !(s1_test_preds[[".pred_Non_Severe"]] > optimal_thresholds$thresh_non_severe |
      s1_test_preds[[".pred_Onset_within_24_hours"]] > optimal_thresholds$thresh_onset_24h)
)

# S2 predictions for uncertain cases
if (length(uncertain_idx_test) > 0) {
  uncertain_raw_test <- raw_test[uncertain_idx_test, ]
  
  # Predictions from specialists
  base_preds_for_s2 <- map_dfc(s2_specialist_fits, ~{
    predict(.x, new_data = uncertain_raw_test, type = "prob")
  })
  
  # Feeding to the meta-learner
  s2_preds <- predict(s2_meta_learner_fit, new_data = base_preds_for_s2, type = "class") %>%
    rename(.pred_class_s2 = .pred_class)
}

# Final predictions data frame ----
final_results_df <- raw_test %>%
  dplyr::select(SFI_5cat) %>% 
  bind_cols(s1_test_preds) %>%
  mutate(s2_invoked = row_number() %in% uncertain_idx_test)

# patching in S2 predictions
if (length(uncertain_idx_test) > 0) {
  final_results_df$.pred_class_s2 <- NA
  final_results_df$.pred_class_s2[uncertain_idx_test] <- s2_preds$.pred_class_s2
}

# Applying safeguard policy
classes_to_freeze <- c("Onset_greater_than_24_hours", "Onset_within_24_hours")

s1_levels <- levels(final_results_df$.pred_class_s1)

final_results_df <- final_results_df %>%
  mutate(
    .pred_class_s2_char = s1_levels[.pred_class_s2],
    
    .pred_class_provisional = factor(
      if_else(
        s2_invoked,
        .pred_class_s2_char,
        as.character(.pred_class_s1)
      ),
      levels = s1_levels
    ),
    
    # safeguard logic
    # Use the character version of the provisional prediction for the final step
    .pred_class_final_char = case_when(
      s2_invoked & as.character(.pred_class_s1) %in% classes_to_freeze ~ as.character(.pred_class_s1),
      
      TRUE ~ as.character(.pred_class_provisional) 
    ),
    
    .pred_class_final = factor(.pred_class_final_char, levels = s1_levels)
  )

toc()

# Saving full model ----
full_cascade_model <- list(
  stage_1_fit = s1_fit,
  calibrator = calibrator,
  optimal_thresholds = optimal_thresholds,
  s2_specialist_fits = s2_specialist_fits,
  s2_meta_learner_fit = s2_meta_learner_fit,
  severe_classes_std = severe_classes_std, 
  s1_levels = s1_levels                 
)

saveRDS(full_cascade_model, "Fully Trained Cascade Model/full_cascade_model.rds")

# Final metrics ----
message("Assembling the final comprehensive results table...")

results_for_metrics <- final_results_df %>%
  mutate(
    SFI_5cat_std = factor(gsub("[- ]", "_", SFI_5cat), levels = levels(.pred_class_final))
  )

# These are based on the final (safeguarded) predictions
main_metrics <- metric_set(bal_accuracy, sens, spec)
main_metrics_results <- results_for_metrics %>%
  main_metrics(truth = SFI_5cat_std, estimate = .pred_class_final) %>%
  dplyr::select(.metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

# Macro PR AUC
pr_auc_macro_val <- results_for_metrics %>%
  pr_auc(
    truth = SFI_5cat_std,
    .pred_Non_Severe,
    .pred_Onset_greater_than_24_hours,
    .pred_Onset_within_24_hours,
    .pred_Probable_Non_Severe,
    .pred_Probable_Severe,
    estimator = "macro"
  ) %>%
  pull(.estimate)

# Individual class pr auc
class_levels <- levels(results_for_metrics$SFI_5cat_std)

# Loop through each class to calculate its PR AUC
per_class_pr_auc <- map_dfr(class_levels, function(current_class) {
  
  prob_col_name <- paste0(".pred_", current_class)
  
  pr_auc_val <- results_for_metrics %>%
    mutate(
      truth_binary = factor(
        if_else(SFI_5cat_std == current_class, current_class, "Other")
      )
    ) %>%
    pr_auc(
      truth = truth_binary,
      !!sym(prob_col_name), 
      event_level = "first"
    )
  
  tibble(
    class = current_class,
    pr_auc = pr_auc_val$.estimate
  )
})

pr_auc_non_severe_val <- per_class_pr_auc[5,2]

# average severe sensitivity
conf_matrix <- results_for_metrics %>%
  conf_mat(truth = SFI_5cat_std, estimate = .pred_class_final)

conf_matrix_safeguard <- results_for_metrics %>%
  conf_mat(truth = SFI_5cat_std, estimate = .pred_class_provisional)

cm_table <- conf_matrix$table
per_class_sens_values <- diag(cm_table) / rowSums(cm_table)
sens_by_class <- tibble(.level = rownames(cm_table), .estimate = per_class_sens_values)
avg_severe_sens_val <- sens_by_class %>%
  filter(.level %in% severe_classes_std) %>%
  summarise(avg_sens = mean(.estimate, na.rm = TRUE)) %>%
  pull(avg_sens)

# Effective Handoff Rates (for Safeguard Policy)
effective_handoffs <- results_for_metrics %>%
  summarise(
    n_effective_handoffs = sum(s2_invoked & !(as.character(.pred_class_s1) %in% severe_classes_std)),
    n_effective_severe_handoffs = sum(s2_invoked & !(as.character(.pred_class_s1) %in% severe_classes_std) & (SFI_5cat_std %in% severe_classes_std)),
    n_total_severe = sum(SFI_5cat_std %in% severe_classes_std),
    n_total = n()
  )
effective_handoff_rate_val <- (effective_handoffs$n_effective_handoffs / effective_handoffs$n_total) * 100
effective_severe_handoff_rate_val <- (effective_handoffs$n_effective_severe_handoffs / effective_handoffs$n_total_severe) * 100

# S2 Net Corrections (With Safeguard)
s2_net_corrections_safeguard_val <- results_for_metrics %>%
  filter(s2_invoked) %>%
  mutate(
    s1_correct = (as.character(SFI_5cat_std) == as.character(.pred_class_s1)),
    s2_correct = (as.character(SFI_5cat_std) == as.character(.pred_class_final))
  ) %>%
  summarise(s2_net_benefit = sum(s2_correct & !s1_correct, na.rm = TRUE) - sum(!s2_correct & s1_correct, na.rm = TRUE)) %>%
  pull(s2_net_benefit)

# S2 Net Corrections (If Safeguard Ignored)
s2_net_corrections_standard_val <- results_for_metrics %>%
  filter(s2_invoked) %>%
  mutate(
    s1_correct = (as.character(SFI_5cat_std) == as.character(.pred_class_s1)),
    s2_correct_standard = (as.character(SFI_5cat_std) == as.character(.pred_class_provisional))
  ) %>%
  summarise(s2_net_benefit_standard = sum(s2_correct_standard & !s1_correct, na.rm = TRUE) - sum(!s2_correct_standard & s1_correct, na.rm = TRUE)) %>%
  pull(s2_net_benefit_standard)

# Handoff Rates
handoff_rates <- results_for_metrics %>%
  summarise(
    `% Handoff Rate` = (sum(s2_invoked) / n()) * 100,
    `% Severe Case Handoff Rate` = (sum(s2_invoked & SFI_5cat_std %in% severe_classes_std) /
                                      sum(SFI_5cat_std %in% severe_classes_std)) * 100
  )

pr_auc_non_severe_val <- per_class_pr_auc[5,2]
  
# Final table
final_summary_table <- tibble(
  `Stage One Model` = "Random Forest - Onset < 24h Specialist",
  `Stage Two Model` = "Meta-Leaner",
  `Policy` = "Safeguard",
  `Balanced Accuracy` = main_metrics_results$bal_accuracy,
  `Avg. Severe Sensitivity` = avg_severe_sens_val,
  `PR AUC (Non-Severe)` = pr_auc_non_severe_val,
  `Macro PR AUC` = pr_auc_macro_val,
  `Macro Sensitivity` = main_metrics_results$sens,
  `Macro Specificity` = main_metrics_results$spec,
  `% Handoff Rate` = handoff_rates$`% Handoff Rate`,
  `% Effective Handoff Rate` = effective_handoff_rate_val,
  `% Severe Case Handoff Rate` = handoff_rates$`% Severe Case Handoff Rate`,
  `% Effective Severe Case Handoff Rate` = effective_severe_handoff_rate_val,
  `S2 Net Corrections (with safeguard)` = s2_net_corrections_safeguard_val,
  `S2 Net Corrections (if safeguard ignored)` = s2_net_corrections_standard_val
)

print(final_summary_table)

saveRDS(final_summary_table, "Visualisations/Tables/final_cascade_final_table_test_data_run.rds")
saveRDS(cm_table, "Visualisations/Tables/conf_matrix_final_cascade.rds")

# SENSITIVITY ANALYSES ----
# When the model predicts a non severe class, how often is the true outcome severe? ----
# severe classes
critical_onset_classes <- c("Onset_greater_than_24_hours", "Onset_within_24_hours")

# Total number of times the model predicted "Non_Severe" (denominator is the same)
predicted_non_severe_total <- sum(cm_table[, "Non_Severe"])

# number of times it predicted "Non_Severe" but was ACTUALLY a critical "Onset" class
predicted_non_severe_but_was_onset <- sum(cm_table[critical_onset_classes, "Non_Severe"])

# Calculate this specific error rate as a percentage
critical_error_rate <- (predicted_non_severe_but_was_onset / predicted_non_severe_total) * 100

sens1 <- print(paste("When the model predicts 'Non-Severe', it is actually a critical 'Onset' class", round(critical_error_rate, 2), "% of the time."))

# When the model predicts a severe class, how often is the true outcome non severe? ----
severe_classes_std <- c("Onset_greater_than_24_hours", "Onset_within_24_hours")

# Total number of times the model predicted any severe class
predicted_severe_total <- sum(cm_table[, severe_classes_std])

# Number of times it predicted a severe class but the truth was "Non_Severe"
predicted_severe_but_was_non_severe <- sum(cm_table["Non_Severe", severe_classes_std])

# Calculate the FDR as a percentage
fdr_severe <- (predicted_severe_but_was_non_severe / predicted_severe_total) * 100

sens2 <- print(paste("When the model predicts any 'Severe' class, it is actually 'Non-Severe'", round(fdr_severe, 2), "% of the time."))


# Of all cases handed off to stage two, what percentage were actually severe? (NON-SAFEGUARD)----
# Filter for only the cases that were handed off
handed_off_cases <- results_for_metrics %>%
  filter(s2_invoked)

# total number of cases handed off
total_handoffs <- nrow(handed_off_cases)

# raw count of how many of those handoffs were truly severe
severe_handoff_count <- sum(handed_off_cases$SFI_5cat_std %in% severe_classes_std)

# percentage (hit rate)
if (total_handoffs > 0) {
  handoff_hit_rate <- (severe_handoff_count / total_handoffs) * 100
} else {
  handoff_hit_rate <- 0
}

sens3 <- print(sprintf(
  "Of all cases handed off to Stage 2, without the safeguard policy (N = %d), %.2f%% (N = %d) were actually severe.",
  total_handoffs,
  handoff_hit_rate,
  severe_handoff_count
))

# Of all cases handed off to stage two that weren't caught by the safeguard policy, what % were actually severe? (SAFEGUARD IN PLACE)----

# Define the hypothetical subgroup of cases handed off
hypothetical_handoff_cases <- results_for_metrics %>%
  filter(
    s2_invoked & !(.pred_class_s1 %in% severe_classes_std)
  )


total_hypothetical_handoffs <- nrow(hypothetical_handoff_cases)

severe_hypothetical_count <- sum(hypothetical_handoff_cases$SFI_5cat_std %in% severe_classes_std)


if (total_hypothetical_handoffs > 0) {
  hypothetical_hit_rate <- (severe_hypothetical_count / total_hypothetical_handoffs) * 100
} else {
  hypothetical_hit_rate <- 0
}

sens4 <- print(sprintf(
  "Of cases handed off that were NOT already called severe by Stage 1 (N = %d), %.2f%% (N = %d) were actually severe.",
  total_hypothetical_handoffs,
  hypothetical_hit_rate,
  severe_hypothetical_count
))

# When stage one missed a severe case and flagged it for stage two processing, how often did stage two correct the mistake? ----
# Filter for the specific group of interest: severe cases S1 missed and handed off
s2_to_the_rescue <- results_for_metrics %>%
  filter(
    s2_invoked,
    SFI_5cat_std %in% severe_classes_std,
    as.character(SFI_5cat_std) != as.character(.pred_class_s1)
  )

# total number of S1 mistakes that S2 had a chance to correct
total_opportunities_to_correct <- nrow(s2_to_the_rescue)

# how many times S2's prediction matched the true severe label
successful_corrections_count <- s2_to_the_rescue %>%
  filter(as.character(SFI_5cat_std) == as.character(.pred_class_provisional)) %>%
  nrow()

if (total_opportunities_to_correct > 0) {
  s2_correction_rate <- (successful_corrections_count / total_opportunities_to_correct) * 100
} else {
  s2_correction_rate <- 0
}

sens5 <- print(sprintf(
  "When S1 missed a severe case and handed it off (N = %d), Stage 2 corrected the mistake %.2f%% (N = %d) of the time.",
  total_opportunities_to_correct,
  s2_correction_rate,
  successful_corrections_count
))

# All sensitivity questions ----
sens1;sens2;sens3;sens4;sens5

# Subgroup performance; by site, sex, and age ----
library(tidyverse)
library(yardstick)

results_with_subgroups <- results_for_metrics %>%
  bind_cols(
    raw_test %>% dplyr::select(site, sex, age.months)
  )

# Helper function 
calculate_subgroup_metrics <- function(subgroup_df, group_name) {
  # Create confusion matrix
  conf_matrix <- subgroup_df %>%
    conf_mat(truth = SFI_5cat_std, estimate = .pred_class_final)
  
  # Calculate balanced accuracy
  bal_acc <- summary(conf_matrix) %>%
    filter(.metric == "bal_accuracy") %>%
    pull(.estimate)
  
  # Calculate avg severe sensitivity
  cm_table <- conf_matrix$table
  per_class_sens <- diag(cm_table) / rowSums(cm_table)
  sens_df <- tibble(.level = rownames(cm_table), .estimate = per_class_sens)
  
  avg_sev_sens <- sens_df %>%
    filter(.level %in% severe_classes_std) %>%
    summarise(avg_sens = mean(.estimate, na.rm = TRUE)) %>%
    pull(avg_sens)
  
  # Return results
  tibble(
    subgroup = group_name,
    bal_accuracy = bal_acc,
    avg_severe_sensitivity = avg_sev_sens
  )
}

# By Site
subgroup_by_site <- results_with_subgroups %>%
  group_by(site) %>%
  group_map(~ calculate_subgroup_metrics(.x, .y$site)) %>%
  bind_rows()

# By Sex
subgroup_by_sex <- results_with_subgroups %>%
  group_by(sex) %>%
  group_map(~ calculate_subgroup_metrics(.x, .y$sex)) %>%
  bind_rows()

# By Age
subgroup_by_age <- results_with_subgroups %>%
  mutate(age_group = if_else(age.months < 12, "Infant (<12m)", "Child (>=12m)")) %>%
  group_by(age_group) %>%
  group_map(~ calculate_subgroup_metrics(.x, .y$age_group)) %>%
  bind_rows()


# results
print("--- Subgroup Performance by Site ---")
print(subgroup_by_site)
print("--- Subgroup Performance by Sex ---")
print(subgroup_by_sex)
print("--- Subgroup Performance by Age ---")
print(subgroup_by_age)

summarize_severe_predictions <- function(data, group_var) {
  # Ensure group_var is treated as a column name
  group_var_sym <- ensym(group_var)
  
  data %>%
    # Filter for only the cases that were actually severe
    filter(SFI_5cat_std %in% severe_classes_std) %>%
    # Group by the subgroup variable (e.g., site, age_group)
    group_by(!!group_var_sym) %>%
    # Count how each of these severe cases was predicted
    count(.pred_class_final, name = "count") %>%
    # Pivot to a wide format for easy reading
    pivot_wider(
      names_from = .pred_class_final,
      values_from = count,
      values_fill = 0 # Fill in 0 for classes that weren't predicted
    ) %>%
    # Add a total column for context
    mutate(
      total_severe_cases = rowSums(across(where(is.numeric))),
      .after = 1
    )
}

# For Site
site_severe_summary <- summarize_severe_predictions(results_with_subgroups, site)

# For Age
age_severe_summary <- results_with_subgroups %>%
  mutate(age_group = if_else(age.months < 12, "Infant (<12m)", "Child (>=12m)")) %>%
  summarize_severe_predictions(age_group)

# For Sex
sex_severe_summary <- summarize_severe_predictions(results_with_subgroups, sex)

# results
print("--- Breakdown of Severe Case Predictions by Site ---")
print(site_severe_summary)

print("--- Breakdown of Severe Case Predictions by Age ---")
print(age_severe_summary)

print("--- Breakdown of Severe Case Predictions by Sex ---")
print(sex_severe_summary)

# saving 
saveRDS(site_severe_summary, "Visualisations/Tables/final_cascade_site_sensitivity_analysis.rds")
saveRDS(age_severe_summary, "Visualisations/Tables/final_cascade_age_sensitivity_analysis.rds")
saveRDS(sex_severe_summary, "Visualisations/Tables/final_cascade_sex_sensitivity_analysis.rds")

