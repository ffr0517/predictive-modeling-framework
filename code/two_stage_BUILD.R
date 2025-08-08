# Set up ----
library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(tidymodels)
library(readr)
library(tictoc)
library(probably) 

full_results_table <- read_csv("Trained Models/all_ml_model_results.csv", show_col_types = FALSE)

raw_train <- readRDS("SpotSepsis Data/Derived/RAW_train_df.rds")
folds <- readRDS("SpotSepsis Data/Derived/folds5x5.rds")

# all stage one recipes  
#recipe_clin_rf <- readRDS("SpotSepsis Data/Derived/ml_recipe_clin.rds") 
#recipe_clin_xgb <- readRDS("SpotSepsis Data/Derived/xgb_ml_recipe_clin.rds") 
#recipe_clin_svm <- readRDS("SpotSepsis Data/Derived/svm_ml_recipe_clin.rds") 
#recipe_clin_glmnet <- readRDS("SpotSepsis Data/Derived/elastic_net_ml_recipe_clin.rds")  

# stage two recipe
#recipe_all <- readRDS("SpotSepsis Data/Derived/xgb_ml_recipe_blueprint.rds") 

# Helper functions ----
# Build model spec from candidate row
build_spec_from <- function(cand) {
  model_type <- str_extract(cand$source_file, "random_forest|xgboost|svm|elastic_net")
  switch(model_type,
         random_forest = rand_forest(mode = "classification", trees = 1000,
                                     mtry = cand$mtry, min_n = cand$min_n) %>% set_engine("ranger"),
         xgboost       = boost_tree(mode = "classification", trees = 1000,
                                    tree_depth = cand$tree_depth, min_n = cand$min_n,
                                    learn_rate = cand$learn_rate) %>% set_engine("xgboost"),
         svm           = svm_rbf(mode = "classification",
                                 cost = cand$cost, rbf_sigma = cand$rbf_sigma) %>% set_engine("kernlab"),
         elastic_net   = multinom_reg(mode = "classification",
                                      penalty = cand$penalty, mixture = cand$mixture) %>% set_engine("glmnet")
  )
}

# Identify uncertain indices given probabilities and thresholds
find_uncertain <- function(prob_df, thresholds) {
  which(
    !(prob_df[[".pred_Non_Severe"]] > thresholds$thresh_non_severe |
        prob_df[[".pred_Onset_within_24_hours"]] > thresholds$thresh_onset_24h)
  )
}

# Determine class with max probability per row
max_prob_class <- function(df) {
  prob_cols <- df %>% dplyr::select(starts_with(".pred_"))
  levels   <- gsub("\\.pred_", "", colnames(prob_cols))
  factor(levels[max.col(prob_cols)], levels = levels)
}

# Merge Stage1 & Stage2 predictions
patch_predictions <- function(base_df, idx_unc, s2_probs, s2_class, levels) {
  df <- base_df
  df$.pred_class_final <- df$.pred_class
  if (length(idx_unc) > 0) {
    df$.pred_class_final[idx_unc] <- s2_class
    prob_names <- colnames(s2_probs)
    df[idx_unc, prob_names] <- s2_probs
  }
  df %>% arrange(.row)
}


# Finding 'stage two' models ----
# Calculate the overall macro PR AUC for all models
results_with_macros <- full_results_table %>%
  rowwise() %>%
  mutate(
    manual_macro_avg_pr_auc = mean(c_across(starts_with("pr_auc_")), na.rm = TRUE)
  ) %>%
  ungroup()

champion_s2_model <- results_with_macros %>%
  filter(str_detect(source_file, "all_feature")) %>%
  slice_max(order_by = manual_macro_avg_pr_auc, n = 1, with_ties = FALSE)

champion_xgb_spec <- boost_tree(
  mode = "classification",
  trees = 1000,
  mtry = champion_s2_model$mtry,
  min_n = champion_s2_model$min_n,
  tree_depth = champion_s2_model$tree_depth,
  learn_rate = champion_s2_model$learn_rate,
  sample_size = champion_s2_model$sample_size
) %>%
  set_engine("xgboost")

# Finding specialist models for stage two 'meta' model ----
outcome_classes <- c(
  "Probable Severe",
  #"Probable Non-Severe", omitted as is the same model as probable severe
  "Onset greater than 24 hours",
  "Non-Severe",
  "Onset within 24 hours"
)

specialist_candidates <- map(outcome_classes, function(oc) {
  full_results_table %>%
    filter(str_detect(source_file, "all_feature")) %>%
    slice_max(order_by = .data[[ paste0("pr_auc_", oc) ]], n = 1, with_ties = FALSE)
})

names(specialist_candidates) <- outcome_classes

champion_xgb_spec <- boost_tree(
  mode = "classification",
  trees = 1000,
  mtry = champion_s2_model$mtry,
  min_n = champion_s2_model$min_n,
  tree_depth = champion_s2_model$tree_depth,
  learn_rate = champion_s2_model$learn_rate,
  sample_size = champion_s2_model$sample_size
) %>%
  set_engine("xgboost")

# Finding 'stage one' clinical only models (identified in all model evaluation script) ----
clinical_models_to_test <- full_results_table %>%
  filter(str_detect(source_file, "clin")) %>%
  mutate(
    avg_severe_sensitivity = (`sens_Onset greater than 24 hours` + `sens_Onset within 24 hours`) / 2
  )

# performance matrix for the domination check
perf_mat <- clinical_models_to_test %>%
  dplyr::select(avg_severe_sensitivity, `pr_auc_Non-Severe`) %>%
  as.matrix()

# check for domination
is_dominated <- function(i, mat) {
  any(apply(mat, 1, function(other) {
    all(other >= mat[i, ]) && any(other > mat[i, ])
  }))
}

# all non-dominated models to create the Pareto front
pareto_mask <- vapply(seq_len(nrow(perf_mat)),
                      FUN = function(i) !is_dominated(i, perf_mat),
                      FUN.VALUE = logical(1))

pareto_models <- clinical_models_to_test[pareto_mask, ]

# simplicity tie-breakers to get a unique model for each performance profile
pareto_unique <- pareto_models %>%
  group_by(avg_severe_sensitivity, `pr_auc_Non-Severe`) %>%
  { if ("min_n" %in% names(.)) slice_max(., order_by = min_n, n = 1, with_ties = FALSE) else . } %>%
  { if ("mtry" %in% names(.)) slice_min(., order_by = mtry, n = 1, with_ties = FALSE) else . } %>%
  ungroup()

# Rank the Pareto models from most efficient to most safe
df_ranked <- pareto_unique %>%
  arrange(avg_severe_sensitivity) %>%
  mutate(rank_pct = (row_number() - 1) / (n() - 1))

targets <- c(0, 0.25, 0.5, 0.75, 1)

pareto_candidates <- map_dfr(targets, function(tgt) {
  df_ranked %>%
    mutate(diff = abs(rank_pct - tgt)) %>%
    slice_min(diff, with_ties = FALSE)
})


# Finding high performing random forest model (max pr auc on within 24 hours outcome) ----
specialist_rf_candidate <- full_results_table %>%
  filter(str_detect(source_file, "clin_random_forest")) %>%
  slice_max(order_by = `pr_auc_Onset within 24 hours`, n = 1, with_ties = FALSE)


# The six stage one candidates ----
all_possible_params <- c(".config", "source_file", "cost_complexity", "tree_depth", "min_n", "mtry", 
                         "learn_rate", "sample_size", "cost", "rbf_sigma", "penalty", "mixture")

candidate_s1_models <- bind_rows(
  pareto_candidates,
  specialist_rf_candidate
) %>%
  distinct(.config, source_file, .keep_all = TRUE) %>%
  dplyr::select(any_of(all_possible_params))

# Saving models ----
saveRDS(champion_s2_model, "Trained Models/Two Stage/champion_s2_model.rds") # Best all feature model
saveRDS(specialist_candidates, "Trained Models/Two Stage/stacked_specialists_s2_models.rds") # The voters for the 'meta' model
saveRDS(candidate_s1_models, "Trained Models/Two Stage/candidate_s1_models.rds") # The 5 pareto models, and the champion within 24 hours model as it used clinical only feature set

# Handoff rule and constraints ----
initial_handoff_thresholds <- list(
  "Non-Severe" = 0.90,
  "Onset within 24 hours" = 0.70
)

# Saving rules/constraints ----
saveRDS(initial_handoff_thresholds, "Trained Models/Two Stage/initial_handoff_thresholds.rds")

# Initialising results and operational log lists ----
final_predictions <- list()

operational_logs <- list()

meta_learner_logs <- list()

pca_logs <- list() 


 
# parallel compute ----
plan(cluster, workers = parallel::detectCores() - 1)

plan()

# Running full loop ----
tic() 

for (i in seq_len(nrow(candidate_s1_models))) {
  s1_candidate <- candidate_s1_models[i, ]
  for (j in seq_len(nrow(folds))) {
    fold    <- folds$splits[[j]]
    fold_id <- folds$id[j]
    task_base <- paste0(s1_candidate$source_file, "_",
                        s1_candidate$.config, "_",
                        "s1_model_", j, "_", 
                        fold_id)
    message("Processing: ", task_base)
    
    # 1. Data splits & relevel
    analysis_set   <- analysis(fold)
    assessment_set <- assessment(fold)
    set.seed((i-1)*100 + j)
    inner_split <- initial_split(analysis_set, prop = 0.8, strata = SFI_5cat)
    train_sub <- training(inner_split)
    cal_sub   <- testing(inner_split)
    full_lvls <- levels(analysis_set$SFI_5cat)
    outcome_levels_vector <- full_lvls
    train_sub$SFI_5cat <- factor(train_sub$SFI_5cat, levels = full_lvls)
    cal_sub$SFI_5cat   <- factor(cal_sub$SFI_5cat,   levels = full_lvls)
    assessment_set$SFI_5cat <- factor(assessment_set$SFI_5cat, levels = full_lvls)
    
    # 2. Stage-1 recipe & spec
    clin_predictors <- c(
      "site", "age.months", "sex", "bgcombyn", "adm.recent", "wfaz",
      "waste", "stunt", "cidysymp", "prior.care", "travel.time.bin",
      "urti", "lrti", "diarrhoeal", "neuro", "auf", "syndrome.resp",
      "syndrome.nonresp", "pneumo", "sev.pneumo", "ensapro", "vomit.all",
      "seiz", "pfacleth", "not.alert", "danger.sign", "hr.all", "rr.all",
      "oxy.ra", "envhtemp", "crt.long", "LqSOFA", "parenteral_screen", "SIRS"
    )
    
    all_predictors_list <- c(
      clin_predictors, "weight", "ANG1", "ANG2", "CHI3L", "CRP", "CXCl10", 
      "IL1ra", "IL6", "IL8", "IL10", "PROC", "TNFR1", "STREM1", "VEGFR1", 
      "supar", "lblac", "lbglu", "enescbchb1"
    )
    
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
    s1_spec <- build_spec_from(s1_candidate)
    s1_fit  <- fit(workflow(s1_recipe, s1_spec), data = train_sub)
    
    # 3. Calibration and threshold tuning
    cal_raw <- predict(s1_fit, new_data = cal_sub, type = "prob")
    
    # 3a. Normalize hyphens & spaces to underscores
    cal_raw <- setNames(cal_raw, gsub("[- ]", "_", colnames(cal_raw)))
    
    # 3b. Fit the calibrator
    calibrator <- cal_estimate_isotonic(
      bind_cols(cal_raw, cal_sub %>% dplyr::select(SFI_5cat)),
      truth = SFI_5cat
    )
    
    #3c. Apply calibrator
    cal_probs <- cal_apply(cal_raw, calibrator)
    
    # 3d. Rebind the true outcome
    cal_probs <- bind_cols(
      cal_probs,
      SFI_5cat = cal_sub$SFI_5cat
    )
    
    # 3e. Add the hard prediction
    cal_probs$.pred_class <- max_prob_class(cal_probs)
    
    cal_probs <- cal_probs %>%
      mutate(
        SFI_5cat = gsub("[- ]", "_", SFI_5cat),
        SFI_5cat = factor(SFI_5cat, levels = levels(.pred_class))
      )
    
    # 3f. Now run your threshold tuning over this enriched cal_probs
    threshold_grid <- expand_grid(
      thresh_non_severe = seq(0.70, 0.95, by = 0.05),
      thresh_onset_24h  = seq(0.70, 0.95, by = 0.05)
    )
    
    tune_res <- threshold_grid %>%
      mutate(
        coverage = map2_dbl(thresh_non_severe, thresh_onset_24h, ~ {
          sum(cal_probs$.pred_Non_Severe > .x |
                cal_probs$.pred_Onset_within_24_hours > .y) /
            nrow(cal_probs)
        }),
        precision = map2_dbl(thresh_non_severe, thresh_onset_24h, ~ {
          df <- cal_probs %>%
            filter(.pred_Non_Severe > .x |
                     .pred_Onset_within_24_hours > .y)
          if (nrow(df) == 0) return(0)
          yardstick::accuracy(df, truth = SFI_5cat, estimate = .pred_class)$.estimate
        })
      )
    
    best_row <- tune_res %>% 
      slice_max(order_by = precision, n = 1, with_ties = FALSE)
    if (nrow(best_row) == 0 || best_row$precision < 0.5) { 
      message("  -> Warning: Could not find optimal threshold. Using defaults.")
      optimal_thresholds <- tibble(thresh_non_severe = 0.95, thresh_onset_24h = 0.95)
    } else {
      optimal_thresholds <- best_row %>%
        dplyr::select(thresh_non_severe, thresh_onset_24h)
    }
    
    # 4. Identify uncertain in train_sub
    train_cal <- cal_apply(
      predict(s1_fit, train_sub, type = "prob") %>%
        setNames(gsub("[- ]", "_", colnames(.))),
      calibrator
    )
    uncertain_idx_train <- find_uncertain(train_cal, optimal_thresholds)
    uncertain_train <- train_sub[uncertain_idx_train, ]
    
    # 5. Build recipe_all on train_sub
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
    
    # 6. Prepare Stage-2 variants
    # List the two S2 strategies you want to run
    s2_variants <- list(
      champion_xgb = NULL,
      stacker = NULL
    )
    
    # 7. Fit specialist models (only needed for the 'stacker' variant)
    model_fits_specialist <- map(outcome_classes, function(oc) {
      
      # Get the hyperparameter set for the current outcome class
      cand <- specialist_candidates[[oc]]
      
      spec <- NULL
      
      # Conditionally build the model spec based on the source_file name
      if (grepl("xgboost", cand$source_file, fixed = TRUE)) {
        
        # Build the XGBoost specification
        spec <- boost_tree(
          mode = "classification", 
          trees = 1000,
          tree_depth = cand$tree_depth,
          min_n      = cand$min_n,
          learn_rate = cand$learn_rate,
          mtry       = cand$mtry,
          sample_size = cand$sample_size
        ) %>%
          set_engine("xgboost")
        
      } else if (grepl("random_forest", cand$source_file, fixed = TRUE)) {
        
        # Build the Random Forest specification
        spec <- rand_forest(
          mode = "classification", 
          trees = 1000,
          min_n = cand$min_n,
          mtry  = cand$mtry
        ) %>%
          set_engine("ranger")
      }
      
      # Fit the workflow using the correctly defined 'spec'
      if (!is.null(spec)) {
        fit(workflow(recipe_all, spec), data = train_sub)
      }
    })
    names(model_fits_specialist) <- outcome_classes
    
    # 8. Loop over each S2 strategy
    for (variant in names(s2_variants)) {
      task_id <- paste0(task_base, "_", variant)
      message("  → Variant: ", variant)
      
      # --- Fit the chosen Stage-2 model ---
      s2_fit <- NULL # Reset s2_fit for each variant
      if (variant == "champion_xgb") {
        # FIT THE CHAMPION XGB MODEL
        # This is fit on the 'uncertain' data using the pre-defined spec
        if (nrow(uncertain_train) > 5) { # Check for minimum data
          s2_fit <- fit(workflow(recipe_all, champion_xgb_spec), data = uncertain_train)
        } else {
          message("   -> Skipping champion_xgb: Not enough uncertain training data.")
        }
        
      } else { # This handles the "stacker" variant
        
        valid_specialists <- purrr::compact(model_fits_specialist)
        
        if (length(valid_specialists) > 0) {
          base_preds_for_stacking <- map_dfc(valid_specialists, function(mf) {
            predict(mf, new_data = cal_sub, type = "prob")
          })
          
          meta_df <- bind_cols(base_preds_for_stacking, SFI_5cat = cal_sub$SFI_5cat)
          
          predictor_names <- setdiff(colnames(meta_df), "SFI_5cat")
          formula_rhs <- paste0("`", predictor_names, "`", collapse = " + ")
          meta_formula <- as.formula(paste("SFI_5cat ~", formula_rhs))
          
          # --- FINAL STACKER WORKFLOW WITH PCA ---
          meta_wf <- workflow(
            recipe(meta_formula, data = meta_df) %>% 
              step_normalize(all_numeric_predictors()) %>%
              step_pca(all_numeric_predictors(), threshold = 0.90), # PCA Step
            multinom_reg() %>% set_engine("nnet", trace = FALSE)
          )
          s2_fit <- fit(meta_wf, data = meta_df)
          
          # --- LOGGING FOR META-LEARNER AND PCA ---
          if(!is.null(s2_fit)) {
            # Log the meta-learner coefficients for PC1, PC2, etc.
            raw_nnet_model <- s2_fit %>% extract_fit_engine()
            coeffs_matrix  <- coef(raw_nnet_model)
            coeffs_tidy <- as_tibble(coeffs_matrix, rownames = "y.level") %>%
              pivot_longer(cols = -y.level, names_to = "term", values_to = "estimate")
            meta_learner_logs[[task_id]] <- coeffs_tidy
            
            # Log the PCA loadings to interpret the components
            prepped_recipe <- prep(meta_wf$pre$actions$recipe$recipe)
            pca_loadings <- tidy(prepped_recipe, number = 2) # Assumes pca is 2nd step
            pca_logs[[task_id]] <- pca_loadings
          }
        } else {
          message("   -> Skipping stacker: No valid specialist models were trained.")
        }
      }
      
      # 9. Final scoring on assessment_set
      assess_cal <- cal_apply(
        predict(s1_fit, assessment_set, type = "prob") %>%
          setNames(gsub("[- ]", "_", colnames(.))),
        calibrator
      )
      assess_df  <- bind_cols(assessment_set, assess_cal) %>%
        mutate(.pred_class = max_prob_class(.), .row = row_number())
      
      uncertain_idx_assess <- find_uncertain(assess_cal, optimal_thresholds)
      
      # 10. Generate S2 predictions and patch them into the S1 results
      final_df <- assess_df # Start with S1 predictions as the default
      
      final_df$.pred_class_final <- final_df$.pred_class
      
      # Only proceed if the S2 model was successfully fit AND there are uncertain cases
      if (!is.null(s2_fit) && length(uncertain_idx_assess) > 0) {
        
        # --- A. Get predictions from the fitted S2 model ---
        uncertain_subset <- assessment_set[uncertain_idx_assess, ]
        
        s2_probs <- NULL
        s2_class <- NULL
        
        if (variant == "champion_xgb") {
          uncertain_subset$SFI_5cat <- factor(uncertain_subset$SFI_5cat, levels = full_lvls)
          s2_probs <- predict(s2_fit, uncertain_subset, type = "prob")
          s2_class <- predict(s2_fit, uncertain_subset, type = "class")$.pred_class
        } else { # Stacker
          valid_specialists <- purrr::compact(model_fits_specialist)
          base_probs <- map_dfc(valid_specialists, function(mf) {
            predict(mf, uncertain_subset, type = "prob")
          })
          s2_probs <- predict(s2_fit, base_probs, type = "prob")
          s2_class <- predict(s2_fit, base_probs, type = "class")$.pred_class
        }
        
        # --- B. Standardize S2 column names to prevent redundancy ---
        s2_class <- gsub("[- ]", "_", s2_class)
        s2_probs_std <- s2_probs %>%
          setNames(gsub("[- ]", "_", colnames(.)))
        
        # --- C. Safely patch the predictions ---
        # Overwrite the uncertain rows with the S2 class prediction
        final_df$.pred_class_final[uncertain_idx_assess] <- s2_class
        
        # Overwrite the probability values in the original columns
        prob_colnames <- colnames(s2_probs_std)
        final_df[uncertain_idx_assess, prob_colnames] <- s2_probs_std
      }
      
      # 1. Define your severe classes with standardized names
      severe_classes <- c("Onset_greater_than_24_hours", "Onset_within_24_hours")
      
      # 2. Create a flag indicating if S2 was invoked
      final_df$s2_invoked <- FALSE
      final_df$s2_invoked[uncertain_idx_assess] <- TRUE
      
      # 3. Add the new safeguard columns
      final_df <- final_df %>%
        mutate(
          # Make sure the prediction columns are standardized characters before comparing
          .pred_class_char = gsub("[- ]", "_", .pred_class),
          .pred_class_final_char = gsub("[- ]", "_", .pred_class_final),
          
          # A. Flag when the safeguard rule was triggered
          s1_safeguard_applied = (s2_invoked == TRUE & .pred_class_char %in% severe_classes),
          
          # B. Create the NEW final prediction column using the standardized character versions
          .pred_class_final_safeguard_char = case_when(
            s1_safeguard_applied == TRUE ~ .pred_class_char,
            TRUE                         ~ .pred_class_final_char
          ),
          
          # C. NOW, convert the final, clean character column to a factor
          .pred_class_final_safeguard = factor(.pred_class_final_safeguard_char, levels = levels(SFI_5cat))
        ) %>%
        # Clean up the intermediate character columns
        dplyr::select(-.pred_class_char, -.pred_class_final_char, -.pred_class_final_safeguard_char)
      
      # 11. Store results
      final_predictions[[task_id]] <- final_df
      operational_logs[[task_id]]  <- tibble(
        task_id  = task_id,
        variant  = variant,
        fold_id  = fold_id,
        coverage = 1 - (length(uncertain_idx_assess) / nrow(assessment_set)),
        optimal_thresh_non_severe = optimal_thresholds$thresh_non_severe,
        optimal_thresh_onset_24h  = optimal_thresholds$thresh_onset_24h
      )
      
    } # End S2 variant loop
  }
}

test_time_info <- toc()

saveRDS(final_predictions, "Trained Models/Two Stage/two_stage_full_run_final_predictions.rds")
saveRDS(operational_logs, "Trained Models/Two Stage/two_stage_full_run_operational_logs.rds")
saveRDS(meta_learner_logs, "Trained Models/Two Stage/two_stage_meta_learner_logs.rds")
saveRDS(pca_logs, "Trained Models/Two Stage/two_stage_pca_logs.rds")