# Set up ----
library(tidyverse) 
library(yardstick)
library(kableExtra)

# All prediction files ----
all_preds_files <- list.files(
  path = "Trained Models/",
  pattern = "_all_preds.csv$",
  recursive = TRUE,
  full.names = TRUE
)
print(all_preds_files)

# Master results data ----
outcome_levels_vector <- c(
  "Non-Severe", "Probable Non-Severe", "Probable Severe",
  "Onset greater than 24 hours", "Onset within 24 hours"
)

col_name_map <- list(
  "Non-Severe" = ".pred_Non-Severe",
  "Probable Non-Severe" = ".pred_Probable Non-Severe",
  "Probable Severe" = ".pred_Probable Severe",
  "Onset greater than 24 hours" = ".pred_Onset greater than 24 hours",
  "Onset within 24 hours" = ".pred_Onset within 24 hours"
)

all_results <- list()
outcome_levels_vector <- names(col_name_map) 

# list of all unique model configurations from the first file 
all_configs <- read_csv(all_preds_files[1], show_col_types = FALSE) %>%
  distinct(.config)

for (file_path in all_preds_files) {
  
  message(paste("Processing:", file_path))
  
  # a. Load the full prediction data for this model type
  preds_data <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      SFI_5cat    = factor(SFI_5cat,    levels = outcome_levels_vector),
      .pred_class = factor(.pred_class, levels = outcome_levels_vector)
    )
  
  # b. Loop through each model configuration within the file
  for (config_id in all_configs$.config) {
    
    # Get the predictions for just this one model
    current_model_preds <- preds_data %>% filter(.config == config_id)
    
    # c. Loop through each class to calculate its metrics
    per_class_results <- lapply(outcome_levels_vector, function(class_level) {
      
      # 1) PR AUC
      prob_col <- col_name_map[[class_level]]
      pr_auc_val <- pr_auc_vec(
        truth    = factor(ifelse(current_model_preds$SFI_5cat == class_level, "event", "other"),
                          levels = c("event","other")),
        estimate = current_model_preds[[prob_col]]
      )
      
      # 1b) ROC AUC
      roc_auc_val <- roc_auc_vec(
        truth    = factor(ifelse(current_model_preds$SFI_5cat == class_level, "event", "other"),
                          levels = c("event","other")),
        estimate = current_model_preds[[prob_col]]
      )
      
      # 2) Binary truth & class prediction
      truth_bin <- factor(ifelse(current_model_preds$SFI_5cat   == class_level, "event", "other"),
                          levels = c("event","other"))
      pred_bin  <- factor(ifelse(current_model_preds$.pred_class == class_level, "event", "other"),
                          levels = c("event","other"))
      
      # 3) Sensitivity & Specificity
      sens_val <- sens_vec(truth = truth_bin, estimate = pred_bin)
      spec_val <- spec_vec(truth = truth_bin, estimate = pred_bin)
      
      # Return a tibble for this class
      tibble(
        class    = class_level,
        pr_auc   = pr_auc_val,
        roc_auc  = roc_auc_val,
        sens     = sens_val,
        spec     = spec_val
      )
    })
    
    # d. Combine & reshape, then compute macros / averages
    model_summary <- bind_rows(per_class_results) %>%
      pivot_wider(
        names_from  = class,
        values_from = c(pr_auc, roc_auc, sens, spec),
        names_sep   = "_"
      ) %>%
      # Compute the four summary columns
      rowwise() %>%
      mutate(
        avg_pr_auc  = mean(c_across(starts_with("pr_auc_")),  na.rm = TRUE),
        avg_roc_auc = mean(c_across(starts_with("roc_auc_")), na.rm = TRUE),
        macro_sens  = mean(c_across(starts_with("sens_")),    na.rm = TRUE),
        macro_spec  = mean(c_across(starts_with("spec_")),    na.rm = TRUE)
      ) %>%
      ungroup() %>%
      # Add identifying columns
      mutate(
        .config     = config_id,
        source_file = basename(file_path)
      )
    
    # Store
    all_results[[paste(file_path, config_id)]] <- model_summary
  }
}

master_detailed_leaderboard <- bind_rows(all_results)

all_possible_params <- c(
  "cost_complexity", "tree_depth", "min_n", "mtry", 
  "learn_rate", "sample_size", "cost", "rbf_sigma", "penalty", "mixture"
)

# mapping from .config to the parameters
all_hyperparams <- map_dfr(all_preds_files, function(file_path) {
  read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::select(.config, any_of(all_possible_params)) %>%
    distinct() %>%
    mutate(source_file = basename(file_path))
})

# joining hyperparameters to results
full_results_table <- left_join(master_detailed_leaderboard, all_hyperparams, by = c(".config", "source_file"))

write.csv(full_results_table, "Trained Models/all_ml_model_results.csv")

# Finding 'Stage one' models for comparison via pareto front comparison (saftey versus efficiency)----
# clinical-only models
clinical_models_to_test <- full_results_table %>%
  filter(str_detect(source_file, "clin")) %>%
  mutate(
    # Create the average sensitivity across the three severe classes
    avg_severe_sensitivity = (`sens_Onset greater than 24 hours` + `sens_Onset within 24 hours`) / 2
  )

perf_mat <- clinical_models_to_test %>%
  dplyr::select(avg_severe_sensitivity, `pr_auc_Non-Severe`) %>%
  as.matrix()

# Function: is row i dominated by any other row?
is_dominated <- function(i, mat) {
  any(apply(mat, 1, function(other) {
    all(other >= mat[i, ]) && any(other > mat[i, ])
  }))
}

# Build logical mask: TRUE if NOT dominated
pareto_mask <- vapply(seq_len(nrow(perf_mat)),
                      FUN = function(i) !is_dominated(i, perf_mat),
                      FUN.VALUE = logical(1))

pareto_models <- clinical_models_to_test[pareto_mask, ]

pareto_unique <- pareto_models %>%
  group_by(avg_severe_sensitivity, `pr_auc_Non-Severe`) %>%
  slice_min(tree_depth,  with_ties = FALSE) %>%
  slice_max(cost_complexity, with_ties = FALSE) %>%
  slice_min(mtry, with_ties = FALSE) %>%
  ungroup()

df_ranked <- pareto_unique %>%
  arrange(avg_severe_sensitivity) %>%
  mutate(rank_pct = (row_number() - 1) / (n() - 1))

# plotting pareto curve/front ----
library(stringr)
library(ggrepel) 

# data for plot
df_ranked_plot <- df_ranked %>%
  mutate(
    model_label = case_when(
      row_number() == 1 ~ "Max Efficiency Profile (XGBoost)",
      row_number() == 2 ~ "",
      row_number() == 3 ~ "",
      row_number() == 4 ~ "",
      row_number() == 5 ~ "Max Safety Profile (SVM)",
      TRUE ~ as.character(row_number()) # A fallback in case there are more rows
    ),
    # We still create the model_name column for coloring the points
    model_name = str_to_title(str_replace(source_file, "_", " ")) %>% 
      str_extract("^[A-Za-z]+")
  )

saveRDS(df_ranked_plot, "Visualisations/df_ranked_plot_data.rds")

plot_data_full <- clinical_models_to_test %>%
  mutate(
    model_name = str_to_title(str_replace(source_file, "_", " ")) %>% 
      str_extract("^[A-Za-z]+")
  )

library(showtext) 

font_add_google("Roboto")
showtext_auto()

tradeoff_p <- ggplot() + 
  geom_point(
    data = plot_data_full, 
    aes(x = `pr_auc_Non-Severe`, y = avg_severe_sensitivity),
    color = "grey80", 
    alpha = 0.5, 
    size = 2
  ) +
  
  stat_ellipse(
    data = plot_data_full %>% filter(source_file != "clin_random_forest_all_preds.csv"),
    aes(x = `pr_auc_Non-Severe`, y = avg_severe_sensitivity, group = source_file, fill = source_file),
    type = "t",
    geom = "polygon",
    color = "black", 
    linetype = "dashed",
    alpha = 0.6,
    show.legend = FALSE
  ) +
  
  stat_ellipse(
    data = plot_data_full %>% filter(source_file == "clin_random_forest_all_preds.csv"),
    aes(x = `pr_auc_Non-Severe`, y = avg_severe_sensitivity, group = source_file, fill = source_file),
    type = "t",
    geom = "polygon",
    color = "black", 
    linetype = "dashed",
    alpha = 0.6,
    show.legend = FALSE
  ) +
  
  geom_line(
    data = df_ranked_plot,
    aes(x = `pr_auc_Non-Severe`, y = avg_severe_sensitivity),
    color = "black", 
    linewidth = 1
  ) +
  
  geom_point(
    data = df_ranked_plot,
    aes(x = `pr_auc_Non-Severe`, y = avg_severe_sensitivity, fill = model_name),
    shape = 21,
    color = "black",
    size = 5,
    stroke = 1,
    show.legend = FALSE
  ) +
  
  geom_text_repel(
    data = df_ranked_plot,
    aes(x = `pr_auc_Non-Severe`, y = avg_severe_sensitivity, label = model_label),
    size = 3.5,
    nudge_y = 0.02,
    direction = "x",
    max.overlaps = 15,
    box.padding = 0.5
  ) +
  
  labs(
    title = "Trade-off Between Safety and Efficiency for Clinical Models",
    subtitle = "Ellipses Represent 95% Confidence by Modelling Approach",
    x = "PR AUC for 'Non-Severe' Class (Screening Efficiency)",
    y = "Average Sensitivity for Severe Classes (Patient Safety)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(colour = 'black', size = 17.5, face = 'bold', family = 'Roboto'),
    plot.subtitle = element_text(colour = 'grey30', size = 12.5, face = 'italic', family = 'Roboto'),
    axis.title = element_text(family = 'Roboto', size = 14),
    axis.text = element_text(size = 12, family = 'Roboto', colour = 'black'),
    legend.position = "none"
  )

ggsave("clinical_trade_off_plot.svg", plot = tradeoff_p, width = 10, height = 5.5, units = "in", dpi = 300)

pareto_front_summary_table <- df_ranked_plot %>%
  
  mutate(
    model_name = str_extract(source_file, "decision_tree|random_forest|xgboost|svm|elastic_net")
  ) %>%
  

  rowwise() %>%
  mutate(
    Hyperparameters = case_when(
      model_name == "decision_tree" ~ paste("Complexity Penalty=", format(cost_complexity, scientific = TRUE, digits = 1),
                                            "; Tree Depth=", tree_depth,
                                            "; Min. Node Size=", min_n, sep = ""),
      model_name == "random_forest" ~ paste("Predictors per Split=", mtry, "; Min. Node Size=", min_n, sep = ""),
      model_name == "xgboost" ~ paste("Tree Depth=", tree_depth, "; Min. Node Size=", min_n, 
                                      "; Learning Rate=", round(learn_rate, 3), sep = ""),
      model_name == "svm" ~ paste("Cost=", round(cost, 1), "; Sigma=", round(rbf_sigma, 3), sep = ""),
      model_name == "elastic_net" ~ paste("Penalty=", format(penalty, scientific = TRUE, digits = 2), 
                                          "; Mixture Ratio=", round(mixture, 2), sep = ""),
      TRUE ~ "NA"
    )
  ) %>%
  ungroup() %>%
  
  mutate(
    `ML Model Type` = case_when(
      model_name == "xgboost"         ~ "Gradient Boosted (XGBoost)",
      model_name == "random_forest"   ~ "Random Forest",
      model_name == "decision_tree"   ~ "Decision Tree",
      model_name == "svm"             ~ "Support Vector Machine",
      model_name == "elastic_net"     ~ "Penalised Regression",
      TRUE ~ model_name
    )
  ) %>%
  
  dplyr::select(
    `ML Model Type`,
    `Non-Severe PR AUC` = `pr_auc_Non-Severe`,
    `Mean 'Severe' Outcome Sensitivity` = avg_severe_sensitivity,
    `Sens (Probable Severe)` = `sens_Probable Severe`,
    `Sens (Onset > 24h)` = `sens_Onset greater than 24 hours`,
    `Sens (Onset < 24h)` = `sens_Onset within 24 hours`,
    `Hyperparameters`
  ) %>%
  
  arrange(desc(`Mean 'Severe' Outcome Sensitivity`))

saveRDS(pareto_front_summary_table, "Visualisations/Tables/Pareto_Table_raw.rds")

# overall champion leaderboard ----
baseline_model_results <- read_csv("Trained Models/DAG_informed/baseline_model_final_metrics.csv", show_col_types = FALSE)

top_ml_champions <- full_results_table %>%
  rowwise() %>%
  mutate(
    manual_macro_avg_pr_auc = mean(c_across(starts_with("pr_auc_")), na.rm = TRUE),
    roc_auc_macro = mean(c_across(starts_with("roc_auc_")), na.rm = TRUE),
    sens_macro = mean(c_across(starts_with("sens_")), na.rm = TRUE),
    spec_macro = mean(c_across(starts_with("spec_")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    feature_set_clean = case_when(
      str_detect(source_file, "all_feature") ~ "All",
      str_detect(source_file, "clin")        ~ "Clinical",
      str_detect(source_file, "bio|lab")     ~ "Laboratory"
    ),
    model_name = str_extract(source_file, "decision_tree|random_forest|xgboost|svm|elastic_net")
  ) %>%
  group_by(feature_set_clean) %>%
  slice_max(order_by = manual_macro_avg_pr_auc, n = 1, with_ties = FALSE) %>%
  ungroup()

final_leaderboard_data <- bind_rows(top_ml_champions, baseline_model_results) %>%
  
  mutate(
    model_name = case_when(
      !is.na(source_file) ~ str_extract(source_file, "decision_tree|random_forest|xgboost|svm|elastic_net"),
      TRUE ~ "Logistic Regression"
    )
  ) %>%
  
  # Step 2: Create the final 'Model' display column
  mutate(
    Model = case_when(
      model_name == "xgboost"         ~ "Gradient Boosted (XGBoost)",
      model_name == "random_forest"   ~ "Random Forest",
      model_name == "decision_tree"   ~ "Decision Tree",
      model_name == "svm"             ~ "Support Vector Machine",
      model_name == "elastic_net"     ~ "Penalised Regression",
      model_name == "Logistic Regression" ~ "Logistic Regression (Baseline)",
      TRUE ~ model_name
    )
  ) %>%
  
  # Step 3: Now that 'Model' exists, create the 'Hyperparameters' column
  rowwise() %>%
  mutate(
    Hyperparameters = case_when(
      model_name == "decision_tree" ~ paste("Complexity Penalty=", format(cost_complexity, scientific = TRUE, digits = 1),
                                            "; Tree Depth=", tree_depth,
                                            "; Min. Node Size=", min_n, sep = ""),
      model_name == "random_forest" ~ paste("Predictors per Split=", mtry, "; Min. Node Size=", min_n, sep = ""),
      model_name == "xgboost" ~ paste("Tree Depth=", tree_depth, "; Min. Node Size=", min_n, 
                                      "; Learning Rate=", round(learn_rate, 3), sep = ""),
      model_name == "svm" ~ paste("Cost=", round(cost, 1), "; Sigma=", round(rbf_sigma, 3), sep = ""),
      model_name == "elastic_net" ~ paste("Penalty=", format(penalty, scientific = TRUE, digits = 2), 
                                          "; Mixture Ratio=", round(mixture, 2), sep = ""),
      Model == "Logistic Regression (Baseline)" ~ "N/A (Not Tuned)",
      TRUE ~ "NA"
    )
  ) %>%
  ungroup() %>%
  
  # Step 4: Create the final 'Feature Set' display column
  mutate(
    `Feature Set` = case_when(
      feature_set_clean == "All" ~ "All",
      feature_set_clean == "Clinical" ~ "Clinical only",
      feature_set_clean == "Laboratory" ~ "Laboratory only",
      feature_set == "Tier 1 Literature" ~ "Features Present in >= 10% of Literature",
      TRUE ~ feature_set
    )
  ) %>%
  
  # Step 5: dplyr::select the final columns for your table
  dplyr::select(
    Model,
    `Feature Set`,
    `Macro PR AUC` = manual_macro_avg_pr_auc,
    `Macro ROC AUC` = roc_auc_macro,
    `Macro Sensitivity` = sens_macro,
    `Macro Specificity` = spec_macro,
    Hyperparameters
  ) %>%
  
  # Arrange by the primary metric to rank the models
  arrange(desc(`Macro PR AUC`))

saveRDS(final_leaderboard_data, "Visualisations/Tables/leaderboard_table_raw.rds")

# by-outcome performance champion leaderboard ----
deep_dive_table <- full_results_table %>%
  # First, pivot the five per-class PR AUC columns into a long format
  pivot_longer(
    cols = starts_with("pr_auc_"),
    names_to = "outcome_metric",
    values_to = "pr_auc_value"
  ) %>%
  # Group by the outcome metric
  group_by(outcome_metric) %>%
  # Find the single row with the highest PR AUC within each group
  slice_max(order_by = pr_auc_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # Add the model_name and Hyperparameter summary columns
  mutate(
    model_name = str_extract(source_file, "decision_tree|random_forest|xgboost|svm|elastic_net")
  ) %>%
  rowwise() %>%
  mutate(
    Hyperparameters = case_when(
      model_name == "decision_tree" ~ paste("Complexity Penalty=", format(cost_complexity, scientific = TRUE, digits = 1),
                                            "; Tree Depth=", tree_depth,
                                            "; Min. Node Size=", min_n, sep = ""),
      model_name == "random_forest" ~ paste("Predictors per Split=", mtry, "; Min. Node Size=", min_n, sep = ""),
      model_name == "xgboost" ~ paste("Tree Depth=", tree_depth, "; Min. Node Size=", min_n, 
                                      "; Learning Rate=", round(learn_rate, 3), sep = ""),
      model_name == "svm" ~ paste("Cost=", round(cost, 1), "; Sigma=", round(rbf_sigma, 3), sep = ""),
      model_name == "elastic_net" ~ paste("Penalty=", format(penalty, scientific = TRUE, digits = 2), 
                                          "; Mixture Ratio=", round(mixture, 2), sep = ""),
      TRUE ~ "NA"
    )
  ) %>%
  ungroup() %>%
  
  # Create clean, publication-ready names
  mutate(
    `Outcome Class` = str_remove(outcome_metric, "pr_auc_") %>% str_replace_all("_", " "),
    `Best Performing Model` = str_to_title(str_replace(model_name, "_", " ")),
    `Feature Set` = case_when(
      str_detect(source_file, "all_feature") ~ "All",
      str_detect(source_file, "clin")        ~ "Clinical only",
      str_detect(source_file, "bio|lab")     ~ "Laboratory only"
    )
  ) %>%
  
  # dplyr::select and rename the final columns for the manuscript
  dplyr::select(
    `Outcome Class`,
    `Best Performing Model`,
    `Feature Set`,
    `Per-Class PR AUC` = pr_auc_value,
    `Hyperparameters` # Add the new column
  ) %>%
  
  # Arrange the table in a logical order
  arrange(`Outcome Class`)

saveRDS(deep_dive_table, "Visualisations/Tables/outcome_leaderboard_raw.rds")



# WORK FROM HERE ----