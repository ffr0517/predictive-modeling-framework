# Set up ----
library(dplyr)
library(yardstick)
library(tidyr)
library(ggnewscale)

final_predictions <- readRDS("Trained Models/Two Stage/two_stage_full_run_final_predictions.rds")
operational_logs <- readRDS("Trained Models/Two Stage/two_stage_full_run_operational_logs.rds")
meta_learner_logs <- readRDS("Trained Models/Two Stage/two_stage_meta_learner_logs.rds")
meta_learner_logs <- readRDS("Trained Models/Two Stage/two_stage_meta_learner_logs.rds")
pca_logs <- readRDS("Trained Models/Two Stage/two_stage_pca_logs.rds")

all_preds_df <- bind_rows(final_predictions, .id = "task_id")
logs_df <- bind_rows(operational_logs)
all_coeffs <- bind_rows(meta_learner_logs, .id = "task_id")
all_loadings <- bind_rows(pca_logs, .id = "task_id")

# Functions used ----
ScreePlot <- function(pca_results, n_components = NULL,
                      title = "Scree Plot of Principal Components", 
                      subtitle = "Bars show individual variance; line shows cumulative variance",
                      x_label = "Principal Component",
                      y_label = "Proportion of Variance Explained",
                      variance_threshold = NULL,
                      bar_alpha = 0.7,
                      label_alpha = 1.0) {
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(showtext)
  
  # Enable showtext to render text
  showtext_auto()
  font_add_google("Roboto")
  
  # --- 1. Extract and prepare data from the PCA object ---
  pca_summary_data <- summary(pca_results)$importance
  
  plot_data <- as.data.frame(t(pca_summary_data)) %>%
    dplyr::rename(
      proportion_of_variance = `Proportion of Variance`,
      cumulative_proportion = `Cumulative Proportion`
    ) %>%
    mutate(component = factor(rownames(.), levels = rownames(.)))
  
  # Optionally, limit to the first N components
  if (!is.null(n_components)) {
    plot_data <- plot_data %>% slice(1:n_components)
  }
  
  # --- 2. Build the plot ---
  plot <- ggplot(plot_data, aes(x = component, group = 1)) +
    
    # Bars for individual variance explained, with manual alpha
    geom_col(aes(y = proportion_of_variance), fill = "#0072B2", alpha = bar_alpha) +
    
    # Line and points for cumulative variance
    geom_line(aes(y = cumulative_proportion), color = "#D55E00", linewidth = 1) +
    geom_point(aes(y = cumulative_proportion), color = "#D55E00", size = 3) +
    
    # Add text labels for the bars, with manual alpha
    geom_text(
      aes(y = proportion_of_variance, label = percent(proportion_of_variance, accuracy = 0.1)),
      vjust = -0.4, size = 3.5, family = "Roboto", fontface = "bold", alpha = label_alpha
    ) +
    
    # Set the primary Y-axis to a 0-100% scale
    scale_y_continuous(name = y_label, labels = percent_format(), limits = c(0, 1)) +
    
    # Add an optional threshold line with the specified styling
    {
      if (!is.null(variance_threshold)) {
        geom_hline(
          yintercept = variance_threshold, 
          linetype = "dashed", 
          color = "red",       # <-- Changed to red
          alpha = 0.4          # <-- Alpha set to 0.4
        )
      }
    } +
    
    # --- 3. Apply Styling ---
    labs(title = title, subtitle = subtitle, x = x_label) +
    theme_classic() +
    theme(
      plot.title = element_text(colour = 'black', size = 17.5, face = 'bold', family = 'Roboto'),
      plot.subtitle = element_text(colour = 'grey30', size = 12.5, face = 'italic', family = 'Roboto'),
      axis.title = element_text(family = 'Roboto', size = 14),
      axis.text = element_text(size = 12, family = 'Roboto', colour = 'black'),
      legend.position = 'none'
    )
  
  return(plot)
}


# Creating combined prediction data ----
all_preds_with_logs <- all_preds_df %>%
  left_join(logs_df, by = "task_id")

final_analysis_df <- all_preds_with_logs %>%
  extract(
    task_id,
    into = c("source_file", ".config", "s1_model_index", "repeat_id", "s2_variant"),
    regex = "^(.*\\.csv)_(.*?)_s1_model_(\\d+)_(Repeat\\d+)_(\\w+)$",
    remove = FALSE
  )

master_levels <- tolower(gsub("[- ]", "_", levels(all_preds_df$SFI_5cat)))

analysis_ready_df <- final_analysis_df %>%
  mutate(
    actual = tolower(gsub("[- ]", "_", SFI_5cat)),
    .pred_class = tolower(gsub("[- ]", "_", .pred_class)),
    .pred_class_final = tolower(gsub("[- ]", "_", .pred_class_final)),
    .pred_class_final_safeguard = case_when(
      s1_safeguard_applied == TRUE ~ .pred_class,
      TRUE                       ~ .pred_class_final
    )
  ) %>%
  mutate(
    across(
      c(actual, .pred_class, .pred_class_final, .pred_class_final_safeguard),
      ~ factor(., levels = master_levels)
    )
  )



# Creating meta learner analysis-ready data ----
# The meanings of each principal component
# 1. Calculate the average loading for every term on every component
avg_loadings <- all_loadings %>%
  group_by(component, terms) %>%
  summarise(avg_loading = mean(value), .groups = "drop")

pc_definitions <- avg_loadings %>%
  # 1. Create a column for the base outcome class from the raw 'terms'
  mutate(
    outcome_class = str_replace_all(terms, c("\\.pred_" = "", "\\.\\.\\.\\d+" = ""))
  ) %>%
  # 2. For each component, calculate the total variance to get percentages
  group_by(component) %>%
  mutate(
    pct_contribution = (avg_loading^2 / sum(avg_loading^2)) * 100
  ) %>%
  # 3. Aggregate the total percentage contribution by the base outcome class
  group_by(component, outcome_class) %>%
  summarise(total_pct_contribution = sum(pct_contribution), .groups = "drop") %>%
  # 4. For each component, create the final descriptive string
  group_by(component) %>%
  # Order the classes by importance within each component
  arrange(desc(total_pct_contribution), .by_group = TRUE) %>%
  summarise(
    # Paste the classes and their percentages into a single string
    component_definition = paste(
      paste0(outcome_class, " (", round(total_pct_contribution), "%)"),
      collapse = " | "
    )
  )

avg_pc_influence <- all_coeffs %>%
  # We only care about the PC terms, not the intercept
  #filter(term != "(Intercept)") %>%
  group_by(term, y.level) %>%
  summarise(avg_estimate = mean(estimate), .groups = "drop")

# final summary table ---
final_summary_table <- avg_pc_influence %>%
  # The 'term' column in the coeffs log corresponds to the 'component' column
  left_join(pc_definitions, by = c("term" = "component")) %>%
  # Select and arrange the columns for a clean report
  dplyr::select(component = term, y.level, avg_estimate, component_definition) %>%
  arrange(component, y.level)

saveRDS(final_summary_table, "Visualisations/Tables/stage_two_meta_learner_pca_results.rds")

# FROM BELOW IS THE OLD ANALYSIS STUFF ----
# ----
# Defining metrics of interest ----
multi_class_metrics <- metric_set(
  accuracy,
  bal_accuracy,
  sens,          # Sensitivity (Recall)
  spec,          # Specificity
  ppv,           # Positive Predictive Value (Precision)
  npv,           # Negative Predictive Value
  roc_auc        # Area under the ROC curve
)

# Performance of each of the (24) pipelines ----
# LONG;
performance_summary <- analysis_ready_df %>%
  pivot_longer(
    cols = c(.pred_class_final, .pred_class_final_safeguard),
    names_to = "policy",
    values_to = "prediction"
  ) %>%
  mutate(policy = if_else(policy == ".pred_class_final", "standard", "safeguard")) %>%
  
  group_by(source_file, .config, s2_variant, policy) %>%
  
  multi_class_metrics(
    truth = actual,
    estimate = prediction,
    .pred_Onset_within_24_hours, .pred_Onset_greater_than_24_hours,
    .pred_Probable_Severe, .pred_Probable_Non_Severe, .pred_Non_Severe
  ) %>%
  ungroup()

# WIDE;
per_class_sensitivity_wide <- analysis_ready_df %>%
  pivot_longer(
    cols = c(.pred_class_final, .pred_class_final_safeguard),
    names_to = "policy",
    values_to = "prediction"
  ) %>%
  mutate(policy = if_else(policy == ".pred_class_final", "standard", "safeguard")) %>%
  group_by(source_file, .config, s2_variant, policy) %>%
  summarise(
    # Add event_level = "second" to every call
    sens_non_severe = sens_vec(
      truth = factor(ifelse(actual == "non_severe", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "non_severe", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    sens_probable_non_severe = sens_vec(
      truth = factor(ifelse(actual == "probable_non_severe", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "probable_non_severe", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    sens_probable_severe = sens_vec(
      truth = factor(ifelse(actual == "probable_severe", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "probable_severe", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    sens_onset_greater_than_24_hours = sens_vec(
      truth = factor(ifelse(actual == "onset_greater_than_24_hours", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "onset_greater_than_24_hours", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    sens_onset_within_24_hours = sens_vec(
      truth = factor(ifelse(actual == "onset_within_24_hours", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "onset_within_24_hours", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    .groups = "drop"
  )

# --- Step 2: Calculate the manual macro-averaged sensitivity ---
manual_macro_sens <- per_class_sensitivity_wide %>%
  rowwise() %>%
  mutate(
    sens = mean(c_across(starts_with("sens_")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  dplyr::select(source_file, .config, s2_variant, policy, sens)

# --- Step 3: Calculate all OTHER metrics ---
performance_summary_base <- analysis_ready_df %>%
  pivot_longer(
    cols = c(.pred_class_final, .pred_class_final_safeguard),
    names_to = "policy",
    values_to = "prediction"
  ) %>%
  mutate(policy = if_else(policy == ".pred_class_final", "standard", "safeguard")) %>%
  group_by(source_file, .config, s2_variant, policy) %>%
  summarise(
    accuracy     = accuracy(cur_data(), truth = actual, estimate = prediction)$.estimate,
    bal_accuracy = bal_accuracy(cur_data(), truth = actual, estimate = prediction)$.estimate,
    spec         = spec(cur_data(), truth = actual, estimate = prediction, average = "macro")$.estimate,
    ppv          = ppv(cur_data(), truth = actual, estimate = prediction, average = "macro")$.estimate,
    npv          = npv(cur_data(), truth = actual, estimate = prediction, average = "macro")$.estimate,
    roc_auc      = roc_auc(cur_data(), truth = actual, 
                           .pred_Onset_within_24_hours, .pred_Onset_greater_than_24_hours, 
                           .pred_Probable_Severe, .pred_Probable_Non_Severe, .pred_Non_Severe)$.estimate,
    .groups = "drop"
  )

# --- Step 4: Join the manual sensitivity back to the main summary table ---
performance_summary_wide <- left_join(
  performance_summary_base,
  manual_macro_sens,
  by = c("source_file", ".config", "s2_variant", "policy")
)

per_fold_summary_df <- performance_summary_wide 

# --- This is the correct final aggregation step ---
final_averaged_results <- per_fold_summary_df %>%
  # Group by the columns that define a unique pipeline AND policy
  group_by(source_file, .config, s2_variant, policy) %>%
  
  # Calculate the mean of the per-fold metrics across all folds
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop" 
  )

# Stage one predictions most often handed off to stage two ----
# --- Step 1: Calculate handoff counts and proportions for EACH FOLD ---
s1_handoff_triggers_per_fold <- analysis_ready_df %>%
  filter(s2_invoked == TRUE) %>%
  # Group by all identifiers for a unique run, including the fold
  group_by(source_file, .config, s1_model_index, repeat_id, s2_variant, .pred_class) %>%
  summarise(handoff_count = n(), .groups = "drop") %>%
  # Now calculate proportions within each fold
  group_by(source_file, .config, s1_model_index, repeat_id, s2_variant) %>%
  mutate(proportion = handoff_count / sum(handoff_count)) %>%
  ungroup()

# --- Step 2: Average the per-fold results to get the final summary ---
s1_handoff_triggers_avg <- s1_handoff_triggers_per_fold %>%
  # Group by the pipeline identifiers and the S1 prediction class
  group_by(source_file, .config, s2_variant, .pred_class) %>%
  # Calculate the average count and proportion across all folds
  summarise(
    avg_handoff_count = mean(handoff_count, na.rm = TRUE),
    avg_proportion = mean(proportion, na.rm = TRUE),
    .groups = "drop"
  )

# --- Step 3: Pivot the averaged results into the final wide format ---
handoff_counts_wide <- s1_handoff_triggers_avg %>%
  pivot_wider(
    id_cols = c(source_file, .config, s2_variant),
    names_from = .pred_class,
    values_from = avg_handoff_count,
    names_prefix = "avg_count_handoff_",
    values_fill = 0
  )

handoff_props_wide <- s1_handoff_triggers_avg %>%
  pivot_wider(
    id_cols = c(source_file, .config, s2_variant),
    names_from = .pred_class,
    values_from = avg_proportion,
    names_prefix = "avg_prop_handoff_",
    values_fill = 0
  )

# Join the two wide tables to create the final summary
s1_handoff_triggers_wide <- left_join(
  handoff_counts_wide,
  handoff_props_wide,
  by = c("source_file", ".config", "s2_variant")
)

# Accuracy of each stage 1 prediction class ----
s1_class_accuracy_per_fold <- analysis_ready_df %>%
  group_by(source_file, .config, s1_model_index, repeat_id, s2_variant, .pred_class) %>%
  summarise(
    n_observations = n(),
    s1_accuracy = mean(actual == .pred_class, na.rm = TRUE),
    .groups = "drop"
  )

s1_class_accuracy_avg <- s1_class_accuracy_per_fold %>%
  group_by(source_file, .config, s2_variant, .pred_class) %>%
  summarise(
    avg_n_observations = mean(n_observations, na.rm = TRUE),
    avg_s1_accuracy = mean(s1_accuracy, na.rm = TRUE),
    .groups = "drop"
  )

s1_class_accuracy_wide <- s1_class_accuracy_avg %>%
  pivot_wider(
    id_cols = c(source_file, .config, s2_variant),
    names_from = .pred_class,
    values_from = c(avg_n_observations, avg_s1_accuracy),
    names_glue = "{.value}_{.pred_class}",
    values_fill = list(avg_n_observations = 0)
  )

# Average impact of the safeguard policy ----
# Pivot the summary table to have standard and safeguard metrics side-by-side
per_fold_performance <- analysis_ready_df %>%
  pivot_longer(
    cols = c(.pred_class_final, .pred_class_final_safeguard),
    names_to = "policy",
    values_to = "prediction"
  ) %>%
  mutate(policy = if_else(policy == ".pred_class_final", "standard", "safeguard")) %>%
  group_by(source_file, .config, s1_model_index, repeat_id, s2_variant, policy) %>%
  summarise(
    bal_accuracy = bal_accuracy_vec(truth = actual, estimate = prediction),
    sens = sens_vec(truth = actual, estimate = prediction),
    spec = spec_vec(truth = actual, estimate = prediction),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(bal_accuracy, sens, spec),
    names_to = ".metric",
    values_to = ".estimate"
  )

per_fold_safeguard_impact <- per_fold_performance %>%
  pivot_wider(
    names_from = policy,
    values_from = .estimate
  ) %>%
  mutate(
    delta = safeguard - standard
  )

averaged_safeguard_impact <- per_fold_safeguard_impact %>%
  group_by(source_file, .config, s2_variant, .metric) %>%
  summarise(
    across(c(standard, safeguard, delta), mean, na.rm = TRUE),
    .groups = "drop"
  )

safeguard_impact_wide <- averaged_safeguard_impact %>%
  pivot_wider(
    id_cols = c(source_file, .config, s2_variant),
    names_from = .metric,
    values_from = c(standard, safeguard, delta),
    names_glue = "{.value}_{.metric}"
  )

# How often S2 was unhelpful versus helpful ----
per_fold_net_benefit <- analysis_ready_df %>%
  filter(s2_invoked == TRUE) %>%
  group_by(source_file, .config, s1_model_index, repeat_id, s2_variant) %>%
  summarise(
    n_helpful_s2 = sum(actual != .pred_class & actual == .pred_class_final, na.rm = TRUE),
    n_unhelpful_s2 = sum(actual == .pred_class & actual != .pred_class_final, na.rm = TRUE),
    .groups = "drop"
  )

s2_net_benefit_avg <- per_fold_net_benefit %>%
  group_by(source_file, .config, s2_variant) %>%
  summarise(
    n_helpful_s2 = mean(n_helpful_s2, na.rm = TRUE),
    n_unhelpful_s2 = mean(n_unhelpful_s2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(s2_net_benefit = n_helpful_s2 - n_unhelpful_s2)

# Joining all wide data ----
join_keys <- c("source_file", ".config", "s2_variant")

# Now, start with your main performance table and join everything to it
combined_analysis_df <- safeguard_impact_wide %>%
  left_join(s2_net_benefit_avg, by = join_keys) %>%
  left_join(s1_handoff_triggers_wide, by = join_keys) %>%
  left_join(s1_class_accuracy_wide, by = join_keys) 

# Adding macro pr auc ----
per_fold_pr_auc <- analysis_ready_df %>%
  group_by(source_file, .config, s1_model_index, repeat_id, s2_variant) %>%
  summarise(
    pr_auc_non_severe = pr_auc_vec(
      truth = factor(ifelse(actual == "non_severe", "event", "other"), levels = c("other", "event")),
      estimate = .data[[".pred_Non_Severe"]],
      event_level = "second"
    ),
    pr_auc_probable_non_severe = pr_auc_vec(
      truth = factor(ifelse(actual == "probable_non_severe", "event", "other"), levels = c("other", "event")),
      estimate = .data[[".pred_Probable_Non_Severe"]],
      event_level = "second"
    ),
    pr_auc_probable_severe = pr_auc_vec(
      truth = factor(ifelse(actual == "probable_severe", "event", "other"), levels = c("other", "event")),
      estimate = .data[[".pred_Probable_Severe"]],
      event_level = "second"
    ),
    pr_auc_onset_greater_than_24_hours = pr_auc_vec(
      truth = factor(ifelse(actual == "onset_greater_than_24_hours", "event", "other"), levels = c("other", "event")),
      estimate = .data[[".pred_Onset_greater_than_24_hours"]],
      event_level = "second"
    ),
    pr_auc_onset_within_24_hours = pr_auc_vec(
      truth = factor(ifelse(actual == "onset_within_24_hours", "event", "other"), levels = c("other", "event")),
      estimate = .data[[".pred_Onset_within_24_hours"]],
      event_level = "second"
    ),
    .groups = "drop"
  )

avg_per_class_pr_auc <- per_fold_pr_auc %>%
  group_by(source_file, .config, s2_variant) %>%
  summarise(
    across(starts_with("pr_auc_"), mean, na.rm = TRUE),
    .groups = "drop"
  )

pipeline_pr_auc <- avg_per_class_pr_auc %>%
  rowwise() %>%
  mutate(
    pr_auc_macro = mean(c_across(starts_with("pr_auc_")), na.rm = TRUE)
  ) %>%
  ungroup()

combined_analysis_df <- combined_analysis_df %>%
  dplyr::select(-any_of("pr_auc_macro")) %>%
  left_join(
    pipeline_pr_auc,
    by = c("source_file", ".config", "s2_variant")
  )

# Creating dataframe for updated pareto-front plot ----
per_fold_sensitivity <- analysis_ready_df %>%
  pivot_longer(
    cols = c(.pred_class_final, .pred_class_final_safeguard),
    names_to = "policy",
    values_to = "prediction"
  ) %>%
  mutate(policy = if_else(policy == ".pred_class_final", "standard", "safeguard")) %>%
  group_by(source_file, .config, s1_model_index, repeat_id, s2_variant, policy) %>%
  summarise(
    sens_non_severe = sens_vec(
      truth = factor(ifelse(actual == "non_severe", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "non_severe", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    sens_probable_non_severe = sens_vec(
      truth = factor(ifelse(actual == "probable_non_severe", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "probable_non_severe", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    sens_probable_severe = sens_vec(
      truth = factor(ifelse(actual == "probable_severe", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "probable_severe", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    sens_onset_greater_than_24_hours = sens_vec(
      truth = factor(ifelse(actual == "onset_greater_than_24_hours", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "onset_greater_than_24_hours", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    sens_onset_within_24_hours = sens_vec(
      truth = factor(ifelse(actual == "onset_within_24_hours", "event", "other"), levels = c("other", "event")),
      estimate = factor(ifelse(prediction == "onset_within_24_hours", "event", "other"), levels = c("other", "event")),
      event_level = "second"
    ),
    .groups = "drop"
  )

per_class_sensitivity_wide <- per_fold_sensitivity %>%
  group_by(source_file, .config, s2_variant, policy) %>%
  summarise(
    across(starts_with("sens_"), mean, na.rm = TRUE),
    .groups = "drop"
  )

# --- Step 2: Calculate the average sensitivity across the severe classes ---
severe_sens_cols <- c("sens_onset_within_24_hours", "sens_onset_greater_than_24_hours")

avg_severe_sens <- per_class_sensitivity_wide %>%
  # Use rowwise() to calculate the mean across the specified columns for each row
  rowwise() %>%
  mutate(
    avg_severe_sensitivity = mean(c_across(all_of(severe_sens_cols)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # dplyr::select only the columns we need for plotting
  dplyr::select(source_file, .config, s2_variant, policy, avg_severe_sensitivity)

avg_severe_sens_wide <- avg_severe_sens %>%
  pivot_wider(
    id_cols = c(source_file, .config, s2_variant),
    names_from = policy,
    values_from = avg_severe_sensitivity,
    names_glue = "avg_severe_sens_{policy}"
  )

plotting_df <- combined_analysis_df %>%
  left_join(avg_severe_sens_wide, by = c("source_file", ".config", "s2_variant")) %>%
  
  mutate(
    model_name = case_when(
      .config == "Preprocessor1_Model45" ~ "XGBoost - Maximum Efficiency Profile",
      .config == "Preprocessor1_Model18" ~ "XGBoost - Efficiency-Balanced Profile",
      .config == "Preprocessor1_Model37" ~ "SVM - Safety-Balanced Profile",
      .config == "Preprocessor1_Model53" ~ "SVM - Maximum Safety Profile",
      .config == "Preprocessor1_Model64" ~ "XGBoost - Balanced Profile",
      .config == "Preprocessor1_Model03" ~ "Random Forest - Onset < 24h Specialist",
      # A fallback to catch any other .config values, just in case
      TRUE                               ~ .config 
    )
  )

plotting_df <- plotting_df %>%
  mutate(
    model_type = str_extract(source_file, "xgboost|svm|elastic_net|random_forest")
  )

# Two-stage tradeoff plot ----
df_ranked_plot <- readRDS("Visualisations/df_ranked_plot_data.rds")

df_ranked_plot_std <- df_ranked_plot %>%
  rename(
    pr_auc_non_severe = `pr_auc_Non-Severe`, 
    avg_severe_sensitivity = avg_severe_sensitivity
  )


tradeoff_plot_with_baseline <- ggplot(plotting_df) +
  # Arrows: standard → safeguard
  geom_segment(
    aes(x = pr_auc_non_severe, y = avg_severe_sens_standard,
        xend = pr_auc_non_severe, yend = avg_severe_sens_safeguard,
        colour = model_name),
    arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7, show.legend = FALSE
  ) +
  # Standard points (hollow)
  geom_point(
    aes(x = pr_auc_non_severe, y = avg_severe_sens_standard,
        colour = model_name, shape = s2_variant),
    size = 4, stroke = 1, fill = "white"
  ) +
  # Safeguard points (solid)
  geom_point(
    aes(x = pr_auc_non_severe, y = avg_severe_sens_safeguard,
        colour = model_name, fill = model_name, shape = s2_variant),
    size = 4, stroke = 1
  ) +
  
  # --- (Rest of your original code) ---
  scale_colour_manual(
    name = "Stage One Profile",
    values = c("XGBoost - Maximum Efficiency Profile" = "#0072B2",
               "XGBoost - Efficiency-Balanced Profile" = "#56B4E9",
               "XGBoost - Balanced Profile" = "#D55E00",
               "SVM - Safety-Balanced Profile" = "#E69F00",
               "SVM - Maximum Safety Profile" = "#009E73",
               "Random Forest - Onset < 24h Specialist" = "#CC79A7")
  ) +
  scale_fill_manual(
    values = c("XGBoost - Maximum Efficiency Profile" = "#0072B2",
               "XGBoost - Efficiency-Balanced Profile" = "#56B4E9",
               "XGBoost - Balanced Profile" = "#D55E00",
               "SVM - Safety-Balanced Profile" = "#E69F00",
               "SVM - Maximum Safety Profile" = "#009E73",
               "Random Forest - Onset < 24h Specialist" = "#CC79A7"),
    guide = "none" 
  ) +
  scale_shape_manual(
    name = "Stage Two Variant",
    values = c(champion_xgb = 21, stacker = 22)
  ) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(fill = "black", size = 4)),
    colour = guide_legend(order = 2)
  ) +
  labs(
    title = "Trade-off Between Safety and Efficiency for Two-Stage Models",
    subtitle = "Solid Shapes Showcase Estimates with Safeguarding Policy Applied",
    x = "Precision-Recall Area Under the Curve for 'Non-Severe' Cases",
    y = "Mean Sensitivity for Severe Cases"
  ) +
  theme_classic(base_family = "Roboto") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic", color = "grey30"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.box = "vertical",
    # --- New lines to make the legend smaller ---
    legend.key.size = unit(0.5, "cm"),      # Adjusts the size of the legend keys
    legend.title = element_text(size = 10), # Adjusts the legend title font size
    legend.text = element_text(size = 8)    # Adjusts the legend item font size
  )


ggsave("WITHLEGENDclinical_trade_off_plot_two_stage.svg", plot = tradeoff_plot_with_baseline, width = 10, height = 5.5, units = "in", dpi = 300)


# Net benefit plot ----
plot_data_net_benefit <- s2_net_benefit_avg %>%
  left_join(
    # dplyr::select just the naming/typing columns from plotting_df to join
    plotting_df %>% distinct(source_file, .config, s2_variant, model_name, model_type),
    by = c("source_file", ".config", "s2_variant")
  )

base_colors <- c(
  "XGBoost - Maximum Efficiency Profile" = "#0072B2",
  "XGBoost - Efficiency-Balanced Profile" = "#56B4E9",
  "XGBoost - Balanced Profile" = "#D55E00",
  "SVM - Safety-Balanced Profile" = "#E69F00",
  "SVM - Maximum Safety Profile" = "#009E73",
  "Random Forest - Onset < 24h Specialist" = "#CC79A7"
)

darker_colors <- c(
  "XGBoost - Maximum Efficiency Profile" = "#004C70", # Darker blue
  "XGBoost - Efficiency-Balanced Profile" = "#347C9D", # Darker light blue
  "XGBoost - Balanced Profile" = "#8F3F00", # Darker orange
  "SVM - Safety-Balanced Profile" = "#A16F00", # Darker yellow-orange
  "SVM - Maximum Safety Profile" = "#00694C", # Darker green
  "Random Forest - Onset < 24h Specialist" = "#8C5373"  # Darker pink
)

final_color_palette <- c(
  setNames(base_colors, paste(names(base_colors), "S1 Champion")),
  setNames(darker_colors, paste(names(darker_colors), "Meta-Learner"))
)


net_benefit_plot_final <- plot_data_net_benefit %>%
  mutate(
    s2_variant_named = case_when(
      s2_variant == "champion_xgb" ~ "S1 Champion",
      s2_variant == "stacker"      ~ "Meta-Learner",
      TRUE                         ~ as.character(s2_variant)
    ),
    plot_label = paste0(model_name, " (", s2_variant_named, ")"),
    model_combo = paste(model_name, s2_variant_named)
  ) %>%
  ggplot(aes(x = reorder(plot_label, s2_net_benefit), y = s2_net_benefit, fill = model_combo)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  
  scale_fill_manual(values = final_color_palette) +
  
  labs(
    title = "Effectiveness of Stage 2 Cascade Models",
    subtitle = "Net benefit is the number of corrected Stage 1 errors minus newly introduced errors",
    x = "Cascade Pipeline Profile",
    y = "Net Corrected Predictions (S2 Net Benefit)"
  ) +
  theme_minimal(base_family = "Roboto") +
  theme(
    plot.title = element_text(size = 18, face = 'bold'),
    plot.subtitle = element_text(size = 12, face = 'italic', color = 'grey30'),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 9),
    legend.position = "none" 
  )

ggsave("net_benefit_plot_two_stage.svg", plot = net_benefit_plot_final, width = 10, height = 5.5, units = "in", dpi = 300)

# Final tables creation ----
severe_obs_cols <- c("avg_n_observations_onset_within_24_hours", "avg_n_observations_onset_greater_than_24_hours")
severe_handoff_cols <- c("avg_count_handoff_onset_within_24_hours", "avg_count_handoff_onset_greater_than_24_hours")

# Add the two new handoff rate columns
plotting_df_with_handoffs <- plotting_df %>%
  rowwise() %>%
  mutate(
    # Total Handoff Rate: Sum of all handoffs / sum of all observations
    `% Handoff Rate` = (sum(c_across(starts_with("avg_count_handoff_"))) / 
                          sum(c_across(starts_with("avg_n_observations_")))) * 100,
    
    # Severe Case Handoff Rate: Sum of severe handoffs / sum of severe observations
    `% Severe Case Handoff Rate` = (sum(c_across(all_of(severe_handoff_cols))) / 
                                      sum(c_across(all_of(severe_obs_cols)))) * 100
  ) %>%
  ungroup()


standard_summary <- plotting_df_with_handoffs %>%
  dplyr::select(
    model_name, s2_variant,
    bal_accuracy = standard_bal_accuracy,
    sens = standard_sens,
    spec = standard_spec,
    pr_auc_macro,
    avg_severe_sensitivity = avg_severe_sens_standard,
    pr_auc_non_severe,
    `% Handoff Rate`, 
    `% Severe Case Handoff Rate`, 
    s2_net_benefit
  ) %>%
  mutate(Policy = "Standard")

safeguard_summary <- plotting_df_with_handoffs %>%
  dplyr::select(
    model_name, s2_variant,
    bal_accuracy = safeguard_bal_accuracy,
    sens = safeguard_sens,
    spec = safeguard_spec,
    pr_auc_macro,
    avg_severe_sensitivity = avg_severe_sens_safeguard,
    pr_auc_non_severe,
    `% Handoff Rate`, 
    `% Severe Case Handoff Rate`, 
    s2_net_benefit
  ) %>%
  mutate(Policy = "Safeguard")

final_results_table <- bind_rows(standard_summary, safeguard_summary) %>%
  rename(
    `Stage One Model` = model_name,
    `Stage Two Model` = s2_variant,
    `Balanced Accuracy` = bal_accuracy,
    `Macro Sensitivity` = sens,
    `Macro Specificity` = spec,
    `Macro PR AUC` = pr_auc_macro,
    `Avg. Severe Sensitivity` = avg_severe_sensitivity,
    `PR AUC (Non-Severe)` = pr_auc_non_severe,
    `S2 Net Corrections` = s2_net_benefit
  ) %>%
  arrange(desc(`Balanced Accuracy`), desc(`Avg. Severe Sensitivity`)) %>%
  mutate(
    `Strategy ID` = row_number()
  ) %>%
  dplyr::select(
    `Strategy ID`,
    `Stage One Model`,
    `Stage Two Model`,
    Policy,
    `Balanced Accuracy`,
    `Avg. Severe Sensitivity`,
    `PR AUC (Non-Severe)`,
    `Macro PR AUC`,
    `Macro Sensitivity`,
    `Macro Specificity`,
    `% Handoff Rate`,
    `% Severe Case Handoff Rate`,
    `S2 Net Corrections`
  ) %>% 
  dplyr::mutate(
    `Stage Two Model` = dplyr::recode(
      `Stage Two Model`,
      "stacker" = "Meta-Leaner",
      "champion_xgb" = "S1 Champion"
    )
  ) 

# View the final table
print(final_results_table)

saveRDS(final_results_table, "Visualisations/Tables/stage_two_final_table.rds")

# PCA ----
# Create the data frame for PCA
pca_input_df <- final_results_table %>%
  dplyr::select(-`Strategy ID`) %>%
  # Select only numeric columns
  dplyr::select(where(is.numeric)) %>%
  # Replace any remaining NAs with 0
  mutate(across(everything(), ~ifelse(is.na(.), 0, .)))

# Store identifiers separately for plotting
model_identifiers <- final_results_table %>%
  dplyr::select(`Strategy ID`, `Stage One Model`, `Stage Two Model`, `Policy`)


# Running PCA
pca_results <- prcomp(pca_input_df, scale. = TRUE, center = TRUE)

summary(pca_results)

pca_results$rotation

# scores
pca_scores <- as.data.frame(pca_results$x) %>%
  bind_cols(model_identifiers)

# plotting
scree_plot <- ScreePlot(
  pca_results = pca_results, 
  variance_threshold = 0.90,
  bar_alpha = 0.2,    
  label_alpha = 0.2,
  title = "Variance Explained in Performance Metrics Across Two-Stage Models"
)

ggsave("scree_plot.svg", plot = scree_plot, width = 10, height = 5.5, units = "in", dpi = 300)

loadings_df <- as.data.frame(pca_results$rotation) %>%
  dplyr::select(PC1, PC2, PC3) %>%
  arrange(-PC1)

pca_scores_with_labels <- pca_scores %>%
  mutate(
    pipeline_label = paste0(`Stage One Model`, " (", `Stage Two Model`, ")", " + ", `Policy`)
  )

pca_subset <- pca_scores_with_labels %>%
  dplyr::select(`Stage One Model`, `Stage Two Model`, `Policy`, PC1, PC2, PC3) %>%
  dplyr::mutate(
    `Stage Two Model` = dplyr::recode(
      `Stage Two Model`,
      "stacker" = "Meta-Leaner",
      "champion_xgb" = "S1 Champion"
    )
  ) %>%
  dplyr::arrange(PC1)

saveRDS(loadings_df, "Visualisations/Tables/stage_two_loadings_df.rds")
saveRDS(pca_subset, "Visualisations/Tables/stage_two_pca_scores_df.rds")


