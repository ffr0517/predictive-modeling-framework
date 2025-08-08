# Set up ----
library(tidyverse)
library(vip) 
library(fastshap)
library(broom) 


# Full model
full_cascade_model <- readRDS("Fully Trained Cascade Model/full_cascade_model.rds")

# Stage one model feature importance ----
# Stage 1 model fit
s1_model_fit <- full_cascade_model$stage_1_fit

# top 20 features for the initial screening model.
vip(s1_model_fit, num_features = 20, geom = "col") +
  labs(
    title = "Feature Importance: Stage 1 'Screener' Model",
    subtitle = "Features used for the initial classification of all cases"
  ) +
  theme_minimal()

# variable importance scores as a data frame
s1_vi_df <- vi(s1_model_fit)

s1_top_feat_df <- s1_vi_df %>%
  arrange(desc(Importance)) 

# View the resulting table
s1_top_feat_df


# Stage two specialist and meta learner feature importance ----
# list of specialist models
s2_specialists <- full_cascade_model$s2_specialist_fits

# Probable severe/probable non-severe specialist ---
probable_specialist <- s2_specialists$Probable_Severe

vip(probable_specialist, num_features = 20, geom = "col") +
  labs(
    title = "Feature Importance: Stage 2 'Probable' Classification(s) Specialist Model",
    subtitle = "Features used by the max performing* model for probable non-severe and probable severe patients handed off from Stage 1",
    caption = "Max performance defined by given outcome PR AUC performance"
  ) +
  theme_minimal()

# variable importance scores as a data frame
prob_specialist_vi_df <- vi(probable_specialist)

prob_specialist_top_feat_df <- prob_specialist_vi_df %>%
  arrange(desc(Importance)) 

# Onset greater than 24 hours specialist ---
greater_than_24_hours_specialist <- s2_specialists$Onset_greater_than_24_hours

vip(greater_than_24_hours_specialist, num_features = 20, geom = "col") +
  labs(
    title = "Feature Importance: Stage 2 'Onset greater than 24 hours' Classification Specialist Model",
    subtitle = "Features used by the max performing* model for severe outcome patients with an onset greater than 24 hours handed off from Stage 1",
    caption = "Max performance defined by given outcome PR AUC performance"
  ) +
  theme_minimal()

# variable importance scores as a data frame
gr_thn_24h_specialist_vi_df <- vi(greater_than_24_hours_specialist)

gr_thn_24h_specialist_top_feat_df <- gr_thn_24h_specialist_vi_df %>%
  arrange(desc(Importance)) 

# Non severe specialist ---
non_severe_specialist <- s2_specialists$Non_Severe

vip(non_severe_specialist, num_features = 20, geom = "col") +
  labs(
    title = "Feature Importance: Stage 2 'non-severe' Classification Specialist Model",
    subtitle = "Features used by the max performing* model for non-severe outcome patients handed off from Stage 1",
    caption = "Max performance defined by given outcome PR AUC performance"
  ) +
  theme_minimal()

# variable importance scores as a data frame
non_severe_specialist_vi_df <- vi(non_severe_specialist)

non_severe_specialist_top_feat_df <- non_severe_specialist_vi_df %>%
  arrange(desc(Importance)) 

# Onset within 24 hours specialist ---
# NOTE: Why post-hoc permutation approach is required for this model:
# The 'ranger' (Random Forest) engine does not calculate feature importance by default. This must be requested with an `importance` flag during the initial model training/fitting stage. The 'xgboost' engine, in contrast, DOES calculate and store feature importance metrics (based on 'Gain') by default during its training process. Since this final Random Forest specialist model was saved without the pre-computed importance scores, we must calculate them now using this post-hoc permutation wrapper. This was not necessary for the XGBoost models.

# predictor sets 
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

# data used by the meta models
train_sub <- readRDS("Final Cascade Meta Learner Datasets/train_sub.rds")

# Onset within 24 hours specialist
onset_within_24_hours_specialist <- full_cascade_model$s2_specialist_fits$Onset_within_24_hours

# original recipe directly from the workflow object
onset_24h_specialist_recipe <- pull_workflow_preprocessor(onset_within_24_hours_specialist)

# extracted recipe using same training data
onset_24h_prepped_recipe <- prep(onset_24h_specialist_recipe, training = train_sub)

# prediction wrapper
pfun_simple <- function(object, newdata) {
  predict(object, new_data = newdata, type = "prob")
}

# baseline performance

# predictions on the original, unshuffled data
baseline_preds <- pfun_simple(onset_within_24_hours_specialist, train_sub)

# combining with truth to calculate the baseline metric
baseline_data <- bind_cols(
  truth = train_sub$SFI_5cat,
  baseline_preds
)

# baseline log-loss
baseline_metric <- mn_log_loss(baseline_data, truth, starts_with(".pred"))

# manual permutation loop
# excluding the outcome variable 'SFI_5cat'
predictor_vars <- setdiff(names(train_sub), "SFI_5cat")

importance_results <- list()


# Loop through each predictor variable
set.seed(55378008)

for (var in predictor_vars) {
  
  # Create a temporary shuffled dataset
  shuffled_data <- train_sub
  
  # Shuffle only the current variable's column
  shuffled_data[[var]] <- sample(shuffled_data[[var]])
  
  # Get predictions on the shuffled data
  shuffled_preds <- pfun_simple(onset_within_24_hours_specialist, shuffled_data)
  
  # Combine with truth
  shuffled_metric_data <- bind_cols(
    truth = shuffled_data$SFI_5cat,
    shuffled_preds
  )
  
  # Calculate the new metric
  shuffled_metric <- mn_log_loss(shuffled_metric_data, truth, starts_with(".pred"))
  
  # Store the results
  importance_results[[var]] <- shuffled_metric$.estimate
  
  cat(".")
}

# tidy data frame
onset_within_24h_specialist_top_feat_df <- enframe(unlist(importance_results), name = "Variable", value = "Shuffled_LogLoss") %>%
  mutate(Importance = Shuffled_LogLoss - baseline_metric$.estimate) %>%
  select(Variable, Importance) %>%
  arrange(desc(Importance))

# Renaming features for clearer interpretation ----
library(knitr) 

# categories for each variable
demographic_predictors <- c("site", "age.months", "sex", "travel.time.bin")
anthropometric_predictors <- c("wfaz", "waste", "stunt")
clinical_symptom_predictors <- c("bgcombyn", "adm.recent", "cidysymp", "prior.care", "urti", "lrti", "diarrhoeal", "neuro", "auf", "syndrome.resp", "syndrome.nonresp", "pneumo", "sev.pneumo", "ensapro", "vomit.all", "seiz", "pfacleth", "not.alert")
clinical_sign_predictors <- c("danger.sign", "hr.all", "rr.all", "oxy.ra", "envhtemp", "crt.long", "LqSOFA", "parenteral_screen", "SIRS")
lab_predictors <- c("ANG1", "ANG2", "CHI3L", "CRP", "CXCl10", "IL1ra", "IL6", "IL8", "IL10", "PROC", "TNFR1", "STREM1", "VEGFR1", "supar", "lblac", "lbglu", "enescbchb1")

# single lookup data frame for variable categories
variable_categories <- tibble(Variable = c(demographic_predictors, anthropometric_predictors, clinical_symptom_predictors, clinical_sign_predictors, lab_predictors)) %>%
  mutate(Category = case_when(
    Variable %in% demographic_predictors ~ "Demographic",
    Variable %in% anthropometric_predictors ~ "Anthropometric",
    Variable %in% clinical_symptom_predictors ~ "Clinical Symptom",
    Variable %in% clinical_sign_predictors ~ "Clinical Sign",
    Variable %in% lab_predictors ~ "Laboratory",
    TRUE ~ "Unknown"
  ))

s1_importance <- s1_top_feat_df

probable_importance <- prob_specialist_top_feat_df

gt24h_importance <- gr_thn_24h_specialist_top_feat_df

nonsevere_importance <- non_severe_specialist_top_feat_df

onset_within_24h_importance <- onset_within_24h_specialist_top_feat_df

unique_variable_names <- c(
  "adm.recent", "adm.recent_X1", "age.months", "ANG1", "ANG2", "auf",
  "auf_X1", "bgcombyn", "bgcombyn_X1", "CHI3L", "cidysymp", "CRP",
  "crt.long", "crt.long_X1", "CXCl10", "danger.sign", "danger.sign_X1",
  "diarrhoeal", "diarrhoeal_X1", "enescbchb1", "ensapro", "ensapro_X1",
  "envhtemp", "hr.all", "IL10", "IL1ra", "IL6", "IL8", "infection",
  "ipdopd", "label", "lbglu", "lblac", "LqSOFA", "lrti", "lrti_X1",
  "na_ind_adm.recent", "na_ind_ANG1", "na_ind_ANG2", "na_ind_CHI3L",
  "na_ind_CRP", "na_ind_CXCl10", "na_ind_danger.sign",
  "na_ind_enescbchb1", "na_ind_IL10", "na_ind_IL1ra", "na_ind_IL6",
  "na_ind_IL8", "na_ind_lbglu", "na_ind_lblac", "na_ind_LqSOFA",
  "na_ind_oxy.ra", "na_ind_pfacleth", "na_ind_PROC", "na_ind_rr.all",
  "na_ind_seiz", "na_ind_SIRS", "na_ind_STREM1", "na_ind_stunt",
  "na_ind_supar", "na_ind_TNFR1", "na_ind_VEGFR1", "na_ind_vomit.all",
  "na_ind_waste", "neuro", "neuro_X1", "not.alert", "not.alert_X1",
  "outcome.binary", "oxy.ra", "parenteral_screen", "parenteral_screen_X1",
  "pfacleth", "pfacleth_X1", "pneumo", "pneumo_X1", "prior.care",
  "prior.care_X1", "PROC", "rr.all", "seiz", "seiz_X1", "sev.pneumo",
  "sev.pneumo_X1", "sex", "sex_X1", "SIRS", "site", "site_id003",
  "site_kh005", "site_la004", "site_la011", "site_vn009", "site_vn010",
  "STREM1", "stunt", "stunt_X1", "supar", "syndrome.nonresp",
  "syndrome.nonresp_X1", "syndrome.resp", "syndrome.resp_X1", "TNFR1",
  "travel.time.bin", "travel.time.bin_X1", "urti", "urti_X1", "VEGFR1",
  "vomit.all", "vomit.all_X1", "waste", "waste_X1", "weight", "weight_bin",
  "wfaz"
)


user_provided_key <- tribble(
  ~base_name,                ~clean_base_name,
  # --- Definitions from your image ---
  "adm.recent",              "Overnight Hospitalisation (last 6 months)",
  "age.months",              "Age (Months)",
  "ANG1",                    "Angiopoietin-1",
  "ANG2",                    "Angiopoietin-2",
  "auf",                     "Undifferentiated Febrile Illness",
  "bgcombyn",                "Comorbidity",
  "CHI3L",                   "CHI3L1 (ng/ml)",
  "cidysymp",                "Duration of Illness (Days)",
  "CRP",                     "C-Reactive Protein",
  "crt.long",                "Capillary Refill Time > 2s",
  "CXCl10",                  "CXCL10",
  "danger.sign",             "WHO Danger Sign",
  "diarrhoeal",              "Diarrhoeal Syndrome",
  "enescbchb1",              "Haemoglobin (mg/dL)",
  "ensapro",                 "Prostration",
  "envhtemp",                "Axillary Temperature (celcius)",
  "hr.all",                  "Heart Rate (bpm)",
  "IL10",                    "Interleukin-10",
  "IL1ra",                   "Interleukin-1 Receptor Antagonist",
  "IL6",                     "Interleukin-6 (pg/ml)",
  "IL8",                     "Interleukin-8",
  "infection",               "Microbiologically Confirmed Infection",
  "ipdopd",                  "Admission Status",
  "lbglu",                   "Glucose (mmol/L)",
  "lblac",                   "Lactate (mmol/L)",
  "LqSOFA",                  "LqSOFA Score",
  "lrti",                    "Lower Respiratory Tract Infection",
  "neuro",                   "Neurological Syndrome",
  "not.alert",               "Not Alert (AVPU < A)",
  "oxy.ra",                  "Oxygen Saturation (%)",
  "parenteral_screen",       "Parenteral Treatment Before Enrolment",
  "pfacleth",                "Lethargy",
  "pneumo",                  "WHO Pneumonia",
  "prior.care",              "Sought Care Prior to Presentation",
  "PROC",                    "Procalcitonin",
  "rr.all",                  "Respiratory Rate (breaths/min)",
  "seiz",                    "Convulsions",
  "sev.pneumo",              "WHO Severe Pneumonia",
  "sex",                     "Sex",
  "SIRS",                    "SIRS Score",
  "site",                    "Clinical Site (Overall)",
  "STREM1",                  "STREM1 (pg/ml)",
  "stunt",                   "Stunting",
  "supar",                   "suPAR (ng/ml)",
  "syndrome.nonresp",        "Non-Respiratory Syndrome",
  "syndrome.resp",           "Respiratory Syndrome",
  "TNFR1",                   "TNFR1 (pg/ml)",
  "travel.time.bin",         "Travel Time to Study Site",
  "urti",                    "Upper Respiratory Tract Infection",
  "VEGFR1",                  "VEGFR1 (pg/ml)",
  "vomit.all",               "Intractable Vomiting",
  "waste",                   "Wasting",
  "wfaz",                    "Weight for Age Z-Score",
  "label",                   "Label",
  "outcome.binary",          "Binary Outcome",
  "weight",                  "Weight (five class)",
  "weight_bin",              "Weight (binary)",
  "site_id003",              "Site: Rumah Sakit Umum Daerah Wates, Indonesia", 
  "site_kh005",              "Site: Angkor Hospital, Cambodia",     
  "site_la004",              "Site: Salavan Provincial Hospital, Laos",
  "site_la011",              "Site: Savannakhet Provincial Hospital, Laos",
  "site_vn009",              "Site: Dong Nai Childrens Hospital, Viet Nam",   
  "site_vn010",              "Site: National Childrens Hospital, Viet Nam"    
)



# adding suffixes for missing data and dummy variables
variable_lookup_table <- tibble(original_name = unique_variable_names) %>%
  mutate(
    base_name = if_else(
      str_starts(original_name, "na_ind_"),
      str_remove(original_name, "^na_ind_"),
      original_name %>% str_remove("_X1$")
    )
  ) %>%
  # Join your key to this table
  left_join(user_provided_key, by = "base_name") %>%
  # Create the final clean_name using a set of rules
  mutate(
    clean_name = case_when(
      str_ends(original_name, "_X1") ~ clean_base_name,
      
      str_starts(original_name, "na_ind_") ~ str_c("Missingness Indicator for: ", clean_base_name),
      
      TRUE ~ clean_base_name
    )
  ) %>%
  # Select and arrange the final columns
  select(original_name, clean_name)

# stage one feature importance final table ----
table1_data <- s1_importance %>%
  # Join with the lookup table to get the clean names
  left_join(variable_lookup_table, by = c("Variable" = "original_name")) %>%
  # Keep only the essential columns
  select(Feature = clean_name, Importance) %>%
  # Remove duplicate features that might arise from dummy variables
  distinct(Feature, .keep_all = TRUE) %>%
  # Add rank and select the top 10
  mutate(Rank = row_number()) %>%
  slice_head(n = 10) %>%
  select(Rank, Feature, Importance)


# Whole cascade synthesised performance ----
s1_imp_renamed <- s1_importance %>% rename(s1_score = Importance)
prob_imp_renamed <- probable_importance %>% rename(prob_score = Importance)
gt24h_imp_renamed <- gt24h_importance %>% rename(gt24h_score = Importance)
nonsevere_imp_renamed <- nonsevere_importance %>% rename(nonsevere_score = Importance)
onset_within_24h_imp_renamed <- onset_within_24h_importance %>% rename(onset_within_24h_score = Importance)

# all importance data frames in one wide data frame
all_importances_wide <- reduce(
  list(s1_imp_renamed, prob_imp_renamed, gt24h_imp_renamed, nonsevere_imp_renamed, onset_within_24h_imp_renamed),
  full_join,
  by = "Variable"
) %>%
  # If a feature wasn't in a model's list, its importance is NA. We replace NA with 0.
  replace(is.na(.), 0) %>%
  # Join with your lookup table to get the clean, interpretable names
  left_join(variable_lookup_table, by = c("Variable" = "original_name")) %>%
  # Keep only the essential columns for the final table
  select(Feature = clean_name, `Stage One` = s1_score, `S2: Probable-type Outcome Specialist` = prob_score,
         `S2: Onset Greater Than 24h Specialist` = gt24h_score, `S2: Non Severe Outcome Specialist` = nonsevere_score,
         `S2: Onset Within 24h Specialist` = onset_within_24h_score)

# mean importance score and rank the features
table3_data_mean_score <- all_importances_wide %>%
  mutate(
    # the average (mean of squares) importance for Stage 2
    stage_2_mean_sq_score = (`S2: Probable-type Outcome Specialist`^2 + `S2: Onset Greater Than 24h Specialist`^2 +
                               `S2: Non Severe Outcome Specialist`^2 + `S2: Onset Within 24h Specialist`^2) / 4,
    
    # final 50/50 weighted RMS score
    # (this squares the S1 score, averages it with the S2 mean-of-squares, and takes the square root)
    `Mean Importance Score` = sqrt(
      (`Stage One`^2 + stage_2_mean_sq_score) / 2
    )
  ) %>%
  select(-stage_2_mean_sq_score) %>%
  # Keep only the essential columns for the final table
  select(`Feature`, `Mean Importance Score`) %>%
  # Remove duplicates that arise from dummy variables mapping to the same clean name
  group_by(Feature) %>%
  summarise(`Mean Importance Score` = max(`Mean Importance Score`)) %>%
  ungroup() %>%
  # Order the features from most to least important
  arrange(desc(`Mean Importance Score`)) %>%
  # Add a final Rank column
  mutate(Rank = row_number()) %>%
  select(Rank, Feature, `Mean Importance Score`)


write.csv(all_importances_wide, "Annexes and Supplementary Materials/All_Feature_Importances_FC.csv")
write.csv(table3_data_mean_score, "Annexes and Supplementary Materials/Weighted_Feature_Importance_FC.rds")
saveRDS(table3_data_mean_score[1:15,], "Visualisations/Tables/weighted_feature_importance_top_15.rds")
saveRDS(variable_lookup_table, "Visualisations/Tables/variable_lookup.rds")




