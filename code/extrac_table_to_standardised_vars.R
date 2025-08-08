# Set up ----
library(tidyverse)
library(janitor)
library(glue)
library(readr)
library(purrr)
library(fuzzyjoin) 
library(rlang)

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

# Reading extraction table data ----
tab  <- read_csv("final_extraction_table.csv") %>%
  clean_names() %>%
  mutate(paper_id = if("x1" %in% names(.)) x1 else row_number())

split_vars <- function(x) {
  x %>% str_split(";\\s*") %>%
    map(~ .x[.x != "" & .x != "No" & !is.na(.x)])
}

required_cols <- c("exposures", "mediators", "unmeasured_biases", "outcomes")
missing_cols <- setdiff(required_cols, names(tab))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in input CSV:", paste(missing_cols, collapse=", ")))
}

var_long <- tab %>%
  transmute(
    paper_id,
    exposures         = split_vars(exposures),
    mediators         = split_vars(mediators),
    unmeasured_biases = split_vars(unmeasured_biases),
    outcomes          = split_vars(outcomes)
  ) %>%
  pivot_longer(-paper_id, names_to = "role", values_to = "vars") %>%
  unnest(vars) %>%
  filter(!is.na(vars) & vars != "") %>%
  mutate(
    vars = str_trim(vars),
    role = recode(role,
                  exposures         = "exposure",
                  mediators         = "mediator",
                  unmeasured_biases = "unobserved_bias",
                  outcomes          = "outcome")
  ) %>%
  select(paper_id, variable_raw = vars, role) %>%
  arrange(paper_id, role, variable_raw)

# Calculating paper-level counts ----
# Summarising and reshaping counts of each variable by role, computing total mentions, and ordering by frequency
variables_summary <- var_long %>%
  count(variable_raw, role, name = "n_papers") %>%
  pivot_wider(names_from = role,
              values_from = n_papers,
              values_fill = 0) %>%
  {
    current_roles <- c("exposure", "mediator", "unobserved_bias", "outcome")
    existing_roles_in_data <- intersect(current_roles, names(.))
    missing_roles_for_sum <- setdiff(current_roles, existing_roles_in_data)
    for(m_role in missing_roles_for_sum) {
      .[m_role] <- 0
    }
    .
  } %>%
  mutate(total_mentions = rowSums(select(., any_of(c("exposure", "mediator", "unobserved_bias", "outcome"))), na.rm = TRUE)) %>%
  arrange(desc(total_mentions))

# Creating a cleaned, lowercase version of the raw variable names
vars <- variables_summary %>%
  mutate(raw_clean = variable_raw %>%
           str_to_lower() %>%
           str_replace_all("[^a-z0-9]+", " ") %>%
           str_squish())

# Manual variable lookup table creation + mapping ----
# Defining a manual lookup table mapping regex patterns to standardized variable names and categories
lookup_manual <- tribble(
  ~pattern, ~variable_std, ~var_category,
  "(?i)crp[\\s/-]?to[\\s/-]?albumin|\\bCAR\\b", "car", "derived_ratio",
  "(?i)crp[\\s/-]?to[\\s/-]?platelet|c[-\\s]?platelet\\s?ratio|\\bCPR\\b", "cpr", "derived_ratio",
  "(?i)lymphocyte.*crp|\\bLCR\\b", "lcr", "derived_ratio",
  "(?i)vital\\s?signs?|vitals\\b|clinical\\s?signs?|other\\s?clinical\\s?signs?", "vital_signs", "clinical_sign",
  "(?i)c-reactive protein|\\bcrp\\b", "crp", "lab_biomarker",
  "(?i)interleukin[ -]?6|\\bil-6\\b", "il_6", "lab_biomarker",
  "(?i)interleukin[ -]?8|\\bil-8\\b", "il_8", "lab_biomarker",
  "(?i)age", "age", "demographic",
  "(?i)sex|gender", "sex", "demographic",
  "(?i)temperature", "temperature", "clinical_sign",
  "(?i)body mass index|\\bbmi\\b", "bmi", "anthropometric",
  "(?i)systolic blood pressure|\\bsbp\\b", "sbp", "clinical_sign",
  "(?i)mortality|death", "mortality", "outcome",
  "(?i)length of stay|\\blos\\b|picu_stay_length", "length_of_stay", "outcome",
  "(?i)selection bias", "selection_bias", "study_bias",
  "(?i)shock", "shock", "clinical_condition",
  "(?i)single[- ]?(centre|center|site|country|hospital|setting)|singlecentre", "single_center_study", "study_bias",
  "(?i)small[- ]?(sample|cohort|size|number)|limited[- ]?(sample|size)|few[- ]?cases|small sbi", "small_sample_size", "study_bias",
  "(?i)hemogram|cbc[- ]?indices", "cbc_indices", "lab_biomarker",
  "(?i)apgar", "apgar_score", "clinical_sign",
  "(?i)blood[- ]?urea[- ]?nitrogen|\\bbun\\b", "bun", "lab_biomarker",
  "(?i)\\b(nt-?)?pro-?bnp\\b|brain natriuretic peptide|\\bbnp\\b", "bnp_ntprobnp", "lab_biomarker",
  "(?i)(l-)?fabp|fatty acid binding protein", "fabp", "lab_biomarker",
  "(?i)thrombomodulin", "thrombomodulin", "lab_biomarker",
  "(?i)bilirubin", "bilirubin", "lab_biomarker",
  "(?i)poor[- ]?(clinical[- ]?)?outcome|poor[- ]?prognosis|poor[- ]?neurological[- ]?outcome", "poor_outcome", "outcome",
  "(?i)incomplete[- ]?(records|data)|partial[- ]?data", "incomplete_data", "study_bias",
  "(?i)limited[- ]?external[- ]?validity|no[- ]?external[- ]?validation", "limited_external_validity", "study_bias",
  "(?i)respiratory[- ]?failure", "respiratory_failure", "clinical_condition",
  "(?i)acute[- ]?kidney[- ]?injury|\\baki\\b", "aki", "clinical_condition",
  "(?i)cytokine(s)?|cytokine[- ]?(panel|response|storm)", "cytokine_general", "lab_biomarker",
  "(?i)procalcitonin|\\bpct\\b", "procalcitonin", "lab_biomarker",
  "(?i)rural[- ]?residence|urban[- ]?residence", "residence_type", "demographic",
  "(?i)central[- ]?line", "central_line", "treatment",
  "(?i)vitamin[- ]?d", "vitamin_d", "lab_biomarker",
  "(?i)westley score", "westley_score", "clinical_score"
)

# Joining cleaned variable names to the lookup table, selecting the first match, and filling missing standardized names and categories
vars1 <- vars %>%
  regex_left_join(lookup_manual, by = c(raw_clean = "pattern")) %>%
  group_by(across(c(names(vars)))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    variable_std = coalesce(variable_std, raw_clean),
    var_category = coalesce(var_category, "unreviewed")
  ) %>%
  select(-pattern)


# Defining manual overrides mapping exact cleaned terms to standardised categories
manual_overrides <- tribble(
  ~term_to_match,              ~var_category,
  "ari",                       "clinical_condition",
  "lar",                       "derived_ratio",
  "par",                       "lab_biomarker",
  "measurement error",         "measurement_bias",
  "hypoxia",                   "clinical_sign",
  "altered consciousness",     "clinical_sign",
  "anc",                       "lab_biomarker",
  "alt",                       "lab_biomarker",
  "ast",                       "lab_biomarker",
  "creatinine",                "lab_biomarker",
  "ldh",                       "lab_biomarker",
  "pdw",                       "lab_biomarker",
  "ph",                        "lab_biomarker",
  "sodium",                    "lab_biomarker",
  "dexamethasone",             "treatment",
  "household crowding",        "demographic",
  "none",                      "none_reported",
  "term birth",                "demographic",
  "acute circulatory failure", "clinical_condition",
  "anaemia",                   "clinical_condition",
  "bradycardia",               "clinical_sign",
  "comorbidity",               "clinical_condition",
  "diarrhoea",                 "clinical_sign",
  "leukopenia",                "lab_biomarker",
  "who danger signs",          "clinical_sign",
  "csf glucose",               "lab_biomarker",
  "dni",                       "lab_biomarker",
  "igg igm",                   "lab_biomarker",
  "il 6 testing implementation","measurement_bias",
  "serum tat",                 "lab_biomarker",
  "cd27",                      "molecular_omics",
  "cd45ro",                    "molecular_omics",
  "host immune response pathways","molecular_omics",
  "snp instruments",           "molecular_omics",
  "ventilatory support",       "treatment",
  "birth order",               "demographic",
  "caregiver education",       "demographic",
  "facility delivery",         "demographic",
  "wealth index",              "demographic",
  "prom",                      "clinical_condition",
  "vur",                       "clinical_condition",
  "convenience sampling",      "study_bias",
  "limited controls",          "study_bias",
  "inclusion criteria",        "study_design_characteristic",
  "data quality",              "study_bias",
  "emr data quality",          "study_bias",
  "mr pleiotropy",             "study_bias",
  "burn etiology",             "clinical_condition",
  "micronutrient supplementation","nutrition",
  "yingyangbao",               "nutrition",
  "helminths",                 "helminth_infection",
  "gnpda2",                    "molecular_omics",
  "gyg1",                      "molecular_omics",
  "mettl7b",                   "molecular_omics",
  "ms4a4a",                    "molecular_omics",
  "nsun7",                     "molecular_omics",
  "olfm4 neutrophils",         "lab_biomarker",
  "stom",                      "clinical_condition",
  "renal perfusion",           "imaging",
  "psp levels",                "lab_biomarker",
  "admission time block",      "temporal_factor",
  "amputations skin grafts",    "treatment",
  "ca 6h",                     "lab_biomarker",
  "caas",                      "lab_biomarker",
  "clinical diagnoses",        "clinical_finding_general",
  "deep breathing",            "clinical_sign",
  "estimated gfr accuracy vs inulin","measurement_bias",
  "evolving protocols",        "study_design_characteristic",
  "hladr",                     "molecular_omics",
  "ifn score",                 "lab_biomarker",
  "immune modulation by vitamin d","pathophysiology",
  "improved vs unimproved housing features","demographic",
  "infant clinical indicators", "clinical_finding_general",
  "inflammatory elevation",    "lab_biomarker_general",
  "inflammatory response to burns","pathophysiology",
  "lab availability",          "study_bias",
  "limited virology",          "study_bias",
  "low leukopenia prevalence", "study_finding",
  "low sensitivity",           "measurement_bias",
  "mrproadm",                  "lab_biomarker",
  "no bleeding endpoints",     "study_design_characteristic",
  "perfusion signs",           "clinical_sign",
  "pic",                       "lab_biomarker",
  "placenta changes",          "pathophysiology",
  "pmnl percent",              "lab_biomarker",
  "possible cons contamination","microbiology_finding",
  "pretermterm mix",           "study_bias",
  "prior episodes",            "clinical_history",
  "regressionbased scores",    "clinical_score",
  "renal ischemia",            "clinical_condition",
  "retn",                      "lab_biomarker",
  "secondary analysis",        "study_design_characteristic",
  "serum adamts13 level",      "lab_biomarker",
  "serum dglucan",             "lab_biomarker",
  "subtype stratification",    "study_design_characteristic",
  "tpaic levels",              "lab_biomarker",
  "unmeasured comorbidities",  "study_bias"
) %>%
  rename(variable_std = term_to_match)

# Defining final mapping of raw variable strings to standardized variable names and categories
final_map <- tribble(
  ~var_space,                                     ~variable_std,                                  ~var_category,
  "loss to follow up",                            "loss_to_follow_up",                            "study_bias",
  "case mix differences",                         "case_mix_differences",                         "study_bias",
  "early late",                                   "early_late",                                   "temporal_factor",
  "i t ratio",                                    "i_t_ratio",                                    "derived_ratio",
  "incomplete follow up",                         "incomplete_follow_up",                         "study_bias",
  "lack of longitudinal follow up",               "lack_of_longitudinal_follow_up",               "study_bias",
  "laz 4",                                        "laz_4",                                        "nutrition",
  "monocytes 6 5",                                "monocytes_6_5",                                "lab_biomarker",
  "pai 1 as mediator of coagulopathy discussed",  "pai_1_mediator_coagulopathy",                  "lab_biomarker",
  "plasma pai 1 level",                           "plasma_pai_1",                                 "lab_biomarker",
  "plasma strem 1 level on diagnosis",            "plasma_strem_1_diagnosis",                     "lab_biomarker",
  "renal urological history",                     "renal_urological_history",                     "clinical_condition",
  "serum apoa v level on admission",              "serum_apoa_v_admission",                       "lab_biomarker",
  "serum ifn",                                    "serum_ifn",                                    "lab_biomarker",
  "acute kidney injury",                          "aki",                                          "clinical_condition",
  "development of nards during admission",        "ards_development",                             "clinical_condition",
  "burn size tbs dba dsi",                        "burn_size",                                    "clinical_sign",
  "early renal resistive index",                  "renal_resistive_index",                        "imaging",
  "eos risk calculator score",                    "eos_risk_score",                               "clinical_score",
  "icu admission",                                "icu_admission",                                "outcome",
  "picu admission",                               "picu_admission",                               "outcome",
  "admission fibrinogen categories",              "fibrinogen_admission_categories",              "lab_biomarker",
  "4 anticonvulsants",                            "anticonvulsants",                              "treatment",
  "circulating micrornas",                        "circulating_mirna",                            "molecular_omics",
  "development of respiratory failure",           "respiratory_failure_development",              "clinical_condition",
  "severe or complicated cap",                    "severe_cap",                                   "clinical_condition",
  "56 plasma protein biomarkers",                 "plasma_protein_panel_56",                      "lab_biomarker",
  "labs",                                         "lab_panel_general",                            "lab_biomarker_general",
  "csf protein 2 g l",                            "csf_protein",                                  "lab_biomarker",
  "pni 52 4",                                     "pni_score",                                    "derived_ratio",
  "serial supar levels",                          "supar_serial",                                 "lab_biomarker",
  "serum ige categories",                         "serum_ige_categories",                         "lab_biomarker",
  "mv within 24h",                                "mechanical_ventilation_24h",                   "treatment",
  "neurocognitive scores over 24 months",         "neurocognitive_score_24m",                     "outcome",
  "neurologic disability at discharge",           "neurologic_disability_discharge",              "outcome",
  "plasma chi3l1 levels daily for 4 days",        "plasma_chi3l1_longitudinal",                   "lab_biomarker",
  "plasma metabolite panel glycerophospholipids sphingolipids", "plasma_metabolite_panel_lipids", "lab_biomarker",
  "serum metabolite panel 55 metabolites",        "serum_metabolite_panel_55",                    "lab_biomarker"
)


# Defining keyword mappings for category assignment and function to determine category based on keywords ----
keyword_map <- list(
  lab_biomarker    = c("crp", "pct", "wbc", "nlr", "plt", "lactate", "il\\d+", "ddimer", "d-dimer",
                       "rdw", "mpv", "plr", "trem1", "ngal", "presepsin", "resistin", "retn",
                       "ferritin", "hepcidin", "stfr", "tau", "amyloid", "hbp", "fabp", "supar",
                       "progranulin", "cystatin", "egfr", "esr", "vegfa", "pai-?1", "strem", "albumin", "platelet",
                       "neutrophil", "bilirubin", "\\bbun\\b", "\\b(nt-?)?pro-?bnp\\b", "thrombomodulin", "procalcitonin",
                       "alt", "ast", "creatinine", "sodium", "ldh", "ige", "glucose", "protein", "caas", "pmnl", "tpaic",
                       "adamts13", "dglucan", "psp", "mrproadm", "ifn", "pic", "biomarker"),
  clinical_sign    = c("fever", "temperature", "hypothermia", "coma", "gcs", "blantyre",
                       "desaturation", "hypoperfusion", "apnea", "vital sign", "vitals", "clinical sign",
                       "hypoxia", "apgar", "bradycardia", "diarrhoea", "danger sign", "breathing", "perfusion"),
  clinical_condition= c("pneumonia", "uti", "sepsis", "shock", "hfmd", "meningitis",
                        "bronchiolitis", "malaria", "covid", "bsi", "\\bnec\\b", "aki",
                        "thrombocytopenia", "mods", "ari", "anaemia", "helminth", "prom", "vur", "stomatitis",
                        "respiratory failure", "circulatory failure", "ischemia", "ards", "cap"),
  outcome            = c("mortality", "death", "length_of_stay", "los", "readmission",
                         "ventilation", "progression", "severity", "hospitalization", "disability",
                         "icu_admission", "picu_admission", "poor_outcome", "prognosis", "neurocognitive"),
  treatment          = c("antibiotic", "ivig", "steroids", "cpap",
                         "oxygen", "ventilation", "transfusion", "central_line", "procedure",
                         "anticonvulsant", "dexamethasone", "support", "graft", "amputation"),
  study_bias         = c("bias", "heterogeneity", "missing", "limitation", "power",
                         "external_validity", "dataset", "misclassification", "selection", "confounding",
                         "follow up", "followup", "single_center", "single_site", "small_sample", "small_cohort",
                         "validation", "control", "sampling", "quality", "pleiotropy", "contamination", "prevalence",
                         "availability", "endpoint"),
  demographic        = c("age", "sex", "gender", "prematurity", "birth_weight", "maternal",
                         "cesarean", "term_birth", "crowding", "education", "residence", "wealth", "order", "housing", "delivery"),
  nutrition          = c("nutrition", "diet", "feeding", "stunting", "wasting", "underweight", "laz", "supplement"),
  molecular_omics  = c("gene", "mir", "transcript", "expression", "omics", "microarray", "mrna", "microrna",
                       "klrg1", "tdrd9", "loc728401", "prkacb", "pros1", "snp", "marker", "pathway", "cd\\d+", "hladr"),
  microbiology_test= c("culture", "pathogen", "viral", "pcr", "bacteremia", "viremia", "ev71", "virology"),
  imaging            = c("ultrasound", "xray", "radiograph", "ct", "mri", "renal_perfusion", "resistive_index"),
  measurement_bias   = c("assay", "timing", "batch", "device", "handling", "cutoff", "error", "implementation", "accuracy", "sensitivity"),
  temporal_factor    = c("season", "time_block", "early/late", "early_late", "admission_time"),
  derived_ratio      = c("ratio", "\\bcar\\b", "\\bcpr\\b", "\\blcr\\b", "\\bpni\\b"),
  anthropometric     = c("weight", "height", "bmi", "muac", "z-score"),
  clinical_score     = c("score", "calculator", "criteria"),
  pathophysiology    = c("response", "modulation", "changes", "cascade", "etiology"),
  study_design_characteristic = c("protocol", "analysis", "stratification", "criteria"),
  clinical_history   = c("episode", "history"),
  microbiology_finding = c("contamination"),
  lab_biomarker_general = c("panel", "labs", "elevation"),
  clinical_finding_general = c("indicator", "diagnoses")
)

cat_by_keyword <- function(x_raw_clean) {
  if (length(x_raw_clean) != 1 || !is.character(x_raw_clean)) {
    return("other")
  }
  for (cat in names(keyword_map)) {
    if (str_detect(x_raw_clean, regex(paste(keyword_map[[cat]], collapse = "|"), ignore_case = TRUE)))
      return(cat)
  }
  "other"
}

# Priority mapping ----
# Creating high-priority overrides from exact cleaned terms in manual_overrides
lookup_specific_overrides <- manual_overrides %>%
  mutate(
    pattern = str_c("^(?:", str_replace_all(fixed(variable_std), "([.^$|()?*+{}\\[\\]])", "\\\\\\1"), ")$"),
    source = "manual_override_raw_clean",
    priority = 1
  )

# Creating high-priority mappings from exact slugs in final_map
lookup_final_map_aliases <- final_map %>%
  transmute(
    pattern = str_c("^(?:", str_replace_all(fixed(variable_std), "([.^$|()?*+{}\\[\\]])", "\\\\\\1"), ")$"),
    variable_std,
    var_category,
    source = "final_map_alias_to_slug",
    priority = 1
  )

# Combining highest-priority direct lookup sources
lookup_high_priority_direct <- bind_rows(lookup_specific_overrides, lookup_final_map_aliases)

# Creating manual regex-based lookup source with priority 2
lookup_manual_src <- lookup_manual %>%
  mutate(source = "manual_regex", priority = 2)

# Identifying unreviewed terms for automatic keyword-based categorization
auto_candidates <- vars1 %>%
  filter(var_category == "unreviewed") %>%
  distinct(raw_clean, variable_raw)

# Building automatic lookup source from keyword matches for priority 3
lookup_auto_src <- auto_candidates %>%
  mutate(
    variable_std_auto = raw_clean,
    var_category_auto = map_chr(raw_clean, cat_by_keyword),
    pattern_auto = str_c("^(?:", str_replace_all(fixed(raw_clean), "([.^$|()?*+{}\\[\\]])", "\\\\\\1"), ")$")
  ) %>%
  filter(var_category_auto != "other" & var_category_auto != "unreviewed") %>%
  select(
    pattern = pattern_auto,
    variable_std = variable_std_auto,
    var_category = var_category_auto
  ) %>%
  mutate(source = "auto_keyword", priority = 3)

# Combining all lookup sources into an intermediate table
combined_lookup_intermediate <- bind_rows(
  lookup_high_priority_direct,
  lookup_manual_src,
  lookup_auto_src
)

# Refining combined lookup by prioritizing and deduplicating patterns
final_lookup <- combined_lookup_intermediate %>%
  group_by(variable_std) %>%
  filter(!(var_category %in% c("other", "unreviewed") & any(!var_category %in% c("other", "unreviewed", NA)))) %>%
  ungroup() %>%
  arrange(pattern, priority, ifelse(var_category %in% c("other", "unreviewed", NA), 1, 0)) %>%
  distinct(pattern, .keep_all = TRUE) %>%
  select(pattern, variable_std, var_category, source, priority)


# Combining cleaned variable names with lookup table, selecting highest-priority match, and creating final slugs and categories
vars_classified <- vars %>%
  select(-any_of(c("variable_std", "var_category", "pattern"))) %>%
  regex_left_join(final_lookup, by = c(raw_clean = "pattern")) %>%
  group_by(across(c(names(vars)))) %>%
  arrange(priority, ifelse(is.na(priority), 1, 0), ifelse(var_category %in% c("other", "unreviewed", NA), 1, 0)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    variable_std_final = coalesce(variable_std, raw_clean),
    var_category_final = coalesce(var_category, "unreviewed"),
    variable_std_final = str_to_lower(str_replace_all(variable_std_final, "[^a-z0-9_]+", "_")) %>%
      str_replace_all("_+", "_") %>%
      str_remove("^_|_$")
  ) %>%
  select(
    variable_raw, raw_clean,
    any_of(c("exposure", "mediator", "unobserved_bias", "outcome", "total_mentions")),
    variable_std = variable_std_final,
    var_category = var_category_final,
    match_source = source,
    match_pattern = pattern
  )

# Aggregating standardized variables into a summary table with counts and metadata ----
variables_std_summary <- vars_classified %>%
  group_by(variable_std, var_category) %>%
  summarise(
    exposure = sum(exposure, na.rm = TRUE),
    mediator = sum(mediator, na.rm = TRUE),
    unobserved_bias = sum(unobserved_bias, na.rm = TRUE),
    outcome = sum(outcome, na.rm = TRUE),
    total_mentions = sum(total_mentions, na.rm = TRUE),
    raw_synonyms = paste(sort(unique(variable_raw)), collapse = "; "),
    match_sources_used = paste(sort(unique(na.omit(match_source))), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_mentions))

# Identifying variables still classified as unreviewed or other for sanity checking ----
remaining_unreviewed <- filter(variables_std_summary, var_category == "unreviewed")
remaining_other <- filter(variables_std_summary, var_category == "other")


# Printing all variables still marked as unreviewed
print(remaining_unreviewed, n = Inf)

# Spot-checking remaining unreviewed, miscategorised, and practical duplicate variables ----
# Defining manual corrections for specific unreviewed variables
hand_patch <- tribble(
  ~variable_std,                                                          ~new_std,                                   ~new_category,
  "admission_time_block",                                                 "admission_time_block",                     "demographic",
  "birth_weight",                                                         "birth_weight",                             "clinical_sign",
  "clinical_diagnoses",                                                   "clinical_diagnoses",                       "clinical_condition",
  "dic_score_components",                                                 "dic_score_components",                     "clinical_sign",
  "eos_risk_calculator_score",                                            "eos_risk_score",                           "clinical_sign",
  "evolving_protocols",                                                   "evolving_protocols",                       "study_bias",
  "fibrogen_admisson_categories",                                         "fibrogen_admisson_categories",             "lab_biomarker",
  "helminths",                                                            "helminths",                                "microbiology_test",
  "il_8_inflammatory_cascade",                                            "il8_inflammatory_cascade",                 "molecular_omics",
  "immune_modulation_by_vitamin_d",                                       "immune_modulation_vitd",                   "treatment",
  "inclusion_criteria",                                                   "inclusion_criteria",                       "study_bias",
  "infant_clinical_indicators",                                           "infant_clinical_indicators",               "clinical_sign",
  "inflammatory_elevation",                                               "inflammatory_elevation",                   "lab_biomarker",
  "inflammatory_response_to_burns",                                       "inflammatory_response_burns",              "clinical_sign",
  "labs",                                                                 "labs",                                     "lab_biomarker",
  "low_leukopenia_prevalence",                                            "low_leukopenia_prevalence",                "study_bias",
  "no_bleeding_endpoints",                                                "no_bleeding_endpoints",                    "study_bias",
  "none",                                                                 "none_reported",                            "study_bias",
  "placenta_changes",                                                     "placenta_changes",                         "clinical_sign",
  "plasma_metabolite_panel_glycerophospholipids_sphingolipids",             "plasma_metab_glyc_sphing",                 "molecular_omics",
  "possible_cons_contamination",                                          "possible_contamination",                   "study_bias",
  "prior_episodes",                                                       "prior_episodes",                           "clinical_condition",
  "regressionbased_scores",                                               "regression_scores",                        "clinical_sign",
  "renal_urological_history",                                             "renal_urological_history",                 "clinical_condition",
  "season",                                                               "season",                                   "demographic",
  "seasonality",                                                          "seasonality",                              "demographic",
  "secondary_analysis",                                                   "secondary_analysis",                       "study_bias",
  "serum_metabolite_panel_55_metabolites",                                "serum_metab_55",                           "molecular_omics",
  "subtype_stratification",                                               "subtype_stratification",                   "study_bias",
  "westley_score",                                                        "westley_score",                            "clinical_sign",
  "admission_serum_lactate",                                              "lactate",                                  "lab_biomarker",
  "admission_serum_lactate_level",                                        "lactate",                                  "lab_biomarker",
  "blood_lactate",                                                        "lactate",                                  "lab_biomarker",
  "6h_lactate_clearance",                                                 "lactate",                                  "lab_biomarker",
  "hscrp",                                                                "crp",                                      "lab_biomarker",
  "admission_hscrp",                                                      "crp",                                      "lab_biomarker",
  "admission_plasma_strem1",                                              "plasma_strem1",                            "lab_biomarker",
  "plasma_strem_1_level_on_diagnosis",                                    "plasma_strem1",                            "lab_biomarker",
  "strem1",                                                               "plasma_strem1",                            "lab_biomarker",
  "admission_presepsin_level",                                            "presepsin",                                "lab_biomarker",
  "plasma_presepsin",                                                     "presepsin",                                "lab_biomarker",
  "cord_presepsin_level",                                                 "presepsin",                                "lab_biomarker",
  "cord_and_peripheral_presepsin",                                        "presepsin",                                "lab_biomarker",
  "pai_1_as_mediator_of_coagulopathy_discussed",                            "plasma_pai_1_level",                       "lab_biomarker",
  "admission_nlr_strata",                                                 "nlr",                                      "lab_biomarker",
  "nlr_5",                                                                "nlr",                                      "lab_biomarker",
  "neutrophiltolymphocyte_ratio",                                         "nlr",                                      "lab_biomarker",
  "admission_nlr_plr_systemic_immune_inflammatory_index_wbc",             "nlr",                                      "lab_biomarker",
  "platelet_neutrophil_and_lymphocyte_counts",                              "plr",                                      "lab_biomarker",
  "plt",                                                                  "platelet_count",                           "lab_biomarker",
  "plt_6h",                                                               "platelet_count",                           "lab_biomarker",
  "monocytes_greater_than_6.5_percent",                                   "monocytes_gt_6_5_percent",                 "lab_biomarker",
  "monocytes_6_5",                                                        "monocytes_gt_6_5_percent",                 "lab_biomarker",
  "ldh_albumin_ratio_on_admission",                                       "ldh_albumin_ratio",                        "derived_ratio",
  "breastfeeding",                                                        "breastfeeding",                            "nutrition",
  "wasting",                                                              "wasting",                                  "nutrition",
  "asthma_history",                                                       "asthma_history",                           "clinical_condition",
  "need_for_picu_vs_ward_outpatient_care_among_sepsis_cases",             "need_for_picu_vs_ward_outpatient_care_among_sepsis_cases", "outcome",
  "pediatric_sepsis_vs_healthy_status",                                   "pediatric_sepsis_vs_healthy_status",       "study_design_characteristic",
  "hospitalized_bronchiolitis_vs_healthy_controls",                       "hospitalized_bronchiolitis_vs_healthy_controls", "study_design_characteristic",
  "mycoplasma_pneumoniae_pneumonia_vs_healthy_infectious_controls",       "mycoplasma_pneumoniae_pneumonia_vs_healthy_infectious_controls", "study_design_characteristic",
  "measurement_bias_of_biomarkers",                                       "measurement_bias_of_biomarkers",           "measurement_bias",
  "biomarker_timing",                                                     "biomarker_timing",                         "measurement_bias",
  "biomarker_timing_variability",                                         "biomarker_timing_variability",             "measurement_bias",
  "biomarker_assay_variability",                                          "biomarker_assay_variability",              "measurement_bias",
  "missing_biomarkers",                                                   "missing_biomarkers",                       "measurement_bias",
  "protein_microarray_of_9345_autoantibodies_8_selected",                 "protein_microarray_autoantibodies",        "molecular_omics",
  "inflammatory_protein_mediators",                                       "inflammatory_protein_mediators",           "molecular_omics",
  "invasive_fungal_disease",                                              "invasive_fungal_disease",                  "clinical_condition",
  "confounding_by_healthcare_access",                                     "confounding_by_healthcare_access",         "study_bias",
  "acute_circulatory_failure",                                            "shock",                                    "clinical_condition",
  "bacterial_vs_viral_pneumonia",                                         "pneumonia",                                "clinical_condition",
  "complicated_pneumonia",                                                "pneumonia",                                "clinical_condition",
  "etiologic_class_bacterial_vs_viral_vs_malarial_pneumonia",             "pneumonia",                                "clinical_condition",
  "pneumonia_with_consolidation",                                         "pneumonia",                                "clinical_condition",
  "radiographic_pneumonia_on_chest_xray",                                 "pneumonia",                                "clinical_condition",
  "severe_or_complicated_cap",                                            "pneumonia",                                "clinical_condition",
  "who_severe_rsv_communityacquired_pneumonia",                           "pneumonia",                                "clinical_condition",
  "incident_malaria_episodes_over_6months",                               "malaria_infection",                        "clinical_condition",
  "malaria_recurrence_over_6_months",                                     "malaria_infection",                        "clinical_condition",
  "disease_class_malaria_vs_bacteremia",                                  "malaria_infection",                        "clinical_condition",
  "malaria_vs_bacterial_bsi_classification",                              "malaria_infection",                        "clinical_condition",
  "clinical_early_onset_neonatal_sepsis_ceons",                           "early_onset_neonatal_sepsis",              "clinical_condition",
  "clinical_earlyonset_neonatal_sepsis",                                  "early_onset_neonatal_sepsis",              "clinical_condition",
  "earlyonset_neonatal_sepsis",                                           "early_onset_neonatal_sepsis",              "clinical_condition",
  "earlyonset_neonatal_sepsis_diagnosis",                                 "early_onset_neonatal_sepsis",              "clinical_condition",
  "early_onset_neonatal_sepsis",                                          "early_onset_neonatal_sepsis",              "clinical_condition",
  "earlyonset_sepsis",                                                    "early_onset_sepsis",                       "clinical_condition",
  "earlyonset_sepsis_diagnosis",                                          "early_onset_sepsis",                       "clinical_condition",
  "earlyonset_sepsis_in_preterm_infants",                                 "early_onset_sepsis",                       "clinical_condition",
  "culture_confirmed_early_onset_neonatal_sepsis_72_h",                   "culture_confirmed_neonatal_sepsis",        "clinical_condition",
  "culture_positive_late_onset_sepsis_workup",                            "culture_confirmed_neonatal_sepsis",        "clinical_condition",
  "culture_proven_neonatal_sepsis",                                       "culture_confirmed_neonatal_sepsis",        "clinical_condition",
  "cultureconfirmed_bacteremic_sepsis_early_vs_late_onset",               "culture_confirmed_neonatal_sepsis",        "clinical_condition",
  "cultureconfirmed_neonatal_sepsis",                                     "culture_confirmed_neonatal_sepsis",        "clinical_condition",
  "culturepositive_earlyonset_sepsis",                                    "culture_confirmed_neonatal_sepsis",        "clinical_condition",
  "culturepositive_sepsis",                                               "culture_confirmed_neonatal_sepsis",        "clinical_condition",
  "cultureproven_neonatal_bacterial_sepsis",                              "culture_confirmed_neonatal_sepsis",        "clinical_condition",
  "neonatal_sepsis_and_severe_sepsis_status",                             "neonatal_sepsis",                          "clinical_condition",
  "neonatal_sepsis_clinical_culture",                                     "neonatal_sepsis",                          "clinical_condition",
  "neonatal_sepsis_cultureproven_and_suspected",                          "neonatal_sepsis",                          "clinical_condition",
  "neonatal_sepsis_diagnosis",                                            "neonatal_sepsis",                          "clinical_condition",
  "neonatal_sepsis_overall_and_culturepositive",                          "neonatal_sepsis",                          "clinical_condition",
  "presence_and_severity_of_neonatal_sepsis",                             "neonatal_sepsis",                          "clinical_condition",
  "presence_of_neonatal_sepsis",                                          "neonatal_sepsis",                          "clinical_condition",
  "occurrence_of_neonatal_sepsis",                                        "neonatal_sepsis",                          "clinical_condition",
  "development_of_sepsis",                                                "sepsis_development",                       "clinical_condition",
  "progression_to_sepsis",                                                "sepsis_development",                       "clinical_condition",
  "assay_variability",                                                    "measurement_variability",                  "measurement_bias",
  "device_measurement_variation",                                         "measurement_variability",                  "measurement_bias",
  "measurementtiming_variability",                                        "measurement_variability",                  "measurement_bias",
  "imaging_variability",                                                  "measurement_variability",                  "measurement_bias",
  "imaging_confirmation_variability",                                     "measurement_variability",                  "measurement_bias",
  "observer_variability_in_radiograph_reading",                           "measurement_variability",                  "measurement_bias",
  "operator_variability_of_ultrasound",                                   "measurement_variability",                  "measurement_bias",
  "sequencing_variability",                                               "measurement_variability",                  "measurement_bias",
  "technical_variability",                                                "measurement_variability",                  "measurement_bias",
  "timing_variability",                                                   "measurement_variability",                  "measurement_bias",
  "variable_assay_timing",                                                "measurement_variability",                  "measurement_bias",
  "variable_cutoffs",                                                     "measurement_variability",                  "measurement_bias",
  "variable_outcome_definitions",                                         "measurement_variability",                  "measurement_bias",
  "variable_reference_standards",                                         "measurement_variability",                  "measurement_bias",
  "variability_in_csf",                                                   "measurement_variability",                  "measurement_bias",
  "endemic_variability",                                                  "measurement_variability",                  "measurement_bias",
  "prettm_variability",                                                   "measurement_variability",                  "measurement_bias",
  "missing_variables",                                                    "measurement_variability",                  "measurement_bias",
  "monocyte_distribution_width",                                          "monocyte_distribution_width",              "lab_biomarker",
  "pathogen_distribution",                                                "pathogen_distribution",                    "microbiology_test",
  "neonatal_sepsis_vs_control_status",                                    "neonatal_sepsis_vs_control_status",        "study_design_characteristic",
  "neonatal_sepsis_vs_control_transcriptomic",                            "neonatal_sepsis_vs_control_transcriptomic", "study_design_characteristic",
  "sepsis_vs_control_vlbw_infants",                                       "sepsis_vs_control_vlbw_infants",           "study_design_characteristic",
  "sepsis_vs_noninfectious_sirs",                                         "sepsis_vs_noninfectious_sirs",             "study_design_characteristic",
  "prom_confounding",                                                     "prom_confounding",                         "study_bias",
  "single_institution",                                                   "single_institution",                       "study_bias",
  "case_fatality_in_cerebral_malaria",                                    "case_fatality_cerebral_malaria",           "outcome",
  "admission_plasma_il8",                                                 "il8",                                      "lab_biomarker",
  "plasma_il27",                                                          "il27",                                     "lab_biomarker",
  "admission_plasma_hbp",                                                 "hbp",                                      "lab_biomarker",
  "serum_hbp_quartiles",                                                  "hbp",                                      "lab_biomarker",
  "retn_gene_expression",                                                 "retn_gene_expression",                     "molecular_omics",
  "plasma_chi3l1_levels",                                                 "chi3l1",                                   "lab_biomarker",
  "plasma_ngal",                                                          "ngal",                                     "lab_biomarker",
  "urine_ngal",                                                           "ngal",                                     "lab_biomarker",
  "serial_supar_levels",                                                  "supar",                                    "lab_biomarker",
  "olfm4_neutrophils",                                                    "olfm4",                                    "lab_biomarker",
  "16_clinical_lab_variables_e_g_ldh_ckmb_glucose_vomiting",              "clinical_lab_variables_panel_16",          "lab_biomarker_general",
  "56_plasma_protein_biomarkers",                                         "protein_panel_56",                         "molecular_omics",
  "admission_lab_biomarkers_haemoglobin_prothrombin_time_routine_biochemistry", "basic_lab_panel",                      "lab_biomarker_general",
  "biomarkers",                                                           "lab_panel_general",                        "lab_biomarker_general",
  "biomarkers_only_at_admission",                                         "lab_panel_general",                        "lab_biomarker_general",
  "labs",                                                                 "lab_panel_general",                        "lab_biomarker_general",
  "ph_lactate_glucose_ast_alt_protein",                                   "blood_gas_biochem_panel",                  "lab_biomarker_general",
  "mpv_at_admission_and_72h",                                             "mpv",                                      "lab_biomarker",
  "wbc_6h",                                                               "wbc",                                      "lab_biomarker",
  "fever_in_past_2weeks",                                                 "recent_fever",                             "clinical_sign",
  "pctguided_antibiotic_algorithms",                                      "pct_guided_antibiotic_algorithm",          "treatment",
  "picu_admission",                                                       "picu_admission",                           "outcome",
  "picu_stay_length",                                                     "picu_stay_length",                         "outcome",
  "burn_size_tbs_dba_dsi",      "burn_size",  "clinical_sign",
  "early_renal_resistive_index",            "renal_resistive_index",              "imaging",
  "icu_admission", "icu_admission", "outcome",
  "admission_fibrinogen_categories", "fibrogen_admisson_categories", "lab_biomarker_general",
  "case_mix_differences", "heterogeneity", "study_bias",
  "early_late", "early_or_late_onset_neonatal_sepsis", "outcome",
  "few_eos_cases", "small_sample_size", "study_bias",
  "il_6", "il6", "lab_biomarker",
  "il_8_levels", "il8", "lab_biomarker",
  "monocytes_6_5", "monocytes_greater_than_6.5_percent", "lab_biomarker",
  "mv_within_24h", "mechanical_ventilation", "outcome",
  "invasive_mechanical_ventilation", "mechanical_ventilation", "outcome",
  "pai_1_as_mediator_of_coagulopathy_discussed", "plasma_pai_1_level", "lab_biomarker",
  "plasma_pai_1_level", "plasma_pai_1_level", "lab_biomarker",
  "sample_size", "small_sample_size", "study_bias",
  "serum_apoa_v_level_on_admission", "serum_apoa_v_level_on_admission", "lab_biomarker",
  "plasma_chi3l1_levels_daily_for_4_days", "plasma_chi3l1_levels", "lab_biomarker"
)

# Applying manual corrections and collapsing any resulting duplicates
variables_std_summary <- variables_std_summary %>%
  left_join(hand_patch, by = "variable_std") %>%
  mutate(
    variable_std = coalesce(new_std, variable_std),
    var_category = coalesce(new_category, var_category),
    match_sources_used = if_else(!is.na(new_std),
                                 "final manual patch",
                                 match_sources_used)
  ) %>%
  select(-new_std, -new_category) %>%
  group_by(variable_std, var_category) %>%
  summarise(
    exposure         = sum(exposure,         na.rm = TRUE),
    mediator         = sum(mediator,         na.rm = TRUE),
    unobserved_bias  = sum(unobserved_bias,  na.rm = TRUE),
    outcome          = sum(outcome,          na.rm = TRUE),
    total_mentions   = sum(total_mentions,   na.rm = TRUE),
    raw_synonyms     = paste(sort(unique(raw_synonyms)),       collapse = "; "),
    match_sources_used = paste(sort(unique(match_sources_used)), collapse = "; "),
    .groups = "drop"
  )

# Printing unreviewed variables after final collapsing step and summarizing category counts
remaining_unreviewed <- filter(variables_std_summary, var_category == "unreviewed")
print(remaining_unreviewed, n = Inf)
table(variables_std_summary$var_category)

write.csv(variables_std_summary, "variables_std_summary.csv")

# Checking each var_category for redundancies/misgroupings ----
cat_list <- variables_std_summary %>% 
  group_by(var_category) %>% 
  summarise(vars = list(variable_std), .groups = "drop") %>% 
  deframe()

# Sanity-check:
sapply(cat_list, length)  

# Calculating percentages for each var ----
N_papers <- nrow(tab) 

variable_processed_df <- variables_std_summary %>%
  mutate(
    total_mentions_in_roles = total_mentions,
    perc_exposure = (exposure / N_papers) * 100,
    perc_mediator = (mediator / N_papers) * 100,
    perc_unobserved_bias = (unobserved_bias / N_papers) * 100
  ) %>%
  arrange(desc(total_mentions_in_roles))


# Summarising by variable category ----
category_summary_df <- variable_processed_df %>%
  group_by(var_category) %>%
  summarise(
    sum_exposure_cat = sum(exposure, na.rm = TRUE),
    sum_mediator_cat = sum(mediator, na.rm = TRUE),
    sum_unobserved_bias_cat = sum(unobserved_bias, na.rm = TRUE),
    total_mentions_cat = sum(total_mentions_in_roles, na.rm = TRUE), 
    unique_vars_in_cat = n_distinct(variable_std),
    .groups = 'drop' 
  ) %>%
  arrange(desc(total_mentions_cat))


# Identifying roles of each variable ----
exposure_threshold_count <- 3 # mentioned as exposure in at least 3 papers
mediator_threshold_count <- 2 # mentioned as mediator in at least 2 papers
bias_threshold_count <- 2     # mentioned as bias in at least 2 papers

variables_with_roles_identified <- variable_processed_df %>%
  mutate(
    is_exposure = ifelse(exposure >= exposure_threshold_count, TRUE, FALSE),
    is_mediator = ifelse(mediator >= mediator_threshold_count, TRUE, FALSE),
    is_unobserved_bias = ifelse(unobserved_bias >= bias_threshold_count, TRUE, FALSE),
    num_roles_identified = is_exposure + is_mediator + is_unobserved_bias
  ) %>%
  filter(num_roles_identified > 0) %>% # only variables meeting at least one role threshold
  select(variable_std, var_category, exposure, mediator, unobserved_bias, 
         perc_exposure, perc_mediator, perc_unobserved_bias,
         is_exposure, is_mediator, is_unobserved_bias, num_roles_identified, total_mentions_in_roles) %>%
  arrange(desc(num_roles_identified), desc(total_mentions_in_roles))

# Creating tiered exposures/mediators/unobserved biases ----
thresholds <- c(tier1 = 0.10,   # ≥ 10 % of papers
                tier2 = 0.05,   # ≥ 5 %
                tier3 = 0.01)   # ≥ 1 %

roles <- list(
  exposures            = "exposure",
  mediators            = "mediator",
  unobserved_biases    = "unobserved_bias"
)

for (role_name in names(roles)) {       
  role_col <- roles[[role_name]]     
  
  for (tier in names(thresholds)) { 
    n_thres <- round(N_papers * thresholds[[tier]])
    
    assign(
      paste0(tier, "_", role_name), 
      variables_with_roles_identified %>% 
        filter(.data[[role_col]] >= n_thres) %>% 
        arrange(desc(.data[[role_col]]))
    )
  }
}

# Saving for supplementary materials ----
write.csv(category_summary_df, file.path(getwd(), "Annexes and Supplementary Materials", "Suppl_Table_S1_Category_Role_Summary.csv"))
write.csv(variable_processed_df, file.path(getwd(), "Annexes and Supplementary Materials", "Suppl_Table_S2_Full_Variable_Roster.csv"))
write.csv(variables_with_roles_identified, file.path(getwd(), "Annexes and Supplementary Materials", "Suppl_Data_S3_Threshold_Filtered_Variables.csv"))

# creating a df of vars meeting baseline role threshold ----
n_total_vars <- nrow(variable_processed_df)

# count of vars exceeding baseline role threshold
n_exposure_vars   <- variables_with_roles_identified %>% filter(is_exposure)        %>% nrow()
n_mediator_vars   <- variables_with_roles_identified %>% filter(is_mediator)        %>% nrow()
n_bias_vars       <- variables_with_roles_identified %>% filter(is_unobserved_bias) %>% nrow()
n_union_vars      <- nrow(variables_with_roles_identified)   # ≥ 1 role threshold

# totals of role specific mentions
tot_exposure_all  <- sum(category_summary_df$sum_exposure_cat,        na.rm = TRUE)
tot_mediator_all  <- sum(category_summary_df$sum_mediator_cat,        na.rm = TRUE)
tot_bias_all      <- sum(category_summary_df$sum_unobserved_bias_cat, na.rm = TRUE)
tot_mentions_all  <- tot_exposure_all + tot_mediator_all + tot_bias_all

# role specific mentions among thresholded set
exp_mentions_keep <- variables_with_roles_identified %>%
  filter(is_exposure) %>% summarise(m = sum(exposure)) %>% pull(m)
med_mentions_keep <- variables_with_roles_identified %>%
  filter(is_mediator) %>% summarise(m = sum(mediator)) %>% pull(m)
bias_mentions_keep<- variables_with_roles_identified %>%
  filter(is_unobserved_bias) %>% summarise(m = sum(unobserved_bias)) %>% pull(m)
union_mentions_keep <- variables_with_roles_identified %>%
  summarise(m = sum(exposure + mediator + unobserved_bias)) %>% pull(m)

# df for vars meeting baseline role thresholds
table_annex_df <- tibble(
  `Tier / Role`                    = c("Exposure ≥ 3 papers",
                                       "Mediator ≥ 2 papers",
                                       "Unobs. bias ≥ 2 papers",
                                       "Union (all roles combined)"),
  `Count of Variables Meeting Threshold` = c(n_exposure_vars,
                                             n_mediator_vars,
                                             n_bias_vars,
                                             n_union_vars),
  `% of All 453 Variables`         = round(c(n_exposure_vars,
                                             n_mediator_vars,
                                             n_bias_vars,
                                             n_union_vars) / n_total_vars * 100, 1),
  `% of Role-Specific Mentions Captured`       = round(c(exp_mentions_keep / tot_exposure_all,
                                                         med_mentions_keep / tot_mediator_all,
                                                         bias_mentions_keep / tot_bias_all,
                                                         union_mentions_keep / tot_mentions_all) * 100, 1)
)

# creating a df for role-tier counts ----
get_vars <- function(role_col, n_threshold) {
  variables_with_roles_identified |>
    filter(.data[[role_col]] >= n_threshold) |>
    pull(variable_std) |>
    unique()
}

tier_grid <- expand_grid(
  tier   = names(thresholds),          # "tier1"  "tier2"  "tier3"
  role   = names(roles)                # "exposures" "mediators" "unobserved_biases"
) |>
  mutate(
    role_col   = roles[role],          # "exposure"  "mediator"  "unobserved_bias"
    n_thresh   = round(N_papers * thresholds[tier]),
    var_set    = map2(role_col, n_thresh, get_vars),
    n_vars     = map_int(var_set, length)
  )

union_counts <- tier_grid |>
  group_by(tier) |>
  summarise(
    role   = "union",
    n_vars = length(reduce(var_set, union)),  
    .groups = "drop"
  )

tier_role_counts_df <- bind_rows(tier_grid |> select(tier, role, n_vars),
                                 union_counts) |>
  mutate(
    role  = recode(role,
                   exposures         = "Exposures",
                   mediators         = "Mediators",
                   unobserved_biases = "Unobserved Biases",
                   union             = "Union"),
    tier  = recode(tier,
                   tier1 = "Tier 1",
                   tier2 = "Tier 2",
                   tier3 = "Tier 3")
  ) |>
  pivot_wider(
    names_from  = tier,
    values_from = n_vars,
    values_fill = 0        
  ) |>
  arrange(match(role, c("Exposures", "Mediators",
                        "Unobserved Biases", "Union")))

# TO DO: Create 3 DAGS; tier 1, tier 2, and tier 3; for each, start with exposures, then add in mediators, and unobserved biases. Draw lines according to intuition and knowledge. Edit Sankey diagram to include children with no outcome.

# Tier data construction ----
# tier 1
tier1 <- tier1_exposures %>%
  full_join(tier1_mediators) %>%
  full_join(tier1_unobserved_biases) %>%
  select(
    variable_std, var_category,
    is_exposure, is_mediator, is_unobserved_bias,
    num_roles_identified,
    perc_exposure, perc_mediator, perc_unobserved_bias
  ) %>%
  mutate(tier = 1L)

# tier 2
tier2 <- tier2_exposures %>%
  full_join(tier2_mediators) %>%
  full_join(tier2_unobserved_biases) %>%
  select(
    variable_std, var_category,
    is_exposure, is_mediator, is_unobserved_bias,
    num_roles_identified,
    perc_exposure, perc_mediator, perc_unobserved_bias
  ) %>%
  mutate(
    tier = if_else(variable_std %in% tier1$variable_std, 1L, 2L)
  )

# tier 3
tier3 <- tier3_exposures %>%
  full_join(tier3_mediators) %>%
  full_join(tier3_unobserved_biases) %>%
  select(
    variable_std, var_category,
    is_exposure, is_mediator, is_unobserved_bias,
    num_roles_identified,
    perc_exposure, perc_mediator, perc_unobserved_bias
  ) %>%
  mutate(
    tier = case_when(
      variable_std %in% tier1$variable_std ~ 1L,
      variable_std %in% tier2$variable_std ~ 2L,
      TRUE                                ~ 3L
    )
  )

# saving for supplementary material
write.csv(tier1, "Suppl_Data_S5_Tier_1_Vars.csv", row.names = FALSE)
write.csv(tier2, "Suppl_Data_S6_Tier_2_Vars.csv", row.names = FALSE)
write.csv(tier3, "Suppl_Data_S7_Tier_3_Vars.csv", row.names = FALSE)


# Tier 1 dag ----
library(dplyr)
library(ggplot2)
library(ggdag)
library(stringr)   # for str_to_lower()
library(scales)    # for nicer % legend labels

tier1_dag <- dagify(
  SFI ~ CRP + Procalcitonin + Age + Underlying_Infection + Single_Centre_Study + Small_Sample_Size,
  CRP ~ Age + Underlying_Infection + Measurement_Variability,
  Procalcitonin ~ Age + Underlying_Infection + Measurement_Variability,
  Underlying_Infection ~ Age,
  labels = c(
    SFI                        = "Severe\nFebrile\nIllness",
    CRP                        = "C-reactive\nProtein",
    Procalcitonin              = "Procalcitonin",
    Age                        = "Age",
    Underlying_Infection       = "Underlying\nInfection",
    Single_Centre_Study        = "Single Centre\nStudy",
    Small_Sample_Size          = "Small Sample\nSize",
    Measurement_Variability    = "Measurement\nVariability"
  ),
  exposure = c("CRP", "Procalcitonin"),
  outcome  = "SFI",
  latent   = "Measurement_Variability"
)


tidy_dagitty(tier1_dag)

tier1_dag_graph <- tidy_dagitty(tier1_dag, layout = "auto")

## helper 
std_name <- function(x) stringr::str_to_lower(x) %>% 
  stringr::str_replace_all("centre", "center")

## 1.  lookup table from the literature 
node_meta <- tier1 %>% 
  mutate(
    name_key  = std_name(variable_std),
    perc_max  = pmax(perc_exposure, perc_mediator, perc_unobserved_bias, na.rm = TRUE),
    source    = "Literature"            # provenance flag
  ) %>% 
  select(name_key, var_category, perc_max, source)

tier1_dag_graph <- tidy_dagitty(tier1_dag, layout = "auto") %>% 
  mutate(name_key = std_name(name))

missing_keys <- setdiff(unique(tier1_dag_graph$name_key), node_meta$name_key)

extra_nodes <- tibble(
  name_key     = missing_keys,
  var_category = case_when(
    missing_keys == "sfi"                  ~ "outcome",
    missing_keys == "underlying_infection" ~ "latent_process",
    TRUE                                   ~ "other"
  ),
  perc_max     = 80/100,   # small but visible
  source       = "Manual"
)

node_meta <- dplyr::bind_rows(node_meta, extra_nodes)

tier1_dag_graph <- tier1_dag_graph %>% 
  dplyr::left_join(node_meta, by = "name_key") %>% 
  dplyr::mutate(
    var_category = tidyr::replace_na(var_category, "other"),
    perc_max     = tidyr::replace_na(perc_max, 50),
    source       = tidyr::replace_na(source, "Manual")
  )

tier1_dag_graph.df <- as.data.frame(tier1_dag_graph)

tier1_dag_graph.df$x <- tier1_dag_graph.df$x * 0.9  # squish width
tier1_dag_graph.df$y <- tier1_dag_graph.df$y * 0.9  # squish height
tier1_dag_graph.df$xend <- tier1_dag_graph.df$xend * 0.9 # squish width
tier1_dag_graph.df$yend <- tier1_dag_graph.df$yend * 0.9  # squish height


t1_p <- ggplot(tier1_dag_graph.df) +
  geom_dag_edges_diagonal(
    aes(x, y, xend = xend, yend = yend),
    curvature   = .003,   # increased to .15
    edge_colour = "grey50",
    edge_width  = 0.5,
    edge_alpha  = 1
  ) +
  geom_dag_point(
    aes(x, y, fill = var_category, colour = var_category, shape = source),  # Alpha mapped to perc_max
    size   = 0,  # Uniform size for nodes
    stroke = 14
  ) +
  geom_dag_text(aes(x, y, label = label), size = 3, colour = "black") +
  scale_shape_manual(values = c(Literature = 21, Manual = 22),
                     name   = "Node source") +
  theme_dag_blank() +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    plot.subtitle= element_text(face = "italic"),
    plot.caption = element_text(face = "italic"),
    axis.text    = element_blank(),
    #legend.position = "none"
  ) +
  labs(
    title    = "Tier 1 DAG",
    subtitle = "Variables present in ≥ 10% of papers",
    caption = "Square nodes represent manual addition.",
    x = NULL, y = NULL
  )

t1_p

ggsave(file = "t1.svg", plot = t1_p, width = 7, height = 5)

# Tier 2 dag ----
library(dplyr)
library(ggplot2)
library(ggdag)
library(stringr)   # for str_to_lower()
library(scales)    # for nicer % legend labels

tier2_dag <- dagify(
  SFI ~ CRP + Procalcitonin + Lactate + NLR + Vital_Signs + Age + Underlying_Infection + Study_Bias,
  CRP ~ Age + Underlying_Infection + Measurement_Variability,
  Procalcitonin ~ Age + Underlying_Infection + Measurement_Variability,
  Lactate ~ Age + Underlying_Infection + Measurement_Variability,
  NLR ~ Age + Underlying_Infection + Measurement_Variability,
  Vital_Signs ~ Age + Underlying_Infection,
  Underlying_Infection ~ Age,
  labels = c(
    SFI                    = "Severe\nFebrile\nIllness",
    CRP                    = "C-reactive\nProtein",
    Procalcitonin          = "Procalcitonin",
    Lactate                = "Lactate",
    NLR                    = "NLR",
    Vital_Signs            = "Vital Signs",
    Age                    = "Age",
    Underlying_Infection   = "Underlying\nInfection",
    Measurement_Variability = "Measurement\nVariability",
    Study_Bias             = "Study\nBias"
  ),
  exposure = c("CRP", "Procalcitonin", "Lactate", "NLR", "Vital_Signs"),
  outcome  = "SFI",
  latent   = "Measurement_Variability"
)




tidy_dagitty(tier2_dag)
ggdag(tier2_dag, layout = "auto", text = FALSE, use_labels = "label") 
ggdag(tier2_dag, layout = "circle")
dag_paths(tier2_dag, from = "CRP")
dag_paths(tier2_dag, from = "Procalcitonin")
ggdag_paths(tier2_dag, from = "CRP")
ggdag_paths(tier2_dag, from = "Procalcitonin")
ggdag_adjustment_set(tier2_dag)
ggdag_parents(tier2_dag, "Procalcitonin")
ggdag_parents(tier2_dag, "CRP")

tier2_dag_graph <- tidy_dagitty(tier2_dag, layout = "auto")

## 1.  lookup table from the literature 
node_meta <- tier2 %>% 
  mutate(
    name_key  = std_name(variable_std),
    perc_max  = pmax(perc_exposure, perc_mediator, perc_unobserved_bias, na.rm = TRUE),
    source    = "Literature"            # provenance flag
  ) %>% 
  select(name_key, var_category, perc_max, source)

tier2_dag_graph <- tidy_dagitty(tier2_dag, layout = "auto") %>% 
  mutate(name_key = std_name(name))

missing_keys <- setdiff(unique(tier2_dag_graph$name_key), node_meta$name_key)

extra_nodes <- tibble(
  name_key     = missing_keys,
  var_category = case_when(
    missing_keys == "sfi"                  ~ "outcome",
    missing_keys == "underlying_infection" ~ "latent_process",
    TRUE                                   ~ "other"
  ),
  perc_max     = 80/100,   # small but visible
  source       = "Manual"
)

node_meta <- dplyr::bind_rows(node_meta, extra_nodes)

tier2_dag_graph <- tier2_dag_graph %>% 
  dplyr::left_join(node_meta, by = "name_key") %>% 
  dplyr::mutate(
    var_category = tidyr::replace_na(var_category, "other"),
    perc_max     = tidyr::replace_na(perc_max, 50),
    source       = tidyr::replace_na(source, "Manual")
  )

tier2_dag_graph.df <- as.data.frame(tier2_dag_graph)

tier2_dag_graph.df$x <- tier2_dag_graph.df$x * 0.9  # squish width
tier2_dag_graph.df$y <- tier2_dag_graph.df$y * 0.9  # squish height
tier2_dag_graph.df$xend <- tier2_dag_graph.df$xend * 0.9 # squish width
tier2_dag_graph.df$yend <- tier2_dag_graph.df$yend * 0.9  # squish height


t2_p <- ggplot(tier2_dag_graph.df) +
  geom_dag_edges_diagonal(
    aes(x, y, xend = xend, yend = yend),
    curvature   = .003,   # increased to .15
    edge_colour = "grey50",
    edge_width  = 0.5,
    edge_alpha  = 1
  ) +
  geom_dag_point(
    aes(x, y, fill = var_category, colour = var_category, shape = source),  # Alpha mapped to perc_max
    size   = 0,  # Uniform size for nodes
    stroke = 14
  ) +
  geom_dag_text(aes(x, y, label = label), size = 3, colour = "black") +
  scale_shape_manual(values = c(Literature = 21, Manual = 22),
                     name   = "Node source") +
  theme_dag_blank() +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    plot.subtitle= element_text(face = "italic"),
    plot.caption = element_text(face = "italic"),
    axis.text    = element_blank(),
    #legend.position = "none"
  ) +
  labs(
    title    = "Tier 2 DAG",
    subtitle = "Variables present in ≥ 5% of papers",
    caption = "Square nodes represent manual addition. \nStudy Bias represents four grouped variables.", # single centre study, small sample size, selection bias, retrospective data
    x = NULL, y = NULL
  )

t2_p

ggsave(file = "t2.svg", plot = t2_p, width = 7, height = 5)



# Tier 3 dag ----
library(dplyr)
library(ggplot2)
library(ggdag)
library(stringr)   # for str_to_lower()
library(scales)    # for nicer % legend labels
library(dagitty)

tier3_dag <- dagify(
  SFI ~ CRP + Procalcitonin + Lactate + NLR + ANC + Bilirubin + BNP_NT_ProBNP + IL6 + MPV +
    WBC + NGAL + STREM1 + Platelet_Count + Presepsin + Vital_Signs +
    Shock + Hypoxia + Hypoperfusion + Endothelial_Activation +
    Immune_Cell_Infiltration + Renal_Perfusion + Cytokines + Plasma_PAI1_Level +
    Age + Underlying_Infection + Study_Bias,
  CRP                  ~ Age + Underlying_Infection + Measurement_Variability + IL6, # Refined: IL6 added as a direct cause
  Procalcitonin        ~ Age + Underlying_Infection + Measurement_Variability,
  Lactate              ~ Age + Underlying_Infection + Measurement_Variability,
  CBC                  ~ Age + Underlying_Infection + Measurement_Variability, # Refined: This is the parent node for blood count components
  NLR                  ~ CBC, # Refined: NLR is derived from CBC results
  ANC                  ~ CBC, # Refined: ANC is derived from CBC results
  MPV                  ~ CBC, # Refined: MPV is derived from CBC results
  WBC                  ~ CBC, # Refined: WBC is derived from CBC results
  Platelet_Count       ~ CBC, # Refined: Platelet Count is derived from CBC results
  Bilirubin            ~ Age + Underlying_Infection + Measurement_Variability,
  BNP_NT_ProBNP        ~ Age + Underlying_Infection + Measurement_Variability,
  IL6                  ~ Age + Underlying_Infection + Measurement_Variability,
  NGAL                 ~ Age + Underlying_Infection + Measurement_Variability,
  STREM1               ~ Age + Underlying_Infection + Measurement_Variability,
  Presepsin            ~ Age + Underlying_Infection + Measurement_Variability,
  Vital_Signs          ~ Shock + Hypoxia + Underlying_Infection + Age, # Refined: Causal parents of Vital Signs now specified
  Cytokines            ~ Underlying_Infection + Age,
  Immune_Cell_Infiltration ~ Underlying_Infection + Cytokines,
  Plasma_PAI1_Level    ~ Underlying_Infection + Cytokines + Measurement_Variability,
  Endothelial_Activation ~ Underlying_Infection + Cytokines + Shock,
  Shock                ~ Underlying_Infection + Age,
  Hypoperfusion        ~ Shock + Endothelial_Activation,
  Renal_Perfusion      ~ Hypoperfusion,
  Hypoxia              ~ Shock + Underlying_Infection + Age,
  Underlying_Infection ~ Age + Breastfeeding + Prematurity + Birth_Weight,
  Breastfeeding        ~ Age,
  Birth_Weight         ~ Prematurity,
  labels = c(
    SFI                  = "Severe\nFebrile\nIllness",
    CRP                  = "C-reactive\nProtein",
    Procalcitonin        = "Procalcitonin",
    Lactate              = "Lactate",
    NLR                  = "NLR",
    ANC                  = "Absolute\nNeutrophil\nCount",
    Bilirubin            = "Bilirubin",
    BNP_NT_ProBNP        = "BNP/NT-proBNP",
    CBC                  = "Complete Blood\nCount",
    IL6                  = "IL-6",
    MPV                  = "Mean Platelet\nVolume",
    WBC                  = "White Blood\nCell Count",
    NGAL                 = "NGAL",
    STREM1               = "Plasma sTREM1",
    Platelet_Count       = "Platelet\nCount",
    Presepsin            = "Presepsin",
    Vital_Signs          = "Vital Signs",
    Shock                = "Shock",
    Hypoxia              = "Hypoxia",
    Hypoperfusion        = "Hypoperfusion",
    Endothelial_Activation = "Endothelial\nActivation",
    Immune_Cell_Infiltration = "Immune Cell\nInfiltration",
    Renal_Perfusion      = "Renal\nPerfusion",
    Cytokines            = "Cytokines",
    Plasma_PAI1_Level    = "Plasma PAI-1\nLevel",
    Age                  = "Age",
    Underlying_Infection = "Underlying\nInfection",
    Breastfeeding        = "Breastfeeding",
    Prematurity          = "Prematurity",
    Birth_Weight         = "Birth Weight",
    Measurement_Variability = "Measurement\nVariability",
    Study_Bias           = "Study\nBias"
  ),
  exposure = c("CRP", "Procalcitonin", "Lactate", "NLR", "ANC", "Bilirubin",
               "BNP_NT_ProBNP", "CBC", "IL6", "MPV", "WBC", "NGAL", "STREM1",
               "Platelet_Count", "Presepsin", "Vital_Signs"),
  outcome  = "SFI",
  latent   = "Measurement_Variability"
)


tidy_dagitty(tier3_dag)
ggdag(tier3_dag, layout = "auto", text = TRUE) 

ggdag_adjustment_set(tier3_dag)

tier3_dag_graph <- tidy_dagitty(tier3_dag, layout = "auto")

## 1.  lookup table from the literature 
node_meta <- tier3 %>% 
  mutate(
    name_key  = std_name(variable_std),
    perc_max  = pmax(perc_exposure, perc_mediator, perc_unobserved_bias, na.rm = TRUE),
    source    = "Literature"            # provenance flag
  ) %>% 
  select(name_key, var_category, perc_max, source)

tier3_dag_graph <- tidy_dagitty(tier3_dag, layout = "auto") %>% 
  mutate(name_key = std_name(name))

missing_keys <- setdiff(unique(tier3_dag_graph$name_key), node_meta$name_key)

extra_nodes <- tibble(
  name_key     = missing_keys,
  var_category = case_when(
    missing_keys == "sfi"                  ~ "outcome",
    missing_keys == "underlying_infection" ~ "latent_process",
    TRUE                                   ~ "other"
  ),
  perc_max     = 80/100,  
  source       = "Manual"
)

node_meta <- dplyr::bind_rows(node_meta, extra_nodes)

tier3_dag_graph <- tier3_dag_graph %>% 
  dplyr::left_join(node_meta, by = "name_key") %>% 
  dplyr::mutate(
    var_category = tidyr::replace_na(var_category, "other"),
    perc_max     = tidyr::replace_na(perc_max, 50),
    source       = tidyr::replace_na(source, "Manual")
  )

tier3_dag_graph.df <- as.data.frame(tier3_dag_graph)

# fixing some by hand
tier3_dag_graph.df <- tier3_dag_graph.df %>%
  mutate(source = ifelse(name %in% c("BNP_NT_ProBNP", "CBC", "Cytokines", "Endothelial_Activation", "Plasma_PAI1_Level", "STREM1"),
                         "Literature", source)) 

tier3_dag_graph.df <- tier3_dag_graph.df %>%
  mutate(var_category = ifelse(name %in% c("BNP_NT_ProBNP", "CBC", "Cytokines", "Plasma_PAI1_Level", "STREM1"),
                               "lab_biomarker", var_category)) 

tier3_dag_graph.df <- tier3_dag_graph.df %>%
  mutate(var_category = ifelse(name %in% c("Endothelial_Activation"),
                               "clinical_condition", var_category))


tier3_dag_graph.df$x <- tier3_dag_graph.df$x * 0.9  # squish width
tier3_dag_graph.df$y <- tier3_dag_graph.df$y * 0.9  # squish height
tier3_dag_graph.df$xend <- tier3_dag_graph.df$xend * 0.9 # squish width
tier3_dag_graph.df$yend <- tier3_dag_graph.df$yend * 0.9  # squish height


t3_p <- ggplot(tier3_dag_graph.df) +
  geom_dag_edges_diagonal(
    aes(x, y, xend = xend, yend = yend),
    curvature   = .003,   # increased to .15
    edge_colour = "grey50",
    edge_width  = 0.5,
    edge_alpha  = 1
  ) +
  geom_dag_point(
    aes(x, y, fill = var_category, colour = var_category, shape = source),  # Alpha mapped to perc_max
    size   = 0,  # Uniform size for nodes
    stroke = 0
  ) +
  geom_dag_text(aes(x, y, colour = var_category, label = label), size = 3) +
  scale_shape_manual(values = c(Literature = 21, Manual = 22),
                     name   = "Node source") +
  theme_dag_blank() +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    plot.subtitle= element_text(face = "italic"),
    plot.caption = element_text(face = "italic"),
    axis.text    = element_blank(),
    #legend.position = "none"
  ) +
  labs(
    title    = "Tier 3 DAG",
    subtitle = "Variables present in ≥ 1% of papers",
    caption = "Square nodes represent manual addition. \nStudy Bias represents multiple grouped variables.",
    x = NULL, y = NULL
  )

t3_p

ggsave(file = "t3.svg", plot = t3_p, width = 16, height = 10)

# Minimal adjustment sets ----
library(dagitty)
tier1_adj <- adjustmentSets(tier1_dag) 
tier2_adj <- adjustmentSets(tier2_dag)
tier3_adj <- adjustmentSets(tier3_dag)

# Tier 1 model formulas ----
expo_tier1    <- exposures(tier1_dag)     
outcome <- outcomes(tier1_dag)[1]   
adjset_tier1  <- adjustmentSets(tier1_dag)[[1]] 
tier1_exp_and_adj_vars  <- c(expo_tier1, adjset_tier1)

# joint formula
tier1_formula <- reformulate(tier1_exp_and_adj_vars, response = outcome)

# one per exposure
tier1_form_by_expo <- map(
  expo_tier1,
  function(e) {
    adj_e <- adjustmentSets(tier1_dag, exposure = e)[[1]]
    reformulate(c(e, adj_e), response = outcome)
  }
)
names(tier1_form_by_expo) <- expo_tier1

# Tier 2 model formulas ----
expo_tier2    <- exposures(tier2_dag)     
outcome <- outcomes(tier2_dag)[1]   
adjset_tier2  <- adjustmentSets(tier2_dag)[[1]] 
tier2_exp_and_adj_vars  <- c(expo_tier2, adjset_tier2)

# joint formula
tier2_formula <- reformulate(tier2_exp_and_adj_vars, response = outcome)

# one per exposure
tier2_form_by_expo <- map(
  expo_tier2,
  function(e) {
    adj_e <- adjustmentSets(tier2_dag, exposure = e)[[1]]
    reformulate(c(e, adj_e), response = outcome)
  }
)
names(tier2_form_by_expo) <- expo_tier2

# Tier 3 model formulas ----
expo_tier3    <- exposures(tier3_dag)     
outcome <- outcomes(tier3_dag)[1]   
adjset_tier3  <- adjustmentSets(tier3_dag)[[1]] 
tier3_exp_and_adj_vars  <- c(expo_tier3, adjset_tier3)

# joint formula
tier3_formula <- reformulate(tier3_exp_and_adj_vars, response = outcome)

# individual formulas
build_formula_if_identifiable <- function(exposure,
                                          dag        = tier3_dag,
                                          outcome_var = outcome) {
  adj_sets <- adjustmentSets(dag, exposure = exposure)
  
  if (length(adj_sets) == 0) return(NULL)
  
  sizes <- lengths(adj_sets)             # vector of set sizes
  adj   <- adj_sets[[ which.min(sizes) ]]  # tie = first smallest
  
  reformulate(c(exposure, adj), response = outcome_var)
}

expo_tier3_named <- set_names(expo_tier3)       

tier3_form_by_expo <- expo_tier3_named |>
  map(build_formula_if_identifiable) |>     
  compact()                                

skipped <- setdiff(names(expo_tier3_named), names(form_by_expo))




# Creating files for supplementary material ----
# Function
write_model_summary <- function(dag, adjset, joint_formula, form_by_expo, file) {
  # DAG structure
  dag_lines <- capture.output(dag)
  
  # minimal adjustment set
  adj_line <- paste0("adjset: ", paste(adjset, collapse = " + "))
  
  # joint formula
  joint_line <- paste(deparse(joint_formula), collapse = "")
  
  # exposure-specific formula
  expo_lines <- vapply(
    names(form_by_expo),
    FUN   = function(ex) {
      f <- form_by_expo[[ex]]
      paste0(ex, ": ", paste(deparse(f), collapse = ""))
    },
    FUN.VALUE = character(1)
  )
  
  # all sections
  output_lines <- c(
    "# --- DAG structure ---",
    dag_lines,
    "",
    "# --- Minimal adjustment set ---",
    adj_line,
    "",
    "# --- Joint formula ---",
    joint_line,
    "",
    "# --- Exposure-specific formulas ---",
    expo_lines
  )
  
  writeLines(output_lines, con = file)
}

write_model_summary(tier1_dag, adjset_tier1, tier1_formula, tier1_form_by_expo, "tier1_summary.txt")
write_model_summary(tier2_dag, adjset_tier2, tier2_formula, tier2_form_by_expo, "tier2_summary.txt")
write_model_summary(tier3_dag, adjset_tier3, tier3_formula, tier3_form_by_expo, "tier3_summary.txt")
