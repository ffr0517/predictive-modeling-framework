# Set up ----
library(tidyverse)
library(janitor) 
library(glue)
library(dagitty)
library(ggdag)
library(readr)
library(purrr)
library(fuzzyjoin)

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

# reading extraction table data ----
tab  <- read_csv("final_extraction_table.csv") %>% 
  clean_names() %>%                     
  mutate(paper_id = row_number())      

split_vars <- function(x) {
  x %>% str_split(";\\s*") %>%               
    map(~ .x[.x != "" & .x != "No"])          
}

# long table of all raw variables and their given roles (according to the extraction table) ----
var_long <- tab %>%                                      
  transmute(
    paper_id = x1,
    exposures          = split_vars(exposures),
    mediators          = split_vars(mediators),
    unmeasured_biases  = split_vars(unmeasured_biases),
    outcomes           = split_vars(outcomes)           
  ) %>%
  pivot_longer(-paper_id, names_to = "role", values_to = "vars") %>%
  unnest(vars) %>%                                      
  mutate(
    vars = str_trim(vars),                               
    role = recode(role,
                  exposures          = "exposure",
                  mediators          = "mediator",
                  unmeasured_biases  = "unobserved_bias",
                  outcomes           = "outcome")
  ) %>%
  select(paper_id, variable_raw = vars, role) %>%
  arrange(paper_id, role, variable_raw)

write_csv(var_long, "sys_review_variables_and_roles_all.csv")

# summary table generation ----
variables_summary <- var_long %>%
  count(variable_raw, role, name = "n_papers") %>%
  pivot_wider(names_from = role,
              values_from = n_papers,
              values_fill = 0) %>%
  arrange(desc(exposure + mediator + unobserved_bias))     

vars <- variables_summary %>%
  mutate(raw_clean = variable_raw %>%
           str_to_lower() %>%
           str_replace_all("[^a-z0-9]+", " ") %>%
           str_squish())

# grouping variables 
lookup <- tribble(
  # ── clinical / vital signs ─────────────────────────────────
  ~pattern,                               ~variable_std,  ~var_category,
  "vital(s?) sign(s?)|vitals|clinical sign(s?)", "vital_signs", "clinical_sign",
  
  # ── lab biomarkers ─────────────────────────────────────────
  "crp|c reactive protein",               "crp",          "lab_biomarker",
  "pct|procalcitonin",                    "pct",          "lab_biomarker",
  "wbc|white blood cell",                 "wbc",          "lab_biomarker",
  "nlr|neutrophil.*lymphocyte.*ratio",    "nlr",          "lab_biomarker",
  "plt|platelet( count| indices)?",       "plt",          "lab_biomarker",
  "lactate( clearance| blood)?",          "lactate",      "lab_biomarker",
  "adamts13",                             "adamts13",     "lab_biomarker",
  "il ?6( testing)?",                     "il6",          "lab_biomarker",
  
  # ── demographic / host factors ─────────────────────────────
  "age",                                  "age",          "demographic",
  "sex",                                  "sex",          "demographic",
  "birth weight|birthweight",             "birth_weight", "demographic",
  "gestational age",                      "gest_age",     "demographic",
  "prematurity|preterm",                  "prematurity",  "demographic",
  
  # ── nutrition / socio-economic ─────────────────────────────
  "minimum dietary diversity",            "min_diet_div", "nutrition",
  "minimum acceptable diet",              "mad",          "nutrition",
  "yingyangbao",                          "yingyangbao",  "nutrition",
  "wealth index|household crowding",      "socio_econ",   "socioeconomic",
  
  # ── study-design biases ────────────────────────────────────
  "single.?centre|single.?site|single hospital", "bias_single_centre", "design_bias",
  "small sample|small cohort|limited sample size", "bias_small_sample", "design_bias",
  "retrospective|retrospective design",   "bias_retrospective", "design_bias",
  "selection bias|convenience sampling",  "bias_selection", "design_bias",
  
  # ── measurement / assay biases ─────────────────────────────
  "measurement error|measurement bias|device measurement variation", "bias_measurement", "measurement_bias",
  "batch effects",                        "bias_batch",  "measurement_bias",
  "assay variability|assay timing|assay cost|assay complexity", "bias_assay", "measurement_bias",
  "timing variability|biomarker timing",  "bias_timing", "measurement_bias",
  
  # ── generic confounding bucket ─────────────────────────────
  "residual confounding|confounding",     "bias_confounding", "confounding_bias"
)

vars2 <- vars %>%
  # left join on REGEX patterns
  fuzzyjoin::regex_left_join(lookup, 
                             by = c(raw_clean = "pattern")) %>%
  mutate(
    variable_std  = coalesce(variable_std, raw_clean),   
    var_category  = coalesce(var_category, "unreviewed")
  ) %>%
  select(variable_raw, variable_std, var_category, everything())

variables_std_summary <- vars2 %>%
  group_by(variable_std, var_category) %>%           # collapse by harmonised name
  summarise(
    across(                                          # sum the role counts
      c(exposure, mediator, unobserved_bias, outcome),
      ~ sum(.x, na.rm = TRUE)
    ),
    raw_synonyms = paste(unique(variable_raw),       # collect the spellings
                         collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(exposure + mediator + unobserved_bias + outcome))

# table of variables to group (both standardised name and category) ----
vars_to_review <- variables_std_summary %>%
  filter(var_category != "unreviewed")