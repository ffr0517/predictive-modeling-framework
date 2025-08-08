# Set up ---- 
# packages 
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)

# extraction table
extrac_tab <- read.csv("extraction_table_complete.csv")

# Renaming cols ----
colnames(extrac_tab)[1:13] <- c(
  "title", "location", "design", "size", "modelling_used", "internal_validation", "ROB", 
  "exposures", "outcomes", "temporal_ordering", "mediators", "unmeasured_biases", "causal_framework"
)

# write.csv(extrac_tab, "final_extraction_table.csv")

# Modelling approaches in the literature ----
ungrouped_modelling_counts <- extrac_tab %>%
  select(modelling_used) %>%
  filter(!is.na(modelling_used)) %>%
  mutate(modelling_used = str_split(modelling_used, ";\\s*")) %>%
  unnest(modelling_used) %>%
  mutate(modelling_used = str_trim(modelling_used)) %>%
  count(modelling_used, sort = TRUE)

# Grouping modelling approaches 
map_modelling <- function(x) {
  case_when(
    str_detect(x, regex("logistic regression",         ignore_case=TRUE)) ~ "Logistic regression",
    str_detect(x, regex("^binary logistic",            ignore_case=TRUE)) ~ "Logistic regression",
    str_detect(x, regex("^multivariable logistic|^multivariate logistic|multiple logistic|conditional logistic|multilevel logistic|stepwise logistic", ignore_case=TRUE)) ~ "Logistic regression",
    
    str_detect(x, regex("^cox",                        ignore_case=TRUE)) ~ "Cox PH model",
    str_detect(x, regex("poisson regression",          ignore_case=TRUE)) ~ "Poisson regression",
    str_detect(x, regex("generalized linear model|^glm|linear regression|mixed-effects? models?", ignore_case=TRUE)) ~ "General/Mixed linear model",
    
    str_detect(x, regex("random forest",               ignore_case=TRUE)) ~ "Random forest",
    str_detect(x, regex("\\bsvm|support vector machine", ignore_case=TRUE)) ~ "SVM",
    
    str_detect(x, regex("classification tree|decision tree|cart", ignore_case=TRUE)) ~ "Decision tree / CART",
    
    str_detect(x, regex("extreme gradient boosting|xgboost", ignore_case=TRUE)) ~ "Gradient-boosted trees",
    str_detect(x, regex("lasso|elasticnet",            ignore_case=TRUE)) ~ "Penalised regression",
    str_detect(x, regex("\\b(pca|fpca)\\b",            ignore_case=TRUE)) ~ "PCA",
    str_detect(x, regex("spline|lowess",               ignore_case=TRUE)) ~ "Spline regression",
    str_detect(x, regex("opls-?da",                    ignore_case=TRUE)) ~ "OPLS-DA",
    str_detect(x, regex("cluster|clustering|consensus clustering|network clustering", ignore_case=TRUE)) ~ "Clustering",
    str_detect(x, regex("wgcna|coexpression network|hubgene|cibersort", ignore_case=TRUE)) ~ "Gene co-expression/network",
    str_detect(x, regex("ssgsea|gsva",                 ignore_case=TRUE)) ~ "Gene-set enrichment",
    str_detect(x, regex("differential expression|\\bdeg\\b", ignore_case=TRUE)) ~ "Differential expression",
    str_detect(x, regex("twosample mr|two-sample mr",  ignore_case=TRUE)) ~ "Two-sample MR",
    
    str_detect(x, regex("randomeffects metaanalysis", ignore_case=TRUE)) ~ "Random-effects meta-analysis",
    str_detect(x, regex("^metaanalysis$", ignore_case=TRUE))            ~ "Random-effects meta-analysis",
    
    str_detect(x, regex("oaxaca decomposition", ignore_case=TRUE))       ~ "Oaxaca decomposition",
    
    str_detect(x, regex("3[- ]?item scoring model|combined EOS risk|nomogram|ISTH DIC scoring|predictive score", ignore_case=TRUE)) 
    ~ "Risk-scoring models",
    
    TRUE ~ "Other (non-modelling) approaches used"
  )
}

# Adding new grouped modelling approach column
extrac_tab <- extrac_tab %>%
  mutate(
    modelling_type = modelling_used %>% 
      str_split(";\\s*") %>%  
      map_chr(~ .x %>% 
                map_chr(map_modelling) %>% 
                unique() %>% 
                paste(collapse = "; "))
  )

# Adding paper_id and obj for total papers
extrac_tab <- extrac_tab %>% mutate(paper_id = row_number()) 
total_papers <- 155

# Table of modelling counts/percentages
grouped_modelling_counts <- extrac_tab %>%
  filter(!is.na(modelling_type)) %>%
  separate_rows(modelling_type, sep = ";\\s*") %>%
  mutate(modelling_type = str_trim(modelling_type)) %>%
  distinct(paper_id, modelling_type) %>%
  count(modelling_type, name = "n", sort = TRUE) %>%
  mutate(pct = n / total_papers * 100)

# List of modelling approaches that feed each type
modelling_raw_by_type <- extrac_tab %>%
  filter(!is.na(modelling_used)) %>%
  separate_rows(modelling_used, sep = ";\\s*") %>%
  mutate(modelling_used = str_trim(modelling_used)) %>%
  mutate(modelling_type = map_modelling(modelling_used)) %>%
  group_by(modelling_type) %>%
  summarize(
    raw_values = list(sort(unique(modelling_used))),
    .groups = "drop"
  ) %>%
  deframe()

# Creating combined df for supplementary showcase
grouped_modelling_counts_with_raw <- grouped_modelling_counts %>%
  mutate(
    raw_values = modelling_raw_by_type[modelling_type],
    raw_values = sapply(raw_values, function(x) paste(x, collapse = ", "))
  )

# Write to CSV
# write.csv(grouped_modelling_counts_with_raw, "Suppl_Data_S4_Modelling_Approaches.csv", row.names = FALSE)

