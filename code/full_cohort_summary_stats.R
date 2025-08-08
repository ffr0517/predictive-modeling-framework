# Set up ----
training_data <- readRDS("SpotSepsis Data/Derived/RAW_train_df.rds")
testing_data <- readRDS("SpotSepsis Data/Derived/RAW_test_df.rds")
variable_lookup_table <- readRDS("Visualisations/Tables/variable_lookup.rds")

full_dataset <- bind_rows(training_data, testing_data)

library(gtsummary)
library(tidyverse)


# Full summary table creation ----
table1_data <- full_dataset %>%
  select(-label, -outcome.binary, -weight, -weight_bin) %>%
  mutate(
    across(c(sex, bgcombyn, adm.recent, waste, stunt, prior.care,
             travel.time.bin, urti, lrti, diarrhoeal, neuro, auf, ensapro, vomit.all,
             seiz, pfacleth, crt.long, not.alert, danger.sign, syndrome.resp,
             syndrome.nonresp, pneumo, sev.pneumo, infection, parenteral_screen), 
           ~ factor(., levels = c(0, 1), labels = c("No", "Yes")))
  ) %>%
  mutate(
    SFI_5cat = fct_recode(SFI_5cat,
                          "Non-Severe" = "Non-Severe",
                          "Probable Non-Severe" = "Probable Non-Severe",
                          "Probable Severe" = "Probable Severe",
                          "Onset > 24h" = "Onset greater than 24 hours",
                          "Onset < 24h" = "Onset within 24 hours"
    )
  )

# continuous vars
summary_continuous <- table1_data %>%
  select(SFI_5cat, where(is.numeric)) %>%
  pivot_longer(-SFI_5cat, names_to = "Variable", values_to = "value") %>%
  group_by(Variable, SFI_5cat) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    p25 = quantile(value, 0.25, na.rm = TRUE),
    p75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(formatted_stat = sprintf("%.1f (%.1f, %.1f)", median, p25, p75)) %>%
  select(Variable, SFI_5cat, formatted_stat)

# categorical
summary_categorical <- table1_data %>%
  select(SFI_5cat, where(is.factor)) %>%
  pivot_longer(-SFI_5cat, names_to = "Variable", values_to = "level") %>%
  filter(!is.na(level)) %>%
  group_by(Variable, SFI_5cat, level) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Variable, SFI_5cat) %>%
  mutate(total_in_group = sum(n)) %>%
  ungroup() %>%
  mutate(
    p = (n / total_in_group) * 100,
    formatted_stat = sprintf("%d (%.1f%%)", n, p)
  ) %>%
  select(Variable, SFI_5cat, level, formatted_stat)

# P-values for Continuous Variables (Kruskal-Wallis Test)
p_values_continuous <- table1_data %>%
  select(where(is.numeric), SFI_5cat) %>%
  pivot_longer(-SFI_5cat, names_to = "Variable", values_to = "value") %>%
  group_by(Variable) %>%
  summarise(
    p.value = kruskal.test(value ~ SFI_5cat)$p.value,
    .groups = "drop"
  )

# chi sq
p_values_categorical <- table1_data %>%
  select(where(is.factor), -site) %>% 
  map_dfr(~ {
    tryCatch({
      test_result <- chisq.test(table(.x, table1_data$SFI_5cat))
      tibble(p.value = test_result$p.value)
    }, error = function(e) {
      tibble(p.value = NA_real_)
    })
  }, .id = "Variable") %>%
  filter(Variable != "SFI_5cat")

p_value_site_overall <- tibble(
  Variable = "site",
  p.value = chisq.test(table(table1_data$site, table1_data$SFI_5cat))$p.value
)


p_values_sites_individual <- levels(table1_data$site) %>%
  map_dfr(~{
    current_site <- .x
    
    temp_data <- table1_data %>%
      mutate(site_group = if_else(site == current_site, "current", "other"))
    
    p_val <- chisq.test(table(temp_data$site_group, temp_data$SFI_5cat))$p.value
    
    tibble(
      Variable = as.character(current_site),
      p.value = p_val
    )
  })

all_p_values <- bind_rows(
  p_values_continuous, 
  p_values_categorical,
  p_value_site_overall,
  p_values_sites_individual
) %>%
  mutate(`p-value` = format.pval(p.value, digits = 2, eps = 0.001))

overall_continuous <- table1_data %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "value") %>%
  group_by(Variable) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    p25 = quantile(value, 0.25, na.rm = TRUE),
    p75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Overall = sprintf("%.1f (%.1f, %.1f)", median, p25, p75)) %>%
  select(Variable, Overall)

overall_dichotomous <- table1_data %>%
  select(where(is.factor), -SFI_5cat, -site) %>% 
  pivot_longer(everything(), names_to = "Variable", values_to = "level") %>%
  filter(!is.na(level), level == "Yes") %>%
  group_by(Variable) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    total = nrow(table1_data),
    Overall = sprintf("%d (%.1f%%)", n, p = (n / total) * 100)
  ) %>%
  select(Variable, Overall)

overall_site <- table1_data %>%
  count(site) %>%
  mutate(
    total = sum(n),
    Overall = sprintf("%d (%.1f%%)", n, p = (n / total) * 100),
    Variable = as.character(site) 
  ) %>%
  select(Variable, Overall)


# Dichotomous variables:
summary_by_group_dichotomous <- summary_categorical %>%
  filter(Variable != "site", level == "Yes") %>%
  select(Variable, SFI_5cat, formatted_stat)

# Polytomous 'site' variable:
summary_by_group_site <- summary_categorical %>%
  filter(Variable == "site") %>%
  mutate(Variable = as.character(level)) %>% # Use site levels as the variable name
  select(Variable, SFI_5cat, formatted_stat)

# binding all categorical summaries together
summary_by_group_categorical_combined <- bind_rows(
  summary_by_group_dichotomous,
  summary_by_group_site
)

# Pivoting the combined categorical summary to wide format
summary_by_group_wide_categorical <- summary_by_group_categorical_combined %>%
  pivot_wider(
    id_cols = Variable,
    names_from = SFI_5cat,
    values_from = formatted_stat,
    values_fill = "0 (0.0%)"
  )

# Bind with the continuous summary
summary_by_group_wide <- bind_rows(
  summary_continuous %>% pivot_wider(id_cols = Variable, names_from = SFI_5cat, values_from = formatted_stat),
  summary_by_group_wide_categorical
)

# And bind the overall summaries together
summary_overall_wide <- bind_rows(
  overall_continuous,
  overall_dichotomous,
  overall_site
)

# final table ----
final_table <- summary_overall_wide %>%
  left_join(summary_by_group_wide, by = "Variable") %>%
  left_join(all_p_values, by = "Variable") %>%
  left_join(variable_lookup_table, by = c("Variable" = "original_name")) %>%
  mutate(
    Characteristic = coalesce(clean_name, Variable)  
  ) %>%
  mutate(
    Characteristic = recode(Characteristic,
                            id003 = "Site: Rumah Sakit Umum Daerah Wates, Indonesia",
                            kh005 = "Site: Angkor Hospital, Cambodia",
                            la004 = "Site: Salavan Provincial Hospital, Laos",
                            la011 = "Site: Savannakhet Provincial Hospital, Laos",
                            vn009 = "Site: Dong Nai Childrens Hospital, Viet Nam",
                            vn010 = "Site: National Childrens Hospital, Viet Nam",
                            bd006 = "Site: Goyalmara Mother and Child Hospital, Bangladesh",
                            .default = Characteristic
    )
  ) %>%
  select(
    Characteristic,
    Overall,
    `Non-Severe`,
    `Probable Non-Severe`,
    `Probable Severe`,
    `Onset > 24h`,
    `Onset < 24h`,
    `p-value`
  ) %>%
  arrange(match(Characteristic, variable_lookup_table$clean_name))

# total number of patients
overall_n <- nrow(table1_data)

# number of patients in each outcome group
group_counts <- table1_data %>%
  count(SFI_5cat)

# Rename the columns of the final table with the counts
final_table_named <- final_table %>%
  rename(
    !!paste0("Overall (N = ", overall_n, ")") := Overall,
    !!paste0("Non-Severe (N = ", group_counts$n[group_counts$SFI_5cat == "Non-Severe"], ")") := `Non-Severe`,
    !!paste0("Probable Non-Severe (N = ", group_counts$n[group_counts$SFI_5cat == "Probable Non-Severe"], ")") := `Probable Non-Severe`,
    !!paste0("Probable Severe (N = ", group_counts$n[group_counts$SFI_5cat == "Probable Severe"], ")") := `Probable Severe`,
    !!paste0("Onset > 24h (N = ", group_counts$n[group_counts$SFI_5cat == "Onset > 24h"], ")") := `Onset > 24h`,
    !!paste0("Onset < 24h (N = ", group_counts$n[group_counts$SFI_5cat == "Onset < 24h"], ")") := `Onset < 24h`
  )

saveRDS(final_table_named, "Visualisations/Tables/full_summary_stats_table.rds")
write.csv(final_table_named, "Annexes and Supplementary Materials/full_summary_stats_table.csv")


