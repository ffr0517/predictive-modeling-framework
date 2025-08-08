# Two-Stage Machine Learning Framework to Support Early Detection of Severe Febrile Illness in Young Children (Code & Derived Artifacts)

**No raw data is included.**  
This repository contains code, trained model artifacts (where permitted), derived tables, and an interactive mockup application for a two-stage cascade machine learning model to support the diagnosis of severe febrile illness.

## Structure
- `code/` — R scripts to build, evaluate, and analyse models; plus utilities for visualisation and reporting.
- `models/` — Trained artifacts and feature-importance outputs (no identifiers).
- `outputs/tables/` — Derived/supplementary tables (no raw or identifiable data).
- `data-templates/` — Input template to reproduce workflows with your own data.
- `docs/` — Protocols and PDFs (if included), redacted for anonymity.
- `cascade_app/` — Shiny app providing an interactive interface to the cascade model (see **App User Guide** below).

## Reproducibility (no data)
To run the modelling workflows locally:
1. Provide data matching `data-templates/input_data_template.csv` columns.
2. Update scripts to use **relative paths**. Example:
   ```r
   # install.packages("here")
   library(here)
   template <- read.csv(here("data-templates", "input_data_template.csv"))

# Cascade Model App – User Guide

### 1. Purpose of the Application
The `cascade_app` folder contains a Shiny application that demonstrates the two-stage cascade model. It enables:
- **Stage 1**: Rapid triage using only minimal clinical data to identify low- and high-risk patients.
- **Stage 2**: More resource-intensive, higher-accuracy modelling applied only to the Stage 1 high-risk subset.

### 2. How to Use the App
**Upload Data**  
- *Clinical-Only Data*: `.csv` file containing `patient_id` and clinical predictors.  
- *Full Lab Data*: `.csv` file containing clinical + laboratory data for the same patients.  
- Use **Download Template** to obtain a `.csv` with all required column headers.

**Run Stage 1 Analysis**  
- Click after uploading the clinical data.  
- Results appear in “Stage 1 Results” with rows flagged (`handoff_recommended = TRUE`) for Stage 2.

**Run Stage 2 Analysis**  
- Available only if patients are flagged in Stage 1.  
- Runs the second model on the flagged subset using full data.  
- Produces a “Final Compiled Report” with definitive predictions and their source.

### 3. Interpreting the Final Report
`prediction_source` indicates how the final result was determined:
- **Stage 1 (Confident)** – Stage 1 model highly confident; no further analysis needed.  
- **Stage 2** – Refined prediction from the second model using full data.  
- **Stage 1 (Safeguard)** – Patient flagged as potentially severe in Stage 1; this prediction is retained regardless of Stage 2 outcome.

