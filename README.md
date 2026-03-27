# Spot Sepsis Cascade

Code and derived artefacts for a two-stage machine learning cascade designed to support early severe-febrile-illness workflows in paediatric research settings.

**No raw data is included.** This repository contains modelling code, trained artefacts where permitted, derived tables, templates, and a Shiny-based interactive mockup application.

## What is here

- `code/` — R scripts for model building, evaluation, analysis, and reporting
- `models/` — trained artefacts and feature-importance outputs without identifiers
- `outputs/tables/` — derived and supplementary outputs without raw data
- `data-templates/` — input templates for reproducing workflows with your own compatible data
- `docs/` — redacted protocols, notes, and supporting PDFs
- `cascade_app/` — Shiny app demonstrating the two-stage cascade workflow

## Reproducibility

To run the modelling workflows locally:

1. Provide data matching the schema in `data-templates/input_data_template.csv`.
2. Update path handling to use relative paths. A typical pattern is:

```r
library(here)

template <- read.csv(here("data-templates", "input_data_template.csv"))
```

## Cascade App Guide

The `cascade_app/` folder contains a Shiny application that demonstrates the two-stage workflow.

### Purpose

- **Stage 1** uses minimal clinical data to identify lower- and higher-risk patients.
- **Stage 2** applies a more resource-intensive model to the Stage 1 subset that should be escalated.

### Typical flow

1. Upload a clinical-only `.csv` with `patient_id` and the required predictors.
2. Upload a full-data `.csv` with clinical and laboratory inputs for the same patients.
3. Run Stage 1 to identify rows where `handoff_recommended = TRUE`.
4. Run Stage 2 on the flagged subset to produce a compiled final report.

### Final report interpretation

- **Stage 1 (Confident)** — Stage 1 was sufficiently confident to stand on its own.
- **Stage 2** — the Stage 2 model produced the final prediction after escalation.
- **Stage 1 (Safeguard)** — Stage 1 flagged a potentially severe case and that safeguard outcome is retained.
