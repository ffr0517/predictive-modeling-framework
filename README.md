# Cascade Modeling (Code & Derived Artifacts)

**No raw data is included.** This repo contains code, trained model artifacts (where permitted), and derived tables used in a project on cascade modeling.

## Structure
- `code/` — R scripts to build, evaluate, and analyze models; plus a simple Shiny app (`cascade_app/app.R`).
- `models/` — Trained artifacts and feature-importance outputs (no identifiers).
- `outputs/tables/` — Derived/supplementary tables (no raw or identifiable data).
- `data-templates/` — Input template to reproduce workflows with your own data.
- `docs/` — Protocols and PDFs (if included), redacted for anonymity.

## Reproducibility (no data)
- Provide data matching `data-templates/input_data_template.csv` columns.
- Update scripts to use **relative paths**. Example:
  ```r
  # install.packages("here")
  library(here)
  template <- read.csv(here("data-templates", "input_data_template.csv"))
