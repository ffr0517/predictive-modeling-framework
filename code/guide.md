# Cascade Model App - User Guide

***

### 1. Purpose of the Application

This application provides an interactive interface for the two-stage cascade model. It allows a user to get an initial risk assessment for a batch of patients using minimal clinical data (Stage 1) and then run a more resource-intensive, accurate model on a smaller, flagged subset of high-risk patients (Stage 2).

### 2. How to Use the App

The workflow is divided into two main stages, controlled by the buttons in the sidebar.

1.  **Upload Data**:
    * **Clinical-Only Data**: Upload a `.csv` file containing the initial data for all patients. This file must contain at least the `patient_id` column and some of the clinical predictors.
    * **Full Lab Data**: Upload a `.csv` file containing the complete dataset (including all lab results) for the same group of patients. This will be accessed during Stage 2.
    * Click **Download Template** to get a `.csv` file with all required column headers.

2.  **Run Stage 1 Analysis**:
    * Click this button after uploading the clinical data.
    * The app will run the first model and display the results in the "Stage 1 Results" table.
    * Rows highlighted in pink (`handoff_recommended = TRUE`) are flagged for the more detailed Stage 2 analysis.

3.  **Run Stage 2 Analysis**:
    * If any patients are flagged, this button will become active.
    * Click it to run the second-stage model on the full lab data for only the flagged patients.
    * A "Final Compiled Report" will appear, showing the definitive prediction for each patient and its source.

### 3. Interpreting the Final Report

The `prediction_source` column tells you how the final result was determined:

* **Stage 1 (Confident)**: The initial model was highly confident in its prediction, and no further analysis was needed.
* **Stage 2**: The patient was flagged in Stage 1, and this is the refined prediction from the second model using full data.
* **Stage 1 (Safeguard)**: A special case where a patient was flagged as potentially severe in Stage 1, and that prediction is maintained as a safety measure regardless of the Stage 2 outcome.