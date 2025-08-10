# ===================================================================
# Cascade Model: Interactive Shiny App (v8 - Documentation Tab)
# ===================================================================

# --- 0. Preamble: Load Libraries and Helper Functions ---
# (This section is unchanged)
library(shiny)
library(tidyverse)
library(tidymodels)
library(hardhat)
library(DT)
library(markdown)

max_prob_class <- function(df, levels) {
  max_col_indices <- max.col(df, ties.method = "first")
  class_names_with_prefix <- colnames(df)
  clean_class_names <- gsub(".pred_", "", class_names_with_prefix, fixed = TRUE)
  factor(clean_class_names[max_col_indices], levels = levels)
}

align_data_to_blueprint <- function(new_data, blueprint) {
  # (Function code is unchanged)
  message("Aligning data to model's blueprint...")
  for (col_name in names(blueprint)) {
    if (!col_name %in% names(new_data)) {
      message(paste(" -> NOTE: Blueprint column '", col_name, "' not found. Adding as NA."))
      expected_class <- class(blueprint[[col_name]])
      if ("factor" %in% expected_class) {
        new_data[[col_name]] <- factor(NA, levels = levels(blueprint[[col_name]]))
      } else if ("integer" %in% expected_class) {
        new_data[[col_name]] <- NA_integer_
      } else if ("numeric" %in% expected_class || "double" %in% expected_class) {
        new_data[[col_name]] <- NA_real_
      } else if ("character" %in% expected_class) {
        new_data[[col_name]] <- NA_character_
      }
      next
    }
    expected_class <- class(blueprint[[col_name]])
    if ("factor" %in% expected_class) {
      expected_levels <- levels(blueprint[[col_name]])
      new_data <- new_data %>% mutate(!!col_name := factor(!!sym(col_name), levels = expected_levels))
    } else if ("integer" %in% expected_class) {
      new_data <- new_data %>% mutate(!!col_name := as.integer(!!sym(col_name)))
    } else if ("numeric" %in% expected_class || "double" %in% expected_class) {
      new_data <- new_data %>% mutate(!!col_name := as.numeric(!!sym(col_name)))
    } else if ("character" %in% expected_class) {
      new_data <- new_data %>% mutate(!!col_name := as.character(!!sym(col_name)))
    }
  }
  message("-> Data alignment complete.")
  return(new_data)
}


# --- Load the Baked-In Model and Define Required Columns ---
# (This section is unchanged)
message("Attempting to load the baked-in model...")
model_path <- "Mock App/full_cascade_model.rds"
if (!file.exists(model_path)) {
  stop(paste("FATAL ERROR: The model file is missing. The app requires '", model_path, "' to be present relative to app.R.", sep=""))
}
full_model <- readRDS(model_path)
model_blueprint <- hardhat::extract_mold(full_model$stage_1_fit)$blueprint$ptypes$predictors
ALL_REQUIRED_COLS <- c("patient_id", names(model_blueprint))
message("-> Model loaded successfully.")

# ===================================================================
# 1. DEFINE THE USER INTERFACE (UI)
# ===================================================================
ui <- navbarPage(
  "Cascade Model Simulator", # Title for the whole app
  
  # --- TAB 1: The Main Application ---
  tabPanel("Run Analysis",
           # We place the entire previous UI inside this first tabPanel
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h4("1. Upload Patient Data"),
               fileInput("clinical_file", "Upload Clinical-Only Data (.csv)", accept = ".csv"),
               fileInput("full_data_file", "Upload Full Lab Data (.csv)", accept = ".csv"),
               downloadButton("download_template", "Download Template (.csv)"),
               hr(),
               h4("2. Run Analysis"),
               actionButton("run_s1", "Run Stage 1 Analysis", class = "btn-primary"),
               uiOutput("s2_button_ui")
             ),
             mainPanel(
               width = 9,
               h3("Stage 1 Results"),
               p("Initial assessment based on clinical data only. Patients recommended for handoff will have their full lab data processed in Stage 2."),
               DT::dataTableOutput("s1_table"),
               uiOutput("final_results_ui")
             )
           ),
           hr(), 
           fluidRow(
             column(width = 12,
                    h4("Analysis Log"),
                    verbatimTextOutput("log_output")
             )
           )
  ),
  
  # --- TAB 2: The User Guide from the Markdown file ---
  tabPanel("User Guide",
           fluidRow(
             column(10, offset = 1, # Center the content for readability
                    # This function reads and renders the markdown file
                    includeMarkdown("guide.md")
             )
           )
  )
)

# ===================================================================
# 2. DEFINE THE SERVER LOGIC
# ===================================================================
# NO CHANGES ARE NEEDED HERE. The server logic is identical to v7.
server <- function(input, output, session) {
  rv <- reactiveValues(log_text = "Welcome! The model is loaded. Please upload patient data and run Stage 1.", s1_results = NULL, patients_for_handoff = NULL)
  output$log_output <- renderText({ rv$log_text })
  
  show_error_modal <- function(message_html, technical_error = NULL) {
    showModal(modalDialog(
      title = tagList(icon("triangle-exclamation"), "Analysis Error"),
      message_html,
      footer = tagList(
        if (!is.null(technical_error)) {
          div(style = "text-align: left; font-size: 0.8em; color: grey;", "Technical details: ", technical_error)
        },
        modalButton("Dismiss")
      ),
      easyClose = TRUE
    ))
  }
  
  output$download_template <- downloadHandler(
    filename = function() { "input_data_template.csv" },
    content = function(file) {
      template_df <- data.frame(matrix(ncol = length(ALL_REQUIRED_COLS), nrow = 0))
      colnames(template_df) <- ALL_REQUIRED_COLS
      write.csv(template_df, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$run_s1, {
    req(input$clinical_file)
    tryCatch({
      clinical_df <- read.csv(input$clinical_file$datapath)
      uploaded_cols <- colnames(clinical_df)
      if (!"patient_id" %in% uploaded_cols) {
        show_error_modal(p("The uploaded clinical data file is missing the essential 'patient_id' column."))
        return()
      }
      if (length(intersect(names(model_blueprint), uploaded_cols)) == 0) {
        show_error_modal(p("The uploaded clinical data file contains no predictor columns recognized by the model. Please check the file against the template."))
        return()
      }
      withProgress(message = 'Running Stage 1...', value = 0, {
        rv$log_text <- "STEP 1: Data validation passed. Loading data...\n"
        incProgress(0.1)
        rv$log_text <- paste(rv$log_text, "-> Loaded data for", nrow(clinical_df), "new patients.\n\n")
        incProgress(0.3, detail = "Aligning data...")
        log_capture <- capture.output({
          s1_input_aligned <- align_data_to_blueprint(new_data = clinical_df, blueprint = model_blueprint)
          s1_preds <- predict(full_model$stage_1_fit, s1_input_aligned, type = "prob") %>%
            setNames(gsub("[- ]", "_", colnames(.))) %>%
            cal_apply(full_model$calibrator) %>%
            mutate(.pred_class_s1 = max_prob_class(., levels = full_model$s1_levels))
          handoff_mask <- !(s1_preds[[".pred_Non_Severe"]] > full_model$optimal_thresholds$thresh_non_severe | s1_preds[[".pred_Onset_within_24_hours"]] > full_model$optimal_thresholds$thresh_onset_24h)
          s1_results <- tibble(patient_id = s1_input_aligned$patient_id, s1_prediction = s1_preds$.pred_class_s1, handoff_recommended = handoff_mask)
          patients_for_handoff <- s1_results %>% filter(handoff_recommended)
          rv$s1_results <- s1_results
          rv$patients_for_handoff <- patients_for_handoff
        }, type = "message")
        rv$log_text <- paste(rv$log_text, "STEP 2: Running Stage 1 Analysis...\n", paste(log_capture, collapse="\n"), "\n")
        incProgress(0.8, detail = "Generating table...")
        output$s1_table <- DT::renderDataTable({ datatable(rv$s1_results, options = list(pageLength = 10), rownames = FALSE) %>% formatStyle('handoff_recommended', target = 'row', backgroundColor = styleEqual(c(TRUE, FALSE), c('#FFF0F0', '#F0FFF0')))})
        rv$log_text <- paste(rv$log_text, "-> Stage 1 complete. Review results.\n")
        if (nrow(rv$patients_for_handoff) > 0) {
          rv$log_text <- paste(rv$log_text, "-> Patients flagged for Stage 2:", paste(rv$patients_for_handoff$patient_id, collapse = ", "), "\n")
        } else {
          rv$log_text <- paste(rv$log_text, "-> No patients flagged for Stage 2.\n")
        }
        incProgress(1)
      })
    }, error = function(e) {
      rv$log_text <- paste(rv$log_text, "\nERROR: Analysis failed. See pop-up for details.")
      error_message_html <- tagList(p(strong("The analysis failed due to an issue with the uploaded data.")),p("This often happens because of data quality problems, such as:"),tags$ul(tags$li("A column that should be numeric contains text (e.g., 'unknown', 'N/A', or '?')."),tags$li("A categorical column contains a value the model was not trained on."),tags$li("The .csv file is malformed (e.g., 'incomplete final line').")),p("Please carefully check your data against the template, fix any issues, and try again."))
      show_error_modal(error_message_html, technical_error = e$message)
    })
  })
  
  output$s2_button_ui <- renderUI({
    req(rv$s1_results)
    if (nrow(rv$patients_for_handoff) > 0) { actionButton("run_s2", "Run Stage 2 for Flagged Patients", class = "btn-success") } 
    else { p(strong("No patients required handoff. Process complete.")) }
  })
  
  observeEvent(input$run_s2, {
    req(rv$s1_results, rv$patients_for_handoff, input$full_data_file)
    tryCatch({
      full_df <- read.csv(input$full_data_file$datapath)
      uploaded_cols_s2 <- colnames(full_df)
      missing_cols_s2 <- setdiff(ALL_REQUIRED_COLS, uploaded_cols_s2)
      if (length(missing_cols_s2) > 0) {
        if (length(missing_cols_s2) > (length(ALL_REQUIRED_COLS) - 3)) {
          error_msg <- p("The uploaded full lab data file appears incorrect.")
        } else {
          error_msg <- tagList(p("The uploaded full lab data file is missing critical columns:"), tags$ul(lapply(missing_cols_s2, tags$li)))
        }
        show_error_modal(error_msg)
        return()
      }
      withProgress(message = 'Running Stage 2...', value = 0, {
        rv$log_text <- paste(rv$log_text, "\nSTEP 4: Data validation passed. Looking up full data and running Stage 2...\n")
        incProgress(0.2, detail = "Loading full data...")
        log_capture <- capture.output({
          if (nrow(rv$patients_for_handoff) > 0) {
            s2_input_raw <- full_df %>% filter(patient_id %in% rv$patients_for_handoff$patient_id)
            s2_input_aligned <- align_data_to_blueprint(new_data = s2_input_raw, blueprint = model_blueprint)
            base_preds_for_s2 <- map_dfc(full_model$s2_specialist_fits, ~{predict(.x, new_data = s2_input_aligned, type = "prob")})
            s2_final_preds <- predict(full_model$s2_meta_learner_fit, new_data = base_preds_for_s2, type = "class")
            s2_results <- tibble(patient_id = s2_input_aligned$patient_id, s2_prediction = s2_final_preds$.pred_class)
          }
        }, type = "message")
        rv$log_text <- paste(rv$log_text, paste(log_capture, collapse="\n"), "\n")
        incProgress(0.7, detail = "Compiling final report...")
        final_results <- rv$s1_results %>%
          left_join(s2_results, by = "patient_id") %>%
          mutate(
            safeguard_activated = handoff_recommended & (as.character(s1_prediction) %in% full_model$severe_classes_std),
            final_prediction = case_when(safeguard_activated ~ s1_prediction, !is.na(s2_prediction) ~ s2_prediction, TRUE ~ s1_prediction),
            prediction_source = case_when(safeguard_activated ~ "Stage 1 (Safeguard)", !is.na(s2_prediction) ~ "Stage 2", TRUE ~ "Stage 1 (Confident)")
          ) %>%
          select(patient_id, final_prediction, prediction_source, s1_prediction, handoff_recommended)
        output$final_results_ui <- renderUI({ tagList(hr(), h3("Final Compiled Report"), p("This table combines predictions."), DT::dataTableOutput("final_table"))})
        output$final_table <- DT::renderDataTable({ datatable(final_results, options = list(pageLength = 10), rownames = FALSE) %>% formatStyle('prediction_source', backgroundColor = styleEqual(c("Stage 1 (Confident)", "Stage 2", "Stage 1 (Safeguard)"), c('#E8F5E9', '#FFF3E0', '#FFEBEE')))})
        rv$log_text <- paste(rv$log_text, "\nSTEP 5: Final results compiled. Process complete.\n")
        incProgress(1)
      })
    }, error = function(e) {
      rv$log_text <- paste(rv$log_text, "\nERROR: Analysis failed. See pop-up for details.")
      error_message_html <- tagList(p(strong("The Stage 2 analysis failed due to an issue with the uploaded data.")), p("Please check the 'Full Lab Data' file for data quality issues and try again."))
      show_error_modal(error_message_html, technical_error = e$message)
    })
  })
}

# ===================================================================
# 3. RUN THE APPLICATION
# ===================================================================
shinyApp(ui = ui, server = server)
