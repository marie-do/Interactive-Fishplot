
# Initialize available drugs (can be extended later with dynamic addition)
rv$available_drugs <- c(
  "Midostaurin",
  "Venetoclax",
  "Cytarabine",
  "Azacitidine"
)

# Drug observer to update drug selection choices
observe({
  if (!is.null(rv$available_drugs) && length(rv$available_drugs) > 0) {
    updateSelectInput(
      session,
      "selected_drug",
      choices = rv$available_drugs,
      selected = rv$available_drugs[1]
    )
  }
})

# Dynamic UI for drug effect inputs, which updates based on the selected drug event and impact mode (global vs per-mutation).
output$drug_numeric_inputs <- renderUI({
  
  req(rv$clones_df, input$patient, input$selected_drug_col, input$drug_mode)
  
  mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
  node_ids <- rownames(mat)
  
  key <- paste(
    input$patient,
    input$selected_drug,
    input$selected_drug_col,
    sep = "_"
  )
  
  effect <- rv$drug_effects[[key]]
  
  if (is.null(effect)) {
    
    effect <- list(
      mode = "global",
      global_value = 1,
      per_mutation = setNames(rep(1, length(node_ids)), node_ids)
    )
    
    rv$drug_effects[[key]] <- effect
  }
  
  if (input$drug_mode == "global") { # If global mode is selected, display a single numeric input to control the overall drug impact on all mutations simultaneously.
    numericInput(
      "drug_global",
      "Global drug impact",
      value = effect$global_value,
      min = 0,
      max = 2,
      step = 0.05
    )
    
  } else {
    
    # If per-mutation mode is selected, display a numeric input for each mutation to allow specific control of drug impact on each clone. 
    node_info <- get_node_labels(input$patient, rv$clones_df)
    
    tagList(
      lapply(node_ids, function(id) {
        
        current_label <- node_info$label[
          match(id, node_info$node_id)
        ]
        
        numericInput( # drug impact per mutation with values between 0 and 2 (0 = complete suppression, 1 = no effect, >1 = expansion), and a step of 0.05 for finer control
          inputId = paste0("drug_", id),
          label = current_label,
          value = effect$per_mutation[[id]],
          min = 0,
          max = 2,
          step = 0.05
        )
      })
    )
  }
})

# Observer for editing drug effect
observeEvent(input$edit_drug_effect, {
  
  req(rv$objects, rv$clones_df, input$patient, input$selected_drug)
  
  mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
  
  # If a sample is a monotimepoint -> drug event cannot be created
  if (is.null(mat) || ncol(mat) < 3) {
    
    showModal(
      modalDialog(
        title = "Drug editing not available",
        icon = icon("exclamation-triangle", class = "text-warning"),
        
        p("This patient has only one real timepoint."),
        br(),
        p("Drug impact can only be edited when at least two real timepoints exist."),
        
        tags$ul(
          tags$li("Drug events are inserted between real timepoints."),
          tags$li("With only one timepoint, no drug event can be created.")
        ),
        
        br(),
        strong("To enable drug editing:"),
        tags$ol(
          tags$li("Click 'Create new timepoint'."),
          tags$li("Then return to 'Edit drug impact'.")
        ),
        
        footer = modalButton("OK"),
        easyClose = TRUE
      )
    )
    
    return()
  }
  
  drug_cols <- seq(2, ncol(mat), by = 2)
  
  if (length(drug_cols) == 0) {
    showNotification("No drug events detected.", type = "warning")
    return()
  }
  
  drug_labels <- paste("Drug event", seq_along(drug_cols))
  
  selected_col <- drug_cols[1]
  
  node_ids <- rownames(mat)
  node_info <- get_node_labels(input$patient, rv$clones_df)
  
  key <- paste(
    input$patient,
    input$selected_drug,
    selected_col,
    sep = "_"
  )
  
  if (is.null(rv$drug_effects[[key]])) {
    
    rv$drug_effects[[key]] <- list(
      mode = "global",
      global_value = 1,
      per_mutation = setNames(rep(1, length(node_ids)), node_ids)
    )
    
  }
  
  showModal(
    modalDialog(
      title = paste("Drug impact –", input$selected_drug),
      size = "l",
      easyClose = FALSE,
      fade = TRUE,
      
      selectInput(
        "selected_drug_col",
        "Select drug event",
        choices = setNames(drug_cols, drug_labels),
        selected = selected_col
      ),
      
      helpText(
        "Drug impact rules:",
        "• 1 = no effect",
        "• <1 = clone shrinkage",
        "• >1 = clone expansion",
        "• Values are automatically rebalanced",
        "to respect clonal hierarchy constraints."
      ),
      
      h4("Mutation tree"),
      div(
        style = "
          height: 500px;
          overflow: auto;
          border: 1px solid #ddd;
          margin-bottom: 20px;
        ",
        grVizOutput("mini_tree", height = "480px")
      ),
      
      hr(),
      
      fluidRow(
        
        column(
          7,
          h4("Drug effect per clone"),
          div(
            style = "
              max-height: 400px;
              overflow-y: auto;
              padding-right: 10px;
            ",
            radioButtons(
              "drug_mode",
              "Impact mode",
              choices = c(
                "Same impact for all mutations" = "global",
                "Specific impact per mutation" = "per_mutation"
              ),
              selected = rv$drug_effects[[key]]$mode,
              inline = FALSE
            ),
            uiOutput("drug_numeric_inputs")
          )
        ),
        
        column(
          5,
          h4("Mutation summary"),
          DTOutput("mini_table")
        )
      ),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_drug_edit",
          "Apply changes",
          class = "btn-success"
        )
      )
    )
  )
  
})

# Observer for confirming drug effect edits, which updates the reactive values storing drug effects based on user input and the selected impact mode.
observeEvent(input$confirm_drug_edit, {
  
  req(input$selected_drug_col)
  
  key <- paste(
    input$patient,
    input$selected_drug,
    input$selected_drug_col,
    sep = "_"
  )
  
  mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
  node_ids <- rownames(mat)
  
  if (is.null(rv$drug_effects[[key]])) {
    rv$drug_effects[[key]] <- list(
      mode = input$drug_mode,
      global_value = 1,
      per_mutation = setNames(rep(1, length(node_ids)), node_ids)
    )
  }
  
  rv$drug_effects[[key]]$mode <- input$drug_mode
  
  if (input$drug_mode == "global") {
    
    rv$drug_effects[[key]]$global_value <- input$drug_global
    
    rv$drug_effects[[key]]$per_mutation <-
      setNames(rep(input$drug_global, length(node_ids)), node_ids)
    
  } else {
    
    for (id in node_ids) {
      rv$drug_effects[[key]]$per_mutation[[id]] <-
        input[[paste0("drug_", id)]]
    }
  }
  
  removeModal()
})

# Add drug observer to allow users to add new drugs to the dataset, which can then be selected for editing drug effects in the fishplot matrix.
observeEvent(input$add_drug, {
  
  showModal(
    modalDialog(
      title = "Add new drug",
      textInput("new_drug_name", "Drug name"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_add_drug", "Add")
      )
    )
  )
})

# Confirm add drug observer
observeEvent(input$confirm_add_drug, {
  
  req(input$new_drug_name)
  
  rv$available_drugs <- unique(
    c(rv$available_drugs, input$new_drug_name)
  )
  
  removeModal()
})

# Drug mode observer to update the drug effect reactive values when the user switches between global and per-mutation impact modes in the drug effect editing modal. 
# This ensures that the correct input values are displayed and stored based on the selected mode.
observeEvent(input$drug_mode, {
  
  req(input$selected_drug_col)
  
  key <- paste(
    input$patient,
    input$selected_drug,
    input$selected_drug_col,
    sep = "_"
  )
  
  if (!is.null(rv$drug_effects[[key]])) {
    rv$drug_effects[[key]]$mode <- input$drug_mode
  }
  
})
observeEvent(
  {
    input$selected_drug
    input$patient
  },
  {
    
    req(rv$objects, input$patient, input$selected_drug)
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    if (is.null(mat) || ncol(mat) < 3) return()
    
    drug_cols <- seq(2, ncol(mat), by = 2)
    node_ids <- rownames(mat)
    
    for (col in drug_cols) {
      
      key <- paste(
        input$patient,
        input$selected_drug,
        col,
        sep = "_"
      )
      
      if (is.null(rv$drug_effects[[key]])) {
        
        rand_val <- generate_random_global_effect()
        
        rv$drug_effects[[key]] <- list(
          mode = "global",
          global_value = rand_val,
          per_mutation = setNames(
            rep(rand_val, length(node_ids)),
            node_ids
          )
        )
      }
    }
  }
)