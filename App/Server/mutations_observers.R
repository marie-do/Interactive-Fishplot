# New mutation observer
observeEvent(input$add_mutation, {
  
  req(rv$clones_df)
  
  df_patient <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == input$patient)
  
  df_tree <- df_patient %>%
    distinct(node_id, parent_id, mutation) %>%
    arrange(as.numeric(node_id))
  
  showModal(
    modalDialog(
      title = paste("Add new mutation – Patient", input$patient),
      size = "l",
      
      div(
        style = "background-color:#f8f9fa; padding:15px; border-radius:8px; margin-bottom:20px;",
        
        h4("Current editing context"),
        
        p(
          strong("IMPORTANT : You are currently editing timepoint: "),
          strong(input$selected_timepoint)
        ),
        
        br(),
        
        p(
          strong("To change timepoint:"),
          br(),
          "Use the timepoint selector in the left sidebar ",
          "and choose another timepoint before clicking 'Create mutation'."
        ),
        
        hr(),
        
        h4("How to create a new mutation"),
        
        tags$ol(
          tags$li(
            strong("Select a parent node: "),
            "Choose the clone from which the new mutation evolves."
          ),
          tags$li(
            strong("Enter the mutation name: "),
            "Example: TP53 R175H."
          ),
          tags$li(
            strong("Set the percentage for the current timepoint: "),
            "This represents the clone size at ",
            strong(input$selected_timepoint), "."
          )
        ),
        
        hr(),
        
        p(
          strong("Important biological rule:"),
          br(),
          "• The mutation will appear at the current timepoint.",
          br(),
          "• It will be automatically set to 0% at all other timepoints.",
          br(),
          "• Parent-child hierarchy constraints will be enforced automatically."
        )
      ),
      
      
      fluidRow(
        column(
          6,
          h4("Current mutation tree"),
          grVizOutput("mini_tree", height = "400px")
        ),
        column(
          6,
          h4("Mutation summary"),
          DTOutput("mini_table")
        )
      ),
      
      hr(),
      h4("New mutation parameters"),
      
      selectInput(
        "new_parent",
        "Select parent node",
        choices = df_tree$node_id
      ),
      
      textInput(
        "new_mutation_name",
        "Mutation name",
        placeholder = "e.g. TP53 R175H"
      ),
      
      numericInput(
        "new_mutation_percent",
        "Percentage (current timepoint)",
        value = 0,
        min = 0,
        max = 100,
        step = 0.1
      ),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_add_mutation",
          "Add mutation",
          class = "btn-success"
        )
      )
    )
  )
})

# Confirm add mutation observer
observeEvent(input$confirm_add_mutation, {
  
  req(rv$clones_df)
  removeModal()
  
  # new node ID generation (incremental based on the maximum existing node ID to ensure uniqueness)
  existing_nodes <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == input$patient) %>%
    distinct(node_id) %>%
    pull(node_id) %>%
    as.numeric()
  
  new_id <- max(existing_nodes) + 1
  new_id <- as.character(new_id)
  
  new_id <- as.character(new_id)
  
  # New row creation
  if (input$new_mutation_name == "") {
    showNotification("Mutation name cannot be empty", type = "error")
    return()
  }
  
  new_row <- rv$clones_df %>%
    filter(sample_id == input$selected_timepoint) %>%
    slice(1) %>%
    mutate(
      node_id = new_id,
      parent_id = input$new_parent,
      mutation = input$new_mutation_name,
      size_percent = input$new_mutation_percent
    )
  
  # Error handling
  if (input$new_parent == new_id) {
    showNotification("Parent cannot be itself", type = "error")
    return()
  }
  
  patient_samples <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == input$patient) %>%
    pull(sample_id) %>%
    unique()
  
  # When a new mutation is added, at the selected timepoint it will have the user-defined percentage, but at all other timepoints it will start with 0%. 
  all_new_rows <- lapply(patient_samples, function(s) {
    new_row %>%
      mutate(
        sample_id = s,
        size_percent = ifelse(
          s == input$selected_timepoint,
          input$new_mutation_percent,
          0
        )
      )
  }) %>%
    bind_rows()
  
  rv$clones_df <- bind_rows(rv$clones_df, all_new_rows)
  # Normalization + objects rebuilding to ensure the new mutation is properly integrated. 
  rv$clones_df <- normalize_timepoint_percentages(rv$clones_df)
  rv$objects <- build_all_objects(rv$clones_df)
})

# Deletion helper text
output$deletion_preview <- renderText({
  
  req(input$mutation_to_delete, input$selected_timepoint)
  
  node <- input$mutation_to_delete
  patient <- input$patient
  current_tp <- input$selected_timepoint
  
  first_appearance <- rv$clones_df %>%
    filter(
      node_id == node,
      get_patient_id(sample_id) == patient,
      size_percent > 0
    ) %>%
    pull(sample_id) %>%
    sort() %>%
    .[1]
  
  if (!is.na(first_appearance) && first_appearance == current_tp) {
    
    paste(
      "This mutation first appears at the current timepoint.",
      "\n→ It will be completely removed from the dataset.",
      "\n→ It will disappear from all timepoints.",
      "\n→ The evolutionary structure will be updated."
    )
    
  } else {
    
    paste(
      "This mutation existed before the current timepoint.",
      "\n→ It cannot biologically disappear.",
      "\n→ It will be reduced to near-zero from this timepoint onward.",
      "\n→ Previous timepoints will remain unchanged."
    )
  }
})

#Observer for mutation deletion
observeEvent(input$delete_mutation, {
  
  req(rv$clones_df, input$patient)
  
  df_patient <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == input$patient)
  
  df_structure <- df_patient %>%
    distinct(node_id, mutation, parent_id)
  
  nodes_with_children <- df_structure$parent_id
  
  active_nodes <- df_patient %>%
    group_by(node_id) %>%
    summarise(total = sum(size_percent, na.rm = TRUE)) %>%
    filter(total > 0) %>%
    pull(node_id)
  
  deletable_clones <- df_structure %>%
    filter(
      parent_id != "root",
      !node_id %in% nodes_with_children,
      node_id %in% active_nodes
    ) %>%
    arrange(as.numeric(node_id))
  
  if (nrow(deletable_clones) == 0) {
    showNotification(
      "No deletable mutations available.\nOnly active leaf mutations (without children) can be removed.",
      type = "warning"
    )
    return()
  }
  
  showModal(
    modalDialog(
      title = "Delete mutation",
      size = "l",
      
      p(
        strong("IMPORTANT : You are currently editing timepoint: "),
        strong(input$selected_timepoint)
      ),
      
      br(),
      
      p(
        strong("To change timepoint:"),
        br(),
        "Use the timepoint selector in the left sidebar ",
        "and choose another sample before clicking 'Create mutation'."
      ),
      
      h4("Deletion rules"),
      tags$ul(
        tags$li("Only leaf mutations can be deleted."),
        tags$li("Root and internal nodes are protected."),
        tags$li("Historical mutations cannot fully disappear."),
        tags$li("New mutations (first appearing now) can be completely removed.")
      ),
      
      hr(),
      
      selectInput(
        "mutation_to_delete",
        "Select mutation",
        choices = setNames(
          deletable_clones$node_id,
          paste0("Node ", deletable_clones$node_id,
                 " | ", deletable_clones$mutation)
        )
      ),
      
      hr(),
      
      h4("What will happen?"),
      verbatimTextOutput("deletion_preview"),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_delete_mutation",
          "Confirm deletion",
          class = "btn-danger"
        )
      )
    )
  )
})

# Observer for confirming mutation deletion with different scenarios based on whether the mutation first appears at the current timepoint or existed before.
observeEvent(input$confirm_delete_mutation, {
  
  req(rv$clones_df, input$mutation_to_delete, input$selected_timepoint)
  
  node_to_delete <- input$mutation_to_delete
  patient <- input$patient
  current_tp <- input$selected_timepoint
  
  patient_timepoints <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == patient) %>%
    pull(sample_id) %>%
    unique() %>%
    sort()
  
  first_appearance <- rv$clones_df %>%
    filter(
      node_id == node_to_delete,
      get_patient_id(sample_id) == patient,
      size_percent > 0
    ) %>%
    pull(sample_id) %>%
    sort() %>%
    .[1]
  
  if (!is.na(first_appearance) && first_appearance == current_tp) {
    
    rv$clones_df <- rv$clones_df %>%
      filter(
        !(node_id == node_to_delete &
            get_patient_id(sample_id) == patient)
      )
    
    showNotification(
      "Mutation completely deleted.",
      type = "message",
      duration = 6
    )
    
  } else {
    
    later_timepoints <- patient_timepoints[
      patient_timepoints >= current_tp
    ]
    
    rv$clones_df <- rv$clones_df %>%
      mutate(
        size_percent = ifelse(
          node_id == node_to_delete &
            get_patient_id(sample_id) == patient &
            sample_id %in% later_timepoints,
          1e-6,
          size_percent
        )
      )
    
    showNotification(
      paste(
        "Mutation biologically suppressed.",
        "It has been reduced to near-zero from this timepoint onward."
      ),
      type = "warning",
      duration = 8
    )
  }
  
  # clean empty clones
  rv$clones_df <- rv$clones_df %>%
    group_by(patient = get_patient_id(sample_id)) %>%
    group_modify(~{
      
      df <- .x
      
      parent_nodes <- unique(df$parent_id)
      
      df %>%
        group_by(node_id) %>%
        filter(
          sum(size_percent, na.rm = TRUE) > 0 |
            node_id %in% parent_nodes |
            parent_id == "root"   
        ) %>%
        ungroup()
      
    }) %>%
    ungroup()

  
  rv$objects <- build_all_objects(rv$clones_df)
  
  removeModal()
})