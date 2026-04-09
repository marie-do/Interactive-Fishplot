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
  
  parent_eps <- 1e-6
  parent_rows_current_tp <- rv$clones_df %>%
    filter(
      sample_id == input$selected_timepoint,
      as.character(node_id) == as.character(input$new_parent)
    )
  
  parent_row_missing <- nrow(parent_rows_current_tp) == 0
  parent_was_zero <- !parent_row_missing &&
    any(is.na(parent_rows_current_tp$size_percent) |
          parent_rows_current_tp$size_percent <= 0)
  
  if (parent_row_missing) {
    parent_template <- rv$clones_df %>%
      filter(
        get_patient_id(sample_id) == input$patient,
        as.character(node_id) == as.character(input$new_parent)
      ) %>%
      slice(1)
    
    if (nrow(parent_template) == 1) {
      parent_template <- parent_template %>%
        mutate(
          sample_id = input$selected_timepoint,
          size_percent = parent_eps
        )
      
      rv$clones_df <- bind_rows(rv$clones_df, parent_template)
    }
  }
  
  if (parent_was_zero) {
    rv$clones_df <- rv$clones_df %>%
      mutate(
        size_percent = ifelse(
          sample_id == input$selected_timepoint &
            as.character(node_id) == as.character(input$new_parent) &
            (is.na(size_percent) | size_percent <= 0),
          parent_eps,
          size_percent
        )
      )
  }
  
  if (parent_row_missing || parent_was_zero) {
    showNotification(
      "Parent was missing or 0% at this timepoint and has been automatically set to a very low value.",
      type = "message",
      duration = 6
    )
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

get_mutation_first_appearance_tp_num <- function(clones_df, target_node_id, patient, patient_timepoints_df) {
  patient_ids <- sub("^(.*)-[^-]+$", "\\1", clones_df$sample_id)
  target_rows <- clones_df[
    clones_df$node_id == target_node_id &
      patient_ids == patient &
      clones_df$size_percent > 0,
  ]
  
  if (nrow(target_rows) == 0) {
    return(NA_real_)
  }
  
  joined_rows <- merge(target_rows, patient_timepoints_df, by = "sample_id")
  
  if (nrow(joined_rows) == 0) {
    return(NA_real_)
  }
  
  min(joined_rows$tp_num, na.rm = TRUE)
}

normalize_root_mutations <- function(df) {
  df$mutation[df$parent_id == "root"] <- "none"
  df
}

# Deletion helper text
output$deletion_preview <- renderText({
  
  req(input$mutation_to_delete, input$selected_timepoint)
  
  node <- input$mutation_to_delete
  patient <- input$patient
  current_tp <- input$selected_timepoint
  
  patient_timepoints_df <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == patient) %>%
    pull(sample_id) %>%
    unique() %>%
    tibble(sample_id = .) %>%
    mutate(tp_num = suppressWarnings(get_suffix_num(sample_id))) %>%
    filter(!is.na(tp_num)) %>%
    arrange(tp_num)
  
  current_tp_num <- patient_timepoints_df %>%
    filter(sample_id == current_tp) %>%
    pull(tp_num) %>%
    dplyr::first()
  
  if (is.na(current_tp_num)) {
    current_tp_num <- suppressWarnings(get_suffix_num(current_tp))
  }
  
  first_appearance_num <- get_mutation_first_appearance_tp_num(
    rv$clones_df,
    node,
    patient,
    patient_timepoints_df
  )
  
  if (is.na(first_appearance_num) || first_appearance_num == current_tp_num) {
    
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
  
  patient_timepoints_df <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == patient) %>%
    pull(sample_id) %>%
    unique() %>%
    tibble(sample_id = .) %>%
    mutate(tp_num = suppressWarnings(get_suffix_num(sample_id))) %>%
    filter(!is.na(tp_num)) %>%
    arrange(tp_num)
  
  current_tp_num <- patient_timepoints_df %>%
    filter(sample_id == current_tp) %>%
    pull(tp_num) %>%
    dplyr::first()
  
  if (is.na(current_tp_num)) {
    current_tp_num <- suppressWarnings(get_suffix_num(current_tp))
  }
  
  first_appearance_num <- get_mutation_first_appearance_tp_num(
    rv$clones_df,
    node_to_delete,
    patient,
    patient_timepoints_df
  )
  
  if (is.na(first_appearance_num) || first_appearance_num == current_tp_num) {
    
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
    
    later_timepoints <- patient_timepoints_df %>%
      filter(tp_num >= current_tp_num) %>%
      pull(sample_id)
    
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
  
  rv$clones_df <- normalize_timepoint_percentages(rv$clones_df)
  rv$clones_df <- normalize_root_mutations(rv$clones_df)
  rv$objects <- build_all_objects(rv$clones_df)
  
  removeModal()
})