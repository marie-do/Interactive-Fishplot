
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
  
  # Only mutations that are active at the current timepoint and do not have children can be deleted, to ensure biological consistency of the tree structure.
  current_active_nodes <- df_patient %>%
    filter(sample_id == input$selected_timepoint,
           size_percent > 0) %>%
    pull(node_id)
  
  deletable_clones <- df_structure %>%
    filter(
      parent_id != "root",
      !node_id %in% nodes_with_children,
      node_id %in% current_active_nodes
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
      filter(node_id != node_to_delete)
    
    # Clean the drug effect associated
    rv$drug_effects <- rv$drug_effects[
      !grepl(paste0("_", node_to_delete, "$"),
             names(rv$drug_effects))
    ]
    
    # Reindexing nodes
    rv$clones_df <- reindex_nodes(rv$clones_df)
    
    showNotification(
      paste(
        "Mutation completely deleted."
      ),
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
  rv$objects <- build_all_objects(rv$clones_df)
  
  removeModal()
})

# Observer for handling cell edits in the metadata table, which updates the reactive metadata dataframe based on user modifications.
observeEvent(input$metadata_table_cell_edit, {
  
  info <- input$metadata_table_cell_edit
  req(info)
  
  df_patient <- rv$metadata_df %>%
    filter(get_patient_id(sample_id) == input$patient)
  
  real_col <- colnames(df_patient)[info$col+1] # +1 because DT indexes start at 0
  
  real_sample_id <- df_patient$sample_id[info$row]
  
  rv$metadata_df[
    rv$metadata_df$sample_id == real_sample_id,
    real_col
  ] <- info$value
})

# Add metadata field observer
observeEvent(input$add_metadata_field, {
  
  showModal(
    modalDialog(
      title = "Add metadata field",
      textInput("new_meta_name", "Field name"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_add_meta", "Add", class = "btn-success")
      )
    )
  )
})

# Confirm add metadata field observer
observeEvent(input$confirm_add_meta, {
  
  req(input$new_meta_name)
  
  new_col <- make.names(input$new_meta_name)
  
  if (new_col %in% colnames(rv$metadata_df)) {
    showNotification("Field already exists",
                     type = "error")
    return()
  }
  
  rv$metadata_df[[new_col]] <- NA
  
  removeModal()
  
  showNotification(
    paste0(
      "Metadata field '", input$new_meta_name,
      "' added.\nClick inside the empty cells of this new column to enter values."
    ),
    type = "message",
    duration = 6
  )
})

# Delete metadata field observer
observeEvent(input$delete_metadata_field, {
  
  showModal(
    modalDialog(
      title = "Delete metadata field",
      selectInput(
        "meta_field_to_delete",
        "Select field",
        choices = setdiff(colnames(rv$metadata_df), "sample_id")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_meta",
                     "Delete",
                     class = "btn-danger")
      )
    )
  )
})

# Confirm delete metadata field observer
observeEvent(input$confirm_delete_meta, {
  
  req(input$meta_field_to_delete)
  
  rv$metadata_df[[input$meta_field_to_delete]] <- NULL
  
  removeModal()
  showNotification("Metadata field deleted",
                   type = "warning")
})