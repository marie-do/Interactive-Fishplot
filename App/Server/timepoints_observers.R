# Update timepoint list
observeEvent(input$patient, {
  
  req(rv$clones_df)
  
  tp <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == input$patient) %>%
    pull(sample_id) %>%
    unique() %>%
    sort()
  
  updateSelectInput(
    session,
    "selected_timepoint",
    choices = tp,
    selected = tp[1]
  )
})

# Adding new timepoint
observeEvent(input$new_timepoint, {
  
  req(rv$clones_df) # Security : only allow adding timepoints if data is loaded
  
  existing <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == input$patient) %>%
    pull(sample_id) %>%
    unique()
  
  nums <- as.numeric(str_extract(existing, "\\d+$"))
  nums <- nums[!is.na(nums)]
  new_idx <- ifelse(length(nums) == 0, 1, max(nums) + 1) # Find the highest existing index and add 1, or start at 1 if none exist
  new_sample <- paste0(input$patient, "-", sprintf("%03d", new_idx))
  
  base_rows <- rv$clones_df %>%
    filter(sample_id == existing[1]) %>%
    mutate(
      sample_id = new_sample,
      size_percent = 0
    ) # Copied from the first existing timepoint of the patient, but with 0 percentages
  
  rv$clones_df <- bind_rows(rv$clones_df, base_rows)
  rv$objects <- build_all_objects(rv$clones_df)
  rand_val <- generate_random_global_effect()
  
  #Initialize drug effects for new drug events, so that 
  mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
  
  if (!is.null(mat) && ncol(mat) >= 3) {
    
    drug_cols <- seq(2, ncol(mat), by = 2)
    node_ids <- rownames(mat)
    
    for (drug in rv$available_drugs) {
      
      for (col in drug_cols) {
        
        key <- paste(input$patient, drug, col, sep = "_")
        
        if (is.null(rv$drug_effects[[key]])) {
          
          rand_val <- generate_random_global_effect()
          
          rv$drug_effects[[key]] <- list( #global mode + per mutation mode are ready to use. Initialized with random values, but can be updated if the user edits drug effects.
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
  }
  
  updateSelectInput(
    session,
    "selected_timepoint",
    choices = c(existing, new_sample),
    selected = new_sample
  )
})

#Observer for deleting timepoint
observeEvent(input$delete_timepoint, {
  
  req(rv$clones_df, input$selected_timepoint)
  
  showModal(
    modalDialog(
      title = "Confirm deletion",
      paste(
        "Are you sure you want to delete timepoint",
        input$selected_timepoint,
        "?"
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_delete_timepoint",
          "Delete",
          class = "btn-danger"
        )
      )
    )
  )
})

# observer for confirming timepoint deletion
observeEvent(input$confirm_delete_timepoint, {
  
  req(rv$clones_df, input$selected_timepoint, input$patient)
  
  tp_to_delete <- input$selected_timepoint
  patient <- input$patient
  
  rv$clones_df <- rv$clones_df %>%
    filter(sample_id != tp_to_delete)
  
  rv$clones_df <- normalize_timepoint_percentages(rv$clones_df)
  rv$objects <- build_all_objects(rv$clones_df)
  
  remaining_timepoints <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == patient) %>%
    pull(sample_id) %>%
    unique() %>%
    sort()
  
  updateSelectInput(
    session,
    "selected_timepoint",
    choices = remaining_timepoints,
    selected = if (length(remaining_timepoints) > 0)
      remaining_timepoints[1]
    else
      NULL
  )
  
  removeModal()
  
  showNotification("Timepoint deleted", type = "warning")
})

observeEvent(input$rename_events, {
  
  req(input$patient, rv$objects)
  
  mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
  req(mat)
  
  current_labels <- colnames(mat)
  
  showModal(
    modalDialog(
      title = paste("Rename events – Patient", input$patient),
      size = "l",
      easyClose = FALSE,
      
      tagList(
        lapply(seq_along(current_labels), function(i) {
          
          textInput(
            inputId = paste0("event_label_", i),
            label = paste("Event", current_labels[i]),
            value = ifelse(
              !is.null(rv$event_labels[[input$patient]]) &&
                length(rv$event_labels[[input$patient]]) >= i,
              rv$event_labels[[input$patient]][i],
              current_labels[i]
            )
          )
        })
      ),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "confirm_event_rename",
          "Apply",
          class = "btn-success"
        )
      )
    )
  )
})

observeEvent(input$confirm_event_rename, {
  
  req(input$patient, rv$objects)
  
  mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
  req(mat)
  
  new_labels <- sapply(seq_len(ncol(mat)), function(i) {
    input[[paste0("event_label_", i)]]
  })
  
  if (any(new_labels == "") || any(duplicated(new_labels))) {
    showNotification(
      "Event names must be unique and non-empty.",
      type = "error"
    )
    return()
  }
  
  rv$event_labels[[input$patient]] <- new_labels
  
  removeModal()
  
  showNotification(
    "Event labels updated.",
    type = "message"
  )
})