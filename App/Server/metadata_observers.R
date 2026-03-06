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

observeEvent(input$Metadata, {
  rv$view <- "metadata"
})