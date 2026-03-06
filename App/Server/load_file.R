# Load file server
observeEvent(input$load_file, {
  
  req(input$load_file)
  
  my_json <- jsonlite::fromJSON(
    input$load_file$datapath,
    simplifyVector = FALSE
  )
  
  # Clone extraction
  base_clones_df <- extract_all_samples(my_json) %>%
    mutate(
      node_id = as.character(node_id),
      parent_id = ifelse(is.na(parent_id), "root", as.character(parent_id)),
      size_percent = as.numeric(size_percent) * 100,
      mutation = ifelse(is.na(mutation), "none", as.character(mutation)),
      sample_id = as.character(sample_id)
    )
  
  #Metadata extraction 
  base_metadata_df <- extract_metadata(my_json)
  
  # Store in reactive values
  rv$metadata_df <- base_metadata_df
  rv$clones_df <- base_clones_df
  rv$objects <- build_all_objects(base_clones_df)
  
  rv$drug_effects <- list()   
  
  for (patient in my_json) {
    
    if (!is.null(patient$drug_effects)) {
      
      rv$drug_effects <- c(
        rv$drug_effects,
        patient$drug_effects
      )
    }
  }
  
  updateSelectInput(
    session,
    "patient",
    choices = sort(unique(get_patient_id(base_clones_df$sample_id))),
    selected = unique(get_patient_id(base_clones_df$sample_id))[1]
  )
  
  showNotification("File loaded successfully", type = "message")
})