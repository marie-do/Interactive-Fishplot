# Save file server
output$save_file <- downloadHandler(
  filename = function() {
    paste0("modified_fishplot_data_", Sys.Date(), ".json")
  },
  
  content = function(file) {
    
    req(rv$clones_df)
    
    df <- rv$clones_df %>%
      mutate(
        patient_id = get_patient_id(sample_id),
        size_percent = size_percent / 100
      )
    
    meta_df <- rv$metadata_df
    
    if (is.null(meta_df)) {
      meta_df <- data.frame(sample_id = unique(df$sample_id))
    }
    
    json_structure <- df %>%
      group_split(patient_id) %>%
      lapply(function(patient_df) {
        
        patient_id_val <- unique(patient_df$patient_id)
        
        patient_drug_effects <- rv$drug_effects[
          startsWith(
            names(rv$drug_effects),
            paste0(patient_id_val, "_")
          )
        ]
        
        if (is.null(patient_drug_effects)) {
          patient_drug_effects <- list()
        }
        
        samples_list <- patient_df %>%
          group_split(sample_id) %>%
          lapply(function(sample_df) {
            
            sample_id_val <- unique(sample_df$sample_id)
            
            meta_row <- meta_df %>%
              filter(sample_id == sample_id_val) %>%
              select(-sample_id)
            
            list(
              sample_id = sample_id_val,
              
              metadata = if (nrow(meta_row) > 0) {
                as.list(meta_row[1, ])
              } else {
                list()
              },
              
              clones = sample_df %>%
                arrange(as.numeric(node_id)) %>%
                select(node_id, parent_id, mutation, size_percent) %>%
                mutate(node_id = as.numeric(node_id)) %>%
                split(seq_len(nrow(.))) %>%
                lapply(as.list)
            )
          })
        
        list(
          patient_id = patient_id_val,
          samples = samples_list,
          drug_effects = patient_drug_effects  
        )
      })
    
    jsonlite::write_json(
      json_structure,
      file,
      pretty = TRUE,
      auto_unbox = TRUE,
      null = "null"
    )
  }
)