#Reactive Matrix : Construction of the cummulative matrix with all the rules applied (drug effects, zero reappearance, parent-child constraints, root constraint)

fishplot_matrix_reactive <- reactive({
  
  req(input$patient, rv$objects)
  
  mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
  if (is.null(mat)) return(NULL)
  
  parents <- rv$objects$all_patient_hierarchies[[input$patient]]$parent_index
  mat_final <- mat
  
  # Fix monotimepoints with single column by duplicating it (fishplot requires at least 2 timepoints)
  if (ncol(mat_final) == 1) {
    
    original_name <- colnames(mat_final)[1]
    
    mat_final <- cbind(mat_final, mat_final)
    
    colnames(mat_final) <- c(
      original_name,
      paste0(original_name, "_duplicated")
    )
  }
  
  # Apply drug effects
  # Only apply drug effects if REAL multiple timepoints exist
  real_timepoints <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == input$patient) %>%
    pull(sample_id) %>%
    unique()
  
  if (length(real_timepoints) >= 2) {
    
    for (col in seq(2, ncol(mat_final) - 1, by = 2)) { #drug event columns every 2 columns, so that each timepoint are separated by a drug event.
      key <- paste(input$patient,
                   input$selected_drug,
                   col,
                   sep = "_")
      
      effect <- rv$drug_effects[[key]] # 2 modes : global or per-mutation. 
      #The observeEvent for drug effect editing will update this reactive value with the new drug effect parameters, which will then be applied here to the matrix.
      
      if (!is.null(effect)) { 
        
        if (effect$mode == "global") {
          
          mat_final[, col] <- mat_final[, col - 1] * effect$global_value # -> global mode : the same multiplier is applied to all mutations
          
        } else {
          
          for (i in seq_len(nrow(mat_final))) {
            
            node_id <- rownames(mat_final)[i]
            if (!node_id %in% names(effect$per_mutation)) {
              mult <- 1
            } else {
              mult <- effect$per_mutation[[node_id]]
            }
            mat_final[i, col] <- mat_final[i, col - 1] * mult # -> per-mutation mode : each mutation has its own multiplier
          }
        }
      }
    }
  }
  
  mat_final <- fix_zero_reappearance(mat_final)
  mat_final <- safe_hierarchy_rebalance(mat_final, parents)
  
  mat_final[mat_final < 0] <- 0
  mat_final[!is.finite(mat_final)] <- 0
  
  return(mat_final)
})