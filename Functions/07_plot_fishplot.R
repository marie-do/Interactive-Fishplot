#1. This function plots a fishplot for a given patient, applying drug effects if specified.
#It handles both mini and full versions of the plot, adjusting margins accordingly. 
#The function checks for the presence of data, applies drug effects to the clone fractions, and ensures that the resulting matrix adheres to fishplot constraints. 
#Finally, it creates the fishplot and adds a legend if it's not a mini version.

plot_fishplot <- function(
    patient,
    concat_by_patient_hierarchy,
    all_patient_hierarchies,
    mutation_lookup,
    clone_palette,
    input,
    rv,
    mini = FALSE,
    event_labels = NULL 
) {
  
  if (mini) {
    par(
      mar = c(3, 3, 2, 3) + 0.1,
      xpd = NA
    )
  } else {
    par(
      mar = c(7, 4, 4, 4) + 0.1,
      xpd = NA
    )
  }
  
  
  mat <- concat_by_patient_hierarchy[[patient]]
  
  if (is.null(mat) || all(is.na(mat))) {
    plot.new()
    text(0.5, 0.5, paste("No fishplot for", patient))
    return()
  }
  
  parents <- all_patient_hierarchies[[patient]]$parent_index
  mat_final <- mat
  was_monotime_duplicated <- FALSE
  
  if (
    ncol(mat_final) == 2 &&
    !is.null(colnames(mat_final)) &&
    grepl("_duplicated$", colnames(mat_final)[2])
  ) {
    was_monotime_duplicated <- TRUE
  }
  
  # FIX MONOTIME BEFORE DRUG
  if (ncol(mat_final) == 1) {
    was_monotime_duplicated <- TRUE
    mat_final <- cbind(mat_final, mat_final)
  }
  
  
  # APPLY DRUG EFFECT
  if (ncol(mat_final) >= 3) {
    
    for (col in seq(2, ncol(mat_final) - 1, by = 2)) {
      
      key <- paste(patient, input$selected_drug, col, sep = "_")
      effect <- rv$drug_effects[[key]]
      
      if (!is.null(effect)) {
        
        if (effect$mode == "global") {
          
          mat_final[, col] <- mat_final[, col - 1] * effect$global_value
          
        } else {
          
          for (i in seq_len(nrow(mat_final))) {
            
            node_id <- rownames(mat_final)[i]
            
            mult <- effect$per_mutation[[node_id]]
            
            if (is.null(mult)) mult <- 1
            
            mat_final[i, col] <- mat_final[i, col - 1] * mult
          }
        }
      }
    }
  }
  
  # Apply parent-child constraint fix iteratively after drug effect until convergence
  if (ncol(mat_final) >= 2) {
    repeat {
      mat_before <- mat_final
      mat_final <- fix_parent_child_conflicts(mat_final, parents)
      mat_final <- enforce_fishplot_constraints(mat_final, parents)
      
      # Check for convergence
      if (isTRUE(all.equal(mat_final, mat_before, tolerance = 1e-10))) {
        break
      }
    }
  }
  
  mat_final <- mat_final |>
    fix_zero_reappearance()
  
  mat_final <- enforce_fishplot_constraints(mat_final, parents)
  
  mat_final[!is.finite(mat_final)] <- 0
  mat_final[mat_final < 0] <- 0
  
  
  clone_ids <- rownames(mat_final)
  fish_colors <- clone_palette[clone_ids]
  timepoints <- seq_len(ncol(mat_final))
  
  # GENERATE CLEAN TIME + DRUG LABELS
  
  n_cols <- ncol(mat_final)
  
  x_labels <- character(n_cols)
  
  if (was_monotime_duplicated && n_cols == 2) {
    x_labels <- c("t1", "duplicated t1")
  } else {
    real_tp_index <- 1
    drug_index <- 1
    
    for (i in seq_len(n_cols)) {
      
      if (i %% 2 == 1) {
        x_labels[i] <- paste0("t", real_tp_index)
        real_tp_index <- real_tp_index + 1
      } else {
        x_labels[i] <- paste0("drug_event_", drug_index)
        drug_index <- drug_index + 1
      }
    }
  }
  
  fishObj <- tryCatch({
    
    suppressWarnings(
      createFishObject(
        frac.table = mat_final,
        parents = parents,
        timepoints = timepoints,
        col = fish_colors
      )
    )
    
  }, error = function(e) {
    
    plot.new()
    text(
      0.5, 0.5,
      paste("Invalid configuration after drug effect.\n", e$message),
      cex = 1
    )
    
    return(NULL)
  })
  
  if (is.null(fishObj)) return()
  
  
  fishObj <- layoutClones(
    fishObj,
    separate.independent.clones = FALSE
  )
  
  if (!mini) {
    
    if (!is.null(event_labels) && length(event_labels) == length(timepoints)) {
      vlab_final <- event_labels
    } else {
      vlab_final <- x_labels
    }
    
  } else {
    vlab_final <- NULL
  }
  
  fishPlot(
    fishObj,
    shape = "spline",
    vlines = if (mini) NULL else timepoints,
    vlab = vlab_final,
    title.btm = if (mini) "" else patient
  )
  
  if (!mini) {
    clone_labels <- sapply(
      clone_ids,
      get_clone_label,
      patient_id = patient,
      mutation_lookup = mutation_lookup
    )
    
    clone_labels <- as.character(clone_labels)
    add_fishplot_legend(fishObj, clone_labels)
  }
}