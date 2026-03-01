set.seed(123)

source("dependencies.R")

invisible(
  lapply(
    list.files(pattern = "^[0-9].*\\.R$"),
    source
  )
)

json_path <- file.path("data", "trees_aml_morita.json")
data_filename <- tools::file_path_sans_ext(basename(json_path))

outdir <- paste0("fishplots_png_", data_filename)
cache_file <- paste0("fishplot_cache_", data_filename, ".RData")

if (file.exists(cache_file)) {
  
  message("Cache found. Loading precomputed objects...")
  load(cache_file)
  
} else {
  
  message("No cache found. Computing everything...")
  
  dir.create(outdir, showWarnings = FALSE)
  
  my_json <- fromJSON(
    json_path,
    simplifyVector = FALSE
  )
  
  clones_df <- extract_all_samples(my_json) |>
    mutate(
      parent_id   = ifelse(is.na(parent_id), "root", parent_id),
      size_percent = suppressWarnings(as.numeric(size_percent)),
      mutation    = ifelse(is.na(mutation), "none", mutation)
    )
  
  global_clone_palette <- build_global_clone_palette(clones_df)
  
  all_matrices <- compute_all_cumulative_matrices(clones_df)
  
  all_matrices <- lapply(all_matrices, function(mat) {
    if (is.matrix(mat) && nrow(mat) >= 1 && all(is.na(mat[1, ]))) {
      mat[1, ] <- 0
    }
    mat
  })
  
  classified <- classify_samples(clones_df)
  
  patient_samples_map <- split(
    unique(clones_df$sample_id),
    get_patient_id(unique(clones_df$sample_id))
  )
  
  all_patient_hierarchies <- lapply(
    patient_samples_map,
    build_patient_hierarchy,
    clones_df = clones_df
  )
  
  concat_by_patient <- lapply(
    names(patient_samples_map),
    function(p) {
      concat_patient_timepoints_hierarchy(
        patient = p,
        multitime_by_patient = patient_samples_map,
        all_matrices = all_matrices,
        all_patient_hierarchies = all_patient_hierarchies
      )
    }
  )
  
  names(concat_by_patient) <- names(patient_samples_map)
  
  mutation_lookup <- build_mutation_lookup(clones_df)
  
  final_fishplot_matrices <- list()
  
  for (patient_id in names(concat_by_patient)) {
    
    mat <- concat_by_patient[[patient_id]]
    
    if (is.null(mat) || all(is.na(mat))) {
      final_fishplot_matrices[[patient_id]] <- NA
      next
    }
    
    parents <- all_patient_hierarchies[[patient_id]]$parent_index
    
    mat_final <- mat |>
      fix_zero_reappearance() |>
      enforce_fishplot_constraints(parents)
    
    if (ncol(mat_final) == 1) {
      mat_final <- cbind(mat_final, mat_final)
    }
    
    mat_final[!is.finite(mat_final)] <- 0
    mat_final[mat_final < 0] <- 0
    
    final_fishplot_matrices[[patient_id]] <- mat_final
    
    clone_ids  <- rownames(mat_final)
    fish_colors <- global_clone_palette[clone_ids]
    timepoints <- seq_len(ncol(mat_final))
    
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
      message("Fishplot failed for patient ", patient_id)
      return(NULL)
    })
    
    if (is.null(fishObj)) next
    
    fishObj <- layoutClones(
      fishObj,
      separate.independent.clones = FALSE
    )
    
    clone_labels <- vapply(
      clone_ids,
      get_clone_label,
      character(1),
      patient_id = patient_id,
      mutation_lookup = mutation_lookup
    )
    
    png(
      file.path(outdir, paste0(patient_id, ".png")),
      width = 2400,
      height = 1600,
      res = 200
    )
    
    par(mar = c(7, 4, 4, 4) + 0.1)
    
    fishPlot(
      fishObj,
      shape = "spline",
      vlines = timepoints,
      vlab = paste0("t", timepoints),
      title.btm = patient_id
    )
    
    add_fishplot_legend(fishObj, clone_labels)
    
    dev.off()
    
    message("Fishplot saved for patient: ", patient_id)
  }
  
  metadata_df <- extract_metadata(my_json)
  
  save(
    final_fishplot_matrices,
    clones_df,
    metadata_df,
    global_clone_palette,
    file = cache_file
  )
  
  message("All objects cached successfully.")
}

final_fishplot_matrices