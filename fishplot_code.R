set.seed(123)

source("dependencies.R")
source("fishplot_set_up.R")

# Load custom functions
invisible(
  lapply(
    list.files(pattern = "^[0-9].*\\.R$"),
    source
  )
)

#1. Load Data
json_path <- "C:/Users/DOGO Marie/Desktop/Interactive-Fishplot/data/trees_aml_morita.json"


my_json <- fromJSON(
  json_path,
  simplifyVector = FALSE
)

# 2. Build & clean clones dataframe
clones_df <- extract_all_samples(my_json) |>
  mutate(
    parent_id   = ifelse(is.na(parent_id), "root", parent_id),
    size_percent = suppressWarnings(as.numeric(size_percent)),
    mutation    = ifelse(is.na(mutation), "none", mutation)
  )

#fishplot colors 
global_clone_palette <- build_global_clone_palette(clones_df)

save(clones_df, file = "clones_df.RData")
write.csv(clones_df, "clones_df.csv", row.names = FALSE)

# 3. Build Cumulative Matrices
all_matrices <- compute_all_cumulative_matrices(clones_df)

# Replace NA first rows with 0 (fishplot safety)
all_matrices <- lapply(all_matrices, function(mat) {
  if (is.matrix(mat) && nrow(mat) >= 1 && all(is.na(mat[1, ]))) {
    mat[1, ] <- 0
  }
  mat
})

# 4. Classify samples
classified <- classify_samples(clones_df)

cat("Monotime samples :", length(classified$monotime), "\n")
cat("Multitime samples:", length(classified$multitime), "\n")

# 5. Build hierarchies per patient
patient_samples_map <- split(
  unique(clones_df$sample_id),
  get_patient_id(unique(clones_df$sample_id))
)

all_patient_hierarchies <- lapply(
  patient_samples_map,
  build_patient_hierarchy,
  clones_df = clones_df
)


# 6. Concatenate timepoints per patient

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


# 7. Generate fishplots

mutation_lookup <- build_mutation_lookup(clones_df)

outdir <- "fishplots_png"
dir.create(outdir, showWarnings = FALSE)

final_fishplot_matrices <- list()

for (patient_id in names(concat_by_patient)) {
  
  mat <- concat_by_patient[[patient_id]]
  
  if (is.null(mat) || all(is.na(mat))) {
    final_fishplot_matrices[[patient_id]] <- NA
    next
  }
  
  parents <- all_patient_hierarchies[[patient_id]]$parent_index
  
  # Apply corrections
  mat_final <- mat |>
    fix_zero_reappearance() |>
    enforce_fishplot_constraints(parents)
  
  if (ncol(mat_final) == 1) {
    mat_final <- cbind(mat_final, mat_final)
  }
  
  mat_final[!is.finite(mat_final)] <- 0
  mat_final[mat_final < 0] <- 0
  
  final_fishplot_matrices[[patient_id]] <- mat_final
  
  # Fishplot inputs
  clone_ids  <- rownames(mat_final)
  fish_colors <- global_clone_palette[clone_ids]
  timepoints <- seq_len(ncol(mat_final))
  
  # Safe fish object creation
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

# 8. Export final matrices and metadata

metadata_df <- extract_metadata(my_json)

save(metadata_df, file = "metadata.RData")

final_fishplot_matrices