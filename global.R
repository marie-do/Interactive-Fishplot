
# Global variables and functions for the fishplot app
set.seed(123)

# Core dependencies
source("dependencies.R")

# If fishplot is installed, silence fishplot_set_up.R. If it's not, install fishplot manually.


if (requireNamespace("fishplot", quietly = TRUE)) {
  library(fishplot)
} else {
  source("fishplot_set_up.R")
}

source("fishplot_code.R")
source("diagram_shiny.R")


#Main function to build all objects needed for the app
build_all_objects <- function(clones_df) {
  
  clone_ids <- clones_df %>%
    distinct(node_id) %>%
    arrange(as.numeric(node_id)) %>%
    pull(node_id)
  
  clone_palette_named <- setNames(
    rep(clone_palette, length.out = length(clone_ids)),
    clone_ids
  )
  
  mutation_lookup <- build_mutation_lookup(clones_df)
  all_matrices <- compute_all_cumulative_matrices(clones_df)
  
  patient_samples_map <- split(
    unique(clones_df$sample_id),
    get_patient_id(unique(clones_df$sample_id))
  )
  
  all_patient_hierarchies <- lapply(
    names(patient_samples_map),
    function(p)
      build_patient_hierarchy(patient_samples_map[[p]], clones_df)
  )
  names(all_patient_hierarchies) <- names(patient_samples_map)
  
  concat_by_patient_hierarchy <- lapply(
    names(patient_samples_map),
    function(p)
      concat_patient_timepoints_hierarchy(
        p,
        patient_samples_map,
        all_matrices,
        all_patient_hierarchies
      )
  )
  names(concat_by_patient_hierarchy) <- names(patient_samples_map)
  
  list(
    mutation_lookup = mutation_lookup,
    all_matrices = all_matrices,
    patient_samples_map = patient_samples_map,
    all_patient_hierarchies = all_patient_hierarchies,
    concat_by_patient_hierarchy = concat_by_patient_hierarchy,
    clone_palette = clone_palette_named
  )
}

generate_random_global_effect <- function() {
  runif(1, 0.2, 1.4)
}

get_node_labels <- function(patient, clones_df) {
  
  clones_df %>%
    filter(get_patient_id(sample_id) == patient) %>%
    distinct(node_id, mutation) %>%
    arrange(as.numeric(node_id)) %>%
    mutate(
      label = paste0(
        "Node ", node_id,
        "  |  ",
        mutation
      )
    ) %>%
    select(node_id, label)
}

# Intro steps from de csv file
intro_steps <- read_delim(
  "fishplot_tour.csv",
  delim = ";",
  col_types = cols()
)

get_step <- function(n) {
  step_text <- intro_steps$text[intro_steps$step == n]
  if (length(step_text) == 0) return("")
  step_text[[1]]
}

normalize_timepoint_percentages <- function(df) {
  
  df %>%
    group_by(sample_id) %>%
    mutate(
      total = sum(size_percent, na.rm = TRUE),
      size_percent = ifelse(
        total > 100,
        size_percent / total * 100,
        size_percent
      )
    ) %>%
    ungroup() %>%
    select(-total)
}