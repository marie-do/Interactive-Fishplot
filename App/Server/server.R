#server

server <- function(input, output, session) {
  server_env <- environment()

  source_server_file <- function(relative_path) {
    candidate_paths <- c(
      relative_path,
      file.path("app", relative_path)
    )

    existing_paths <- candidate_paths[file.exists(candidate_paths)]

    if (length(existing_paths) == 0) {
      stop(paste("Unable to find server module:", relative_path))
    }

    source(existing_paths[[1]], local = server_env)
  }
  
  rv <- reactiveValues(
    clones_df = NULL,
    metadata_df = NULL,
    objects = NULL,
    view = "fish",
    drug_effects = list(),
    zoom_level = 1,
    event_labels = list()
  )
  # Load modules
  source_server_file("server/drugs_observers.R")
  source_server_file("server/edit_sizes_observers.R")
  source_server_file("server/fishplot_matrix_reactive.R")
  source_server_file("server/get_depth_first_order.R")
  source_server_file("server/load_file.R")
  source_server_file("server/metadata_observers.R")
  source_server_file("server/mutations_observers.R")
  source_server_file("server/reindex_nodes.R")
  source_server_file("server/save_load.R")
  source_server_file("server/safe_hierarchy_rebalance.R")
  source_server_file("server/timepoints_observers.R")
  source_server_file("server/views.R")
}