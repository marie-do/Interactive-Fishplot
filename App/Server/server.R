#server

server <- function(input, output, session) {
  
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
  source("Server/drugs_observers.R", local = TRUE)
  source("server/edit_sizes_observers.R", local = TRUE)
  source("server/fishplot_matrix_reactive.R", local = TRUE)
  source("server/get_depth_first_order.R", local = TRUE)
  source("server/load_file.R", local = TRUE)
  source("server/metadata_observers.R", local = TRUE)
  source("server/reindex_nodes.R", local = TRUE)
  source("server/save_load.R", local = TRUE)
  source("server/safe_hierarchy_rebalance.R", local = TRUE)
  source("server/timepoints_observers.R", local = TRUE)
  source("server/views.R", local = TRUE)
}