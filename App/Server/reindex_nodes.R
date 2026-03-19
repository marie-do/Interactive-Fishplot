# Helper function to reindex node IDs after full deletion of a mutation, to ensure that node IDs remain consecutive and consistent.
# This is important for the fishplot construction and to avoid issues with missing node IDs.
reindex_nodes_per_patient <- function(df) {
  
  df %>%
    group_by(patient = get_patient_id(sample_id)) %>%
    group_modify(~{
      
      unique_nodes <- .x %>%
        distinct(node_id) %>%
        arrange(as.numeric(node_id)) %>%
        pull(node_id)
      
      new_ids <- as.character(seq_along(unique_nodes))
      id_map <- setNames(new_ids, unique_nodes)
      
      .x %>%
        mutate(
          node_id = id_map[node_id],
          parent_id = case_when(
            parent_id == "root" ~ "root",
            parent_id %in% names(id_map) ~ id_map[parent_id],
            TRUE ~ parent_id   
          )
        )
    }) %>%
    ungroup()
}