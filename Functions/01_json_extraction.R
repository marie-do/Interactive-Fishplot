# 1. Nodes extraction
# Recursively traverses a tree node from the JSON and converts it into a flat table
extract_nodes <- function(node, parent = NA) {
  mutation_label <- NA_character_
  
  if (is.null(node$gene_events)) { # Mutation label extraction
    NA
  } else {
    mutation_label <- node$gene_events %>%
      imap(~ paste(.y, .x$SNV)) %>%   # gene + SNV
      unlist() %>%
      paste(collapse = " | ")
  }
  
  size_p <- if (is.null(node$size_percent)) NA else node$size_percent # Size_percent extraction
  
  # Creation of one tibble row for the current node
  df <- tibble(
    node_id = as.character(node$node_id),
    parent_id = as.character(parent),
    size_percent = size_p,
    mutation = mutation_label
  )
  
  # If the node has children -> recursively extract them and bind to the current tibble
  # steps : calls extract_nodes on each child, passes current node_id as parent, binds all child tibbles together.
  if (!is.null(node$children)) {
    child_list <- node$children
    
    if (is.data.frame(child_list)) {
      child_list <- split(child_list, seq(nrow(child_list)))
    }
    
    df_children <- bind_rows(
      map(child_list, extract_nodes, parent = node$node_id)
    )
    
    df <- bind_rows(df, df_children)
  }
  
  df
} # Result -> a tibble with columns node_id, parent_id, size_percent, mutation for all nodes in the tree
# = A flat dataframe representing the entire subtree

#2.All samples extraction
#For each samples, checks if it has a tree, if yes -> extract nodes from its tree + add sample_id column.
extract_all_samples <- function(json_data) {
  bind_rows(
    lapply(names(json_data), function(sample) {
      if (is.null(json_data[[sample]]$tree)) return(NULL)
      
      extract_nodes(json_data[[sample]]$tree) %>%
        mutate(sample_id = sample)
    })
  )
} # Result -> a tibble with columns node_id, parent_id, size_percent, mutation, sample_id for all nodes across all samples

#3. Metadata extraction
extract_metadata <- function(json_data) {
  bind_rows(
    lapply(names(json_data), function(sample) {
      meta <- json_data[[sample]]$metadata
      if (is.null(meta)) return(NULL)
      
      tibble(
        sample_id = sample,
        !!!meta
      )
    })
  )
}

#4. Creates mapping : grouped and deduplicated
build_mutation_lookup <- function(clones_df) {
  clones_df %>%
    filter(mutation != "none") %>%
    group_by(sample_id, node_id) %>%
    summarise(
      mutation = first(mutation),
      .groups = "drop"
    )
}