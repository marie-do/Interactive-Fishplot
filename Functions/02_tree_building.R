#1. Tree building
# Reconstructs a hierarchical tree structure from the flat dataframe
build_tree <- function(df_sample, node_id) {
  node_row <- df_sample[df_sample$node_id == node_id, , drop = FALSE]
  if (nrow(node_row) == 0) return(NULL)
  
  children_ids <- df_sample$node_id[df_sample$parent_id == node_id]
  children_list <- lapply(children_ids, function(ch) {
    build_tree(df_sample, ch)
  })
  children_list <- compact(children_list)
  
  list(
    node_id = as.character(node_row$node_id[1]),
    mutation = node_row$mutation[1],
    size_percent = suppressWarnings(as.numeric(node_row$size_percent[1])),
    freq = suppressWarnings(as.numeric(node_row$size_percent[1]) * 100),
    children = children_list
  )
}

#2. Mapping function
# For a given sample, identifies root nodes (those with parent_id "root")
# and builds a tree for each root using build_tree. Returns a list of trees for the sample.
mapping <- function(df_sample) {
  roots <- df_sample$node_id[df_sample$parent_id == "root"]
  roots <- as.character(roots)
  
  map(roots, ~ build_tree(df_sample, .x)) %>%
    compact()
} # Result -> a list of tree structures (nested lists) for each root node in the sample

#3. Nodes extraction from tree
# Extracts all nodes from a built tree structure, adding depth
extract_nodes_tree <- function(tree, depth = 1) {
  tibble(
    node_id = tree$node_id,
    freq = tree$freq,
    depth = depth
  ) %>%
    bind_rows(
      map_dfr(tree$children, extract_nodes_tree, depth = depth + 1)
    )
}

#4. Descendants collection
#  Returns all descendant node IDs of a node (including itself)
collect_descendants <- function(tree) {
  c(
    tree$node_id,
    unlist(map(tree$children, collect_descendants))
  )
}

#5. Descendants mapping
# Creates a table mapping each node to all its descendants
get_desc_map <- function(tree) {
  tibble(
    node_id = tree$node_id,
    descendants = list(collect_descendants(tree))
  ) %>%
    bind_rows(
      map_dfr(tree$children, get_desc_map)
    )
}