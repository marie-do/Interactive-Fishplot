#1. Cumulative matrix construction
# Builds a cumulative frequentcy matrix for a sample,
# where each entry (i, j) represents the cumulative frequency of clone i at timepoint j
# (here we have only one timepoint per sample).
# The cumulative frequency is calculated by summing the frequencies of all descendant clones of i.
build_matrix <- function(df_sample) {
  
  trees <- mapping(df_sample)
  if (length(trees) == 0)
    return(matrix(NA, nrow = 0, ncol = 0))
  
  nodes <- map_dfr(trees, extract_nodes_tree)
  desc_map <- map_dfr(trees, get_desc_map)
  
  df <- nodes %>%
    filter(node_id != "root")
  
  mat <- matrix(
    0,
    nrow = nrow(df),
    ncol = 1
  )
  
  rownames(mat) <- df$node_id
  colnames(mat) <- "t1"
  
  freq_map <- setNames(df$freq, df$node_id)
  
  for (i in seq_len(nrow(df))) {
    
    node <- df$node_id[i]
    desc <- desc_map$descendants[
      desc_map$node_id == node
    ][[1]]
    
    vals <- freq_map[desc]
    vals[is.na(vals)] <- 0
    
    mat[i, 1] <- sum(vals)
  }
  
  mat
}

#2. Matrix normalization
# Rescales matrix columns if values exceed 100
normalize_matrix <- function(mat) {
  if (!is.matrix(mat) || ncol(mat) == 0) return(mat)
  
  for (j in seq_len(ncol(mat))) {
    col_values <- mat[, j]
    max_val <- max(col_values, na.rm = TRUE)
    
    if (is.finite(max_val) && max_val > 100) {
      mat[, j] <- col_values * (100 / max_val)
    }
  }
  
  mat
}

#3. Cumulative matrices for all samples
compute_all_cumulative_matrices <- function(clones_df) {
  samples <- unique(clones_df$sample_id)
  
  res <- map(samples, function(s) {
    df_s <- clones_df %>% filter(sample_id == s)
    mat <- build_matrix(df_s)
    
    if (is.matrix(mat) && nrow(mat) > 0) {
      mat <- normalize_matrix(mat)
    }
    
    return(mat)
  })
  
  names(res) <- samples
  res
} #Steps : split data by sample -> build matrix -> normalize -> store in list
