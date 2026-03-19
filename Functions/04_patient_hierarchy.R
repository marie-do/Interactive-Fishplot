#14. Sample classification
get_patient_id <- function(sample_id) {
  sub("^(.*)-[^-]+$", "\\1", sample_id)
}

#15. Classify samples into monotime and multitime based on patient_id
classify_samples <- function(clones_df) {
  sample_ids <- unique(clones_df$sample_id)
  patient_ids <- get_patient_id(sample_ids)
  
  df <- tibble(
    sample_id = sample_ids,
    patient_id = patient_ids
  )
  
  counts <- df %>% count(patient_id, name = "n_samples") # computes how many samples each patient has
  
  patients_monotime <- counts %>%
    filter(n_samples == 1) %>%
    pull(patient_id)
  
  patients_multitime <- counts %>%
    filter(n_samples > 1) %>%
    pull(patient_id)
  
  list(
    monotime = sort(
      df %>%
        filter(patient_id %in% patients_monotime) %>%
        pull(sample_id)
    ),
    multitime = sort(
      df %>%
        filter(patient_id %in% patients_multitime) %>%
        pull(sample_id)
    ),
    multitime_by_patient = df %>%
      filter(patient_id %in% patients_multitime) %>%
      group_by(patient_id) %>%
      summarize(
        samples = list(sort(sample_id)),
        .groups = "drop"
      ) %>%
      deframe()
  )
} # This distinguishes patient observed once and Patients with longitudinal follow-up

#16. get_suffix_num function : extract the timepoint number from a sample ID
# (e.g., from "AML-02-001" it would extract 1, from "AML-02-002" it would extract 2, etc.)
get_suffix_num <- function(x) {
  as.numeric(sub(".*-", "", x))
}

#17. Build a cumulative clonal hierarchy across all timepoints of a patient
build_patient_hierarchy <- function(patient_samples, clones_df) {
  patient_samples <- patient_samples[
    order(get_suffix_num(patient_samples)) # ensures chronological reconstruction
  ]
  
  cumulative_nodes <- c() #stores all unique clones observed
  parent_map <- list() # stores parent relationship
  
  for (s in patient_samples) {
    df_s <- clones_df %>% filter(sample_id == s)
    nodes_s <- as.character(df_s$node_id)
    parents_s <- as.character(df_s$parent_id)
    
    new_nodes <- nodes_s[!nodes_s %in% cumulative_nodes]
    cumulative_nodes <- c(cumulative_nodes, new_nodes)
    
    for (nd in new_nodes) {
      parent_map[[as.character(nd)]] <-
        parents_s[df_s$node_id == nd]
    }
  }# for each sample, it extracts nodes and parents, identify new clones not see before and adds them to cumulative list -> This reconstructs evolutionary accumulation
  
  
  cn <- as.character(cumulative_nodes)
  
  parent_indices <- sapply(cn, function(nd) {
    
    p <- parent_map[[nd]]
    
    if (is.null(p) || p == "root") {
      return(0)
    }
    
    idx <- match(p, cumulative_nodes)
    
    if (is.na(idx)) {
      return(0)
    }
    
    idx - 1
  })
  
  # remove artificial root node (the one whose mutation == "root")
  root_nodes <- names(parent_map)[unlist(parent_map) == "root"]
  
  keep <- !(cn %in% root_nodes)
  
  list(
    final_nodes = cumulative_nodes[keep],
    parent_index = parent_indices[keep]
  )
} # This reconstructs the complete clonal architecture of a patient across time.

# 18. Concatenate cumulative matrices of all timepoints for a patient according to the hierarchy
concat_patient_timepoints_hierarchy <- function(
    patient,
    multitime_by_patient,
    all_matrices,
    all_patient_hierarchies
) {
  
  hier <- all_patient_hierarchies[[patient]]
  if (is.null(hier)) return(NULL)
  
  final_rows <- as.character(hier$final_nodes)
  patient_samples <- multitime_by_patient[[patient]]
  patient_samples <- patient_samples[order(get_suffix_num(patient_samples))]
  
  # monotimepoint case
  # steps : extracts last column of clones matrix, aligns with final hierarchy and returns a single-column matrix.
  if (length(patient_samples) == 1) {
    
    s <- patient_samples[1]
    mat <- all_matrices[[s]]
    
    if (is.null(mat)) {
      message("Matrix NULL for sample: ", s)
    }
    
    if (!is.null(mat) && ncol(mat) == 0) {
      message("Matrix with 0 columns for sample: ", s)
    }
    
    if (is.null(mat) || ncol(mat) == 0 || nrow(mat) == 0) {
      return(NULL)
    }
    
    if (is.null(mat) || ncol(mat) == 0 || nrow(mat) == 0) {
      return(NULL)
    }
    
    last <- mat[, ncol(mat), drop = FALSE]
    
    t1 <- matrix(
      0,
      nrow = length(final_rows),
      ncol = 1,
      dimnames = list(final_rows, "t1")
    )
    
    common <- intersect(rownames(last), final_rows)
    t1[common, 1] <- last[common, 1]
    
    return(t1)
  }
  
  #multitimepoint case
  #steps : extracts last column for each timepoint, aligns clones, builds final time matrix with real timepoints and drug timepoints interleaved.
  # This produces a longitudinal clonal evolution matrix per patient.
  
  lastcols <- lapply(patient_samples, function(s) {
    
    mat <- all_matrices[[s]]
    
    if (is.null(mat)) {
      message("Matrix NULL for sample: ", s)
    }
    
    if (!is.null(mat) && ncol(mat) == 0) {
      message("Matrix with 0 columns for sample: ", s)
    }
    
    if (is.null(mat) || ncol(mat) == 0 || nrow(mat) == 0) {
      return(NULL)
    }
    
    last <- mat[, ncol(mat), drop = FALSE]
    
    aligned <- matrix(
      0,
      nrow = length(final_rows),
      ncol = 1,
      dimnames = list(final_rows, NULL)
    )
    
    common <- intersect(rownames(last), final_rows)
    aligned[common, 1] <- last[common, 1]
    
    aligned
  })
  
  lastcols <- lastcols[!sapply(lastcols, is.null)]
  
  final_mat <- matrix(
    0,
    nrow = length(final_rows),
    ncol = 0,
    dimnames = list(final_rows, NULL)
  )
  
  col_idx <- 1
  
  for (i in seq_along(lastcols)) {
    
    real_tp <- lastcols[[i]]
    
    final_mat <- cbind(final_mat, real_tp)
    colnames(final_mat)[ncol(final_mat)] <- paste0("t", col_idx)
    col_idx <- col_idx + 1
    
    if (i < length(lastcols)) {
      
      drug_tp <- real_tp
      
      final_mat <- cbind(final_mat, drug_tp)
      colnames(final_mat)[ncol(final_mat)] <- paste0("t", col_idx)
      col_idx <- col_idx + 1
    }
  }
  
  final_mat
}