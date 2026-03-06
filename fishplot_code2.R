library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

set.seed(123)

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

#3. JSON loading
json_path <- file.path("data", "trees_aml_morita.json")
data_filename <- tools::file_path_sans_ext(basename(json_path))

#4. Clones_df creation + cleaning
clones_df <- extract_all_samples(my_json) %>%
  mutate(
    parent_id = ifelse(is.na(parent_id), "root", parent_id), # replaces NA parent_id with "root"
    size_percent = suppressWarnings(as.numeric(size_percent)), # converts size_percent to numeric
    mutation = ifelse(is.na(mutation), "none", mutation) #replaces missing mutations with "none"
  )

print(clones_df)

save(clones_df, file = "clones_df.RData")
write.csv(clones_df, "clones_df.csv", row.names = FALSE)

#5. Tree building
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

#6. Mapping function
# For a given sample, identifies root nodes (those with parent_id "root")
# and builds a tree for each root using build_tree. Returns a list of trees for the sample.
mapping <- function(df_sample) {
  roots <- df_sample$node_id[df_sample$parent_id == "root"]
  roots <- as.character(roots)
  
  map(roots, ~ build_tree(df_sample, .x)) %>%
    compact()
} # Result -> a list of tree structures (nested lists) for each root node in the sample

#7. Nodes extraction from tree
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

#8. Descendants collection
#  Returns all descendant node IDs of a node (including itself)
collect_descendants <- function(tree) {
  c(
    tree$node_id,
    unlist(map(tree$children, collect_descendants))
  )
}

#9. Descendants mapping
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

#10. Cumulative matrix construction
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

#11. Matrix normalization
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

#12. Cumulative matrices for all samples
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

all_matrices <- compute_all_cumulative_matrices(clones_df)

# Final check: if any matrix has NA in the first row, replace with 0 (to avoid issues in fishplot)
for (s in names(all_matrices)) {
  mat <- all_matrices[[s]]
  
  if (is.matrix(mat) && nrow(mat) >= 1 && all(is.na(mat[1, ]))) {
    mat[1, ] <- 0
  }
  
  cat("\n====================\n")
  cat("Sample:", s, "\n")
}

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

classified <- classify_samples(clones_df)

cat("Number samples monotime :", length(classified$monotime), "\n")
cat("Number samples multitime :", length(classified$multitime), "\n")

cat("\n--- Samples monotimepoint ---\n")
print(classified$monotime)

cat("\n--- Samples multitimepoints ---\n")
print(classified$multitime)

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
    nodes_s <- df_s$node_id
    parents_s <- df_s$parent_id
    
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
    if (is.null(p) || p == "root") return(0)
    match(p, cumulative_nodes) - 1
  })
  
  root_id <- cn[parent_indices == 0][1]
  keep <- cn != root_id
  
  list(
    final_nodes = cumulative_nodes[keep],
    parent_index = parent_indices[keep]
  )
} # This reconstructs the complete clonal architecture of a patient across time.


patient_samples_map <- split(
  unique(clones_df$sample_id),
  get_patient_id(unique(clones_df$sample_id))
)

all_patient_hierarchies <- list()

for (patient in names(patient_samples_map)) {
  samples <- patient_samples_map[[patient]]
  hier <- build_patient_hierarchy(samples, clones_df)
  all_patient_hierarchies[[patient]] <- hier
  
  cat("\n====================================\n")
  cat("Patient :", patient, "\n")
  cat("Samples :", paste(samples, collapse = ", "), "\n\n")
  cat("Final node_id order:\n")
  print(hier$final_nodes)
  cat("\nParent index vector:\n")
  print(hier$parent_index)
}

# 18. Concatenate cumulative matrices of all timepoints for a patient according to the hierarchy
concat_patient_timepoints_hierarchy <- function(
    patient,
    multitime_by_patient,
    all_matrices,
    all_patient_hierarchies
) {
  
  hier <- all_patient_hierarchies[[patient]]
  if (is.null(hier)) return(NULL)
  
  final_rows <- hier$final_nodes
  patient_samples <- multitime_by_patient[[patient]]
  patient_samples <- patient_samples[order(get_suffix_num(patient_samples))]
  
  # monotimepoint case
  # steps : extracts last column of clones matrix, aligns with final hierarchy and returns a single-column matrix.
  if (length(patient_samples) == 1) {
    
    s <- patient_samples[1]
    mat <- all_matrices[[s]]
    if (is.null(mat)) return(NULL)
    
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
    if (is.null(mat)) return(NULL)
    
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


all_patients <- unique(get_patient_id(clones_df$sample_id))

#19. Apply the concatenation function to all patients and print results
concat_by_patient_hierarchy <- lapply(
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

names(concat_by_patient_hierarchy) <- names(patient_samples_map)

for (p in names(concat_by_patient_hierarchy)) {
  cat("\n====================\n")
  cat("Patient:", p, "\n")
  print(concat_by_patient_hierarchy[[p]])
}

## Numerical corrections

#20. Fix zero reappearance: ensures that once a clone appears (value > 0), it cannot have zero or NA values in subsequent timepoints.
# This is important for fishplot visualization to maintain the continuity of clones over time, even if they become very small.
# steps : cleans numerical issues, detects first appearance and prevents disapperance after appearance
fix_zero_reappearance <- function(mat, eps = 1e-4) {
  
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0 #cleans numerical issues
  
  for (i in seq_len(nrow(mat))) {
    
    row_vals <- mat[i, ]
    
    first_pos <- which(row_vals > 0)[1] #detects first appearance
    
    if (!is.na(first_pos) && first_pos > 1) {
      mat[i, first_pos - 1] <- eps
    }
    
    appeared <- FALSE
    
    for (t in seq_len(ncol(mat))) {
      
      val <- mat[i, t]
      
      if (!is.na(val) && val > 0) {
        appeared <- TRUE
      }
      
      if (appeared && (is.na(val) || val == 0)) { #prevents disapperance after appearance
        mat[i, t] <- eps
      }
    }
  }
  
  mat
} # Clones do not truly disappear between measurements, they may fall below detection but still exists biologically
# Necessary because fishplot requires : continuous trajectories and no abrupt zero-positive transitions

#21. Fix parent-child conflicts: ensures that at any timepoint, the sum of child clone frequencies does not exceed the parent clone frequency.
# Ensures biological consistency
fix_parent_child_conflicts <- function(mat, parents) {
  
  n_time <- ncol(mat)
  n_clones <- nrow(mat)
  
  for (t in seq_len(n_time)) {
    
    for (i in seq_len(n_clones)) {
      
      children <- which(parents == i)
      if (length(children) == 0) next # if children exceed parent -> rescale children to fit parent frequency
      
      parent_val <- mat[i, t]
      sum_children <- sum(mat[children, t])
      
      if (sum_children > parent_val && sum_children > 0) {
        
        ratio <- parent_val / sum_children
        mat[children, t] <- mat[children, t] * ratio
        
      }
    }
  }
  mat
} # A subclone cannot exceed its ancestral clone -> This enforces evolutionary logic

#22. Fix independent clones: ensures that the sum of frequencies of independent clones does not exceed 100% at any timepoint.
# If multiple clones have no parent (parent ==0) -> rescale them proportionally
fix_independent_clones <- function(mat, parents) {
  
  root_idx <- which(parents == 0)
  if (length(root_idx) <= 1) return(mat)
  
  for (t in seq_len(ncol(mat))) {
    s <- sum(mat[root_idx, t])
    
    if (s > 100 && s > 0) {
      mat[root_idx, t] <- mat[root_idx, t] * (100 / s)
    }
  }
  
  mat
} # Total tumor fraction cannot exceed 100%

#23. Enforce all fishplot constraints: applies both parent-child and independent clone corrections iteratively until no violations remain.
# Iteratively enforce all fishplot mathematical constraints
enforce_fishplot_constraints <- function(mat, parents, tol = 1e-8) {
  
  mat <- as.matrix(mat)
  n_time <- ncol(mat)
  n_clones <- nrow(mat)
  
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0
  
  for (t in seq_len(n_time)) {
    
    repeat {
      
      changed <- FALSE
      
      # PARENT / CHILD constraint
      for (i in seq_len(n_clones)) {
        
        children <- which(parents == i)
        if (length(children) == 0) next
        
        parent_val <- mat[i, t]
        sum_children <- sum(mat[children, t])
        
        if (sum_children > parent_val + tol && sum_children > 0) {
          
          ratio <- parent_val / sum_children
          mat[children, t] <- mat[children, t] * ratio
          changed <- TRUE
        }
      }
      
      #  ROOT constraint
      root_idx <- which(parents == 0)
      if (length(root_idx) > 0) {
        
        root_sum <- sum(mat[root_idx, t])
        
        if (root_sum > 100 + tol && root_sum > 0) {
          
          mat[root_idx, t] <- mat[root_idx, t] * (100 / root_sum)
          changed <- TRUE
        }
      }
      
      if (!changed) break
    }
  }
  
  mat
}

# Assigns consitent colors to clones globally : same clone = same color across patients
clone_palette <- c(
  "lightcoral", "skyblue3", "sandybrown", "paleturquoise3",
  "thistle", "darkolivegreen3", "gold", "orchid",
  "steelblue", "salmon", "khaki3", "plum",
  "lightseagreen", "rosybrown", "dodgerblue3"
)

global_clone_palette <- clones_df %>%
  distinct(node_id) %>%
  arrange(as.numeric(node_id)) %>%
  pull(node_id) %>%
  setNames(
    rep(clone_palette, length.out = length(.)),
    .
  )

#24. Creates mapping : grouped and deduplicated
build_mutation_lookup <- function(clones_df) {
  clones_df %>%
    filter(mutation != "none") %>%
    group_by(sample_id, node_id) %>%
    summarise(
      mutation = first(mutation),
      .groups = "drop"
    )
}

mutation_lookup <- build_mutation_lookup(clones_df)

#25. Generates readable legend label
get_clone_label <- function(
    node_id,
    patient_id,
    mutation_lookup,
    max_show = 2
) {
  muts <- mutation_lookup %>%
    filter(grepl(patient_id, sample_id)) %>%
    filter(node_id == !!node_id) %>%
    pull(mutation) %>%
    unique()
  
  if (length(muts) == 0) {
    return(paste0("clone ", node_id))
  }
  
  if (length(muts) <= max_show) {
    return(paste(muts, collapse = ", "))
  }
  
  paste0(
    paste(muts[1:max_show], collapse = ", "),
    " + ",
    length(muts) - max_show,
    " muts"
  )
}

#26. Adds a legend to the fishplot with clone labels
add_fishplot_legend <- function(fish, clone_labels) {
  stopifnot(
    is.character(clone_labels),
    length(clone_labels) == length(fish@col)
  )
  
  par(xpd = TRUE)
  legend(
    "bottom",
    legend = clone_labels,
    fill = fish@col,
    border = "black",
    ncol = 2,
    cex = 0.8,
    bty = "n",
    inset = -0.15
  )
  par(xpd = FALSE)
}

outdir <- "fishplots_png"
if (!dir.exists(outdir)) dir.create(outdir)


all_patients <- names(concat_by_patient_hierarchy)
final_fishplot_matrices <- list()

for (patient_id in all_patients) {
  
  cat("\n===============================\n")
  cat("PATIENT :", patient_id, "\n")
  cat("===============================\n")
  
  mat <- concat_by_patient_hierarchy[[patient_id]]
  
  if (is.null(mat) || all(is.na(mat))) {
    final_fishplot_matrices[[patient_id]] <- NA
    next
  }
  
  parents <- all_patient_hierarchies[[patient_id]]$parent_index
  
  mat_final <- mat |>
    fix_zero_reappearance()
  
  mat_final <- enforce_fishplot_constraints(mat_final, parents)
  
  if (ncol(mat_final) == 1) {
    mat_final <- cbind(mat_final, mat_final)
  }
  
  final_fishplot_matrices[[patient_id]] <- mat_final
  
  clone_ids <- rownames(mat_final)
  fish_colors <- global_clone_palette[clone_ids]
  timepoints <- seq_len(ncol(mat_final))
  
  mat_final[!is.finite(mat_final)] <- 0
  mat_final[mat_final < 0] <- 0
  
  
  # SAFE createFishObject
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
    
    plot.new()
    text(
      0.5, 0.5,
      paste(
        "Fishplot cannot be generated.\n",
        "Drug configuration invalid.\n\n",
        e$message
      ),
      cex = 1
    )
    
    return(NULL)
  })
  
  if (is.null(fishObj)) return()
  
  
  fishObj <- layoutClones(
    fishObj,
    separate.independent.clones = FALSE
  )
  
  
  clone_labels <- sapply(
    clone_ids,
    get_clone_label,
    patient_id = patient_id,
    mutation_lookup = mutation_lookup
  )
  
  clone_labels <- vapply(clone_labels, as.character, character(1))
  
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
  
  cat("Fishplot recorded :", patient_id, "\n")
}

cat("\n=======================================\n")
cat(" ALL FINAL MATRICES USED FOR FISHPLOTS \n")
cat("=======================================\n")

for (p in names(final_fishplot_matrices)) {
  cat("\n--- Patient :", p, "---\n")
  print(final_fishplot_matrices[[p]])
}

final_fishplot_matrices


plot_fishplot <- function(
    patient,
    concat_by_patient_hierarchy,
    all_patient_hierarchies,
    mutation_lookup,
    clone_palette,
    input,
    rv,
    mini = FALSE
) {
  
  if (mini) {
    par(
      mar = c(3, 3, 2, 3) + 0.1,
      xpd = NA
    )
  } else {
    par(
      mar = c(7, 4, 4, 4) + 0.1,
      xpd = NA
    )
  }
  
  
  mat <- concat_by_patient_hierarchy[[patient]]
  
  if (is.null(mat) || all(is.na(mat))) {
    plot.new()
    text(0.5, 0.5, paste("No fishplot for", patient))
    return()
  }
  
  parents <- all_patient_hierarchies[[patient]]$parent_index
  mat_final <- mat
  
  # FIX MONOTIME BEFORE DRUG
  if (ncol(mat_final) == 1) {
    mat_final <- cbind(mat_final, mat_final)
  }
  
  
  # APPLY DRUG EFFECT
  if (ncol(mat_final) >= 3) {
    
    for (col in seq(2, ncol(mat_final) - 1, by = 2)) {
      
      key <- paste(patient, input$selected_drug, col, sep = "_")
      effects <- rv$drug_effects[[key]]
      
      if (is.null(effects)) {
        effects <- generate_random_drug_effects(rownames(mat_final))
        rv$drug_effects[[key]] <- effects
      }
      
      for (i in seq_len(nrow(mat_final))) {
        mat_final[i, col] <- mat_final[i, col - 1] * effects[i]
      }
    }
  }
  
  
  mat_final <- mat_final |>
    fix_zero_reappearance()
  
  mat_final <- enforce_fishplot_constraints(mat_final, parents)
  
  mat_final[!is.finite(mat_final)] <- 0
  mat_final[mat_final < 0] <- 0
  
  
  clone_ids <- rownames(mat_final)
  fish_colors <- clone_palette[clone_ids]
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
    
    plot.new()
    text(
      0.5, 0.5,
      paste("Invalid configuration after drug effect.\n", e$message),
      cex = 1
    )
    
    return(NULL)
  })
  
  if (is.null(fishObj)) return()
  
  
  fishObj <- layoutClones(
    fishObj,
    separate.independent.clones = FALSE
  )
  
  fishPlot(
    fishObj,
    shape = "spline",
    vlines = if (mini) NULL else timepoints,
    vlab = if (mini) NULL else paste0("t", timepoints),
    title.btm = if (mini) "" else patient
  )
  
  if (!mini) {
    clone_labels <- sapply(
      clone_ids,
      get_clone_label,
      patient_id = patient,
      mutation_lookup = mutation_lookup
    )
    
    clone_labels <- as.character(clone_labels)
    add_fishplot_legend(fishObj, clone_labels)
  }
}


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

metadata_df <- extract_metadata(my_json)
print(metadata_df)

save(metadata_df, file = "metadata.RData")

