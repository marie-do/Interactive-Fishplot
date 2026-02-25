# Assigns consistent colors to clones globally : same clone = same color across patients
clone_palette <- c(
  "lightcoral", "skyblue3", "sandybrown", "paleturquoise3",
  "thistle", "darkolivegreen3", "gold", "orchid",
  "steelblue", "salmon", "khaki3", "plum",
  "lightseagreen", "rosybrown", "dodgerblue3"
)

# 1. Creates a global clone palette based on the unique node_ids across all patients
# (NOW wrapped in a function to avoid using clones_df at source time)
build_global_clone_palette <- function(clones_df) {
  
  node_ids <- clones_df %>%
    distinct(node_id) %>%
    arrange(as.numeric(node_id)) %>%
    pull(node_id)
  
  setNames(
    rep(clone_palette, length.out = length(node_ids)),
    node_ids
  )
}

# 2. Generates readable legend label
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

# 3. Adds a legend to the fishplot with clone labels
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