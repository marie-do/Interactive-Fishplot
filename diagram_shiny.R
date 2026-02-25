
#core dot builder
build_mutation_dot <- function(df, title, samples,
                               global_clone_palette,
                               width_scale = 15) {
  
  clones <- df %>%
    distinct(node_id, parent_id, mutation) %>%
    arrange(as.numeric(node_id))
  
  clone_ids <- clones$node_id
  is_root   <- clones$mutation == "none"
  
  clone_colors <- global_clone_palette[clone_ids]
  
  mutation_nodes <- paste0(
    clone_ids,
    '[style=empty, label="',
    ifelse(is_root, "Root", clones$mutation),
    '",
    URL="javascript:void(0)",
    fontcolor=dimgray,
    penwidth=2];'
  )
  
  mutation_edges <- clones %>%
    filter(parent_id != "root") %>%
    transmute(
      paste0(
        parent_id, " -> ", node_id,
        ' [color=dimgray, weight=4, penwidth=5];'
      )
    ) %>% pull()
  
  time_nodes <- c()
  time_edges <- c()
  
  for (t in seq_along(samples)) {
    
    s <- samples[t]
    
    df_t <- df %>%
      filter(sample_id == s) %>%
      select(node_id, size_percent)
    
    df_t <- tibble(node_id = clone_ids) %>%
      left_join(df_t, by = "node_id") %>%
      mutate(
        size_percent = replace_na(size_percent, 0),
        width = pmax(size_percent / width_scale, 0.01),
        label = ifelse(is_root, "", paste0(round(size_percent, 1), "%"))
      )
    
    time_nodes <- c(
      time_nodes,
      paste0(
        "t", t, "_", df_t$node_id[!is_root],
        '[label="', df_t$label[!is_root],
        '", width=', df_t$width[!is_root],
        ', fixedsize=true, forcelabels=true, shape=circle, style=filled, color=',
        clone_colors[!is_root], '];'
      )
    )
    
    from <- if (t == 1) df_t$node_id else paste0("t", t - 1, "_", df_t$node_id)
    to   <- paste0("t", t, "_", df_t$node_id)
    
    time_edges <- c(
      time_edges,
      paste0(
        from[!is_root], " -> ", to[!is_root],
        '[arrowhead=none, style=dashed, penwidth=5, color=',
        clone_colors[!is_root], '];'
      )
    )
  }
  
  cluster_children <- lapply(which(!is_root), function(i) {
    cid <- clone_ids[i]
    paste0(
      "subgraph cluster__", i, " {\n",
      "label=\"\";\n",
      "rank=same;\n",
      paste0("t", seq_along(samples), "_", cid, ";", collapse = "\n"),
      "\npenwidth=0;\n}"
    )
  })
  
  time_label <- paste0("Timepoints\n", paste(samples, collapse = "      "))
  
  dot <- paste(
    "digraph g {",
    "rankdir=LR;",
    "nodesep=1.5;",
    'labelloc="t";',
    'fontname="Helvetica bold";',
    "fontsize=28;",
    paste0('label="', title, '";'),
    "",
    "node [color=dimgray, fontsize=24, style=filled, fontcolor=black, penwidth=5];",
    "",
    paste(mutation_nodes, collapse = "\n"),
    "",
    paste(mutation_edges, collapse = "\n"),
    "",
    paste(time_nodes, collapse = "\n"),
    "",
    paste(time_edges, collapse = "\n"),
    "",
    "subgraph cluster_timepoints {",
    paste0('label="', time_label, '";'),
    "rank=same;",
    "penwidth=0;",
    paste(cluster_children, collapse = "\n"),
    "}",
    "}"
  )
  
  grViz(dot)
}


#Shiny wraper
plot_mutation_tree_shiny <- function(patient_id, clones_df,
                                     global_clone_palette) {
  
  samples <- clones_df %>%
    filter(str_starts(sample_id, paste0(patient_id, "-"))) %>%
    distinct(sample_id) %>%
    pull(sample_id) %>%
    sort()
  
  if (length(samples) == 0) return(NULL)
  
  df_tree <- clones_df %>%
    filter(sample_id %in% samples)
  
  title <- if (length(samples) == 1) samples else patient_id
  
  build_mutation_dot(
    df      = df_tree,
    title   = title,
    samples = samples,
    global_clone_palette = global_clone_palette
  )
}

