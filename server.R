#server

server <- function(input, output, session) {
  
  #Reactive Matrix : Construction of the cummulative matrix with all the rules applied (drug effects, zero reappearance, parent-child constraints, root constraint)
  safe_hierarchy_rebalance <- function(mat, parents) {
    
    for (col in seq_len(ncol(mat))) {
      
      root_index <- which(parents == 0)
      mat[root_index, col] <- 100
      
      repeat {
        
        changed <- FALSE
        
        for (i in seq_along(parents)) {
          
          parent <- parents[i]
          if (parent == 0) next
          
          children <- which(parents == parent)
          children_sum <- sum(mat[children, col])
          
          if (children_sum > mat[parent, col] && children_sum > 0) {
            
            scale_factor <- mat[parent, col] / children_sum
            mat[children, col] <- mat[children, col] * scale_factor
            
            changed <- TRUE
          }
        }
        
        if (!changed) break
      }
      
    }
    
    mat
  }
  
  fishplot_matrix_reactive <- reactive({
    
    req(input$patient, rv$objects)
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    if (is.null(mat)) return(NULL)
    
    parents <- rv$objects$all_patient_hierarchies[[input$patient]]$parent_index
    mat_final <- mat
    
    # Fix monotimepoints with single column by duplicating it (fishplot requires at least 2 timepoints)
    if (ncol(mat_final) == 1) {
      
      original_name <- colnames(mat_final)[1]
      
      mat_final <- cbind(mat_final, mat_final)
      
      colnames(mat_final) <- c(
        original_name,
        paste0(original_name, "_duplicated")
      )
    }
    
    # Apply drug effects
    # Only apply drug effects if REAL multiple timepoints exist
    real_timepoints <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient) %>%
      pull(sample_id) %>%
      unique()
    
    if (length(real_timepoints) >= 2) {
      
      for (col in seq(2, ncol(mat_final) - 1, by = 2)) { #drug event columns every 2 columns, so that each timepoint are separated by a drug event.
        key <- paste(input$patient,
                     input$selected_drug,
                     col,
                     sep = "_")
        
        effect <- rv$drug_effects[[key]] # 2 modes : global or per-mutation. 
        #The observeEvent for drug effect editing will update this reactive value with the new drug effect parameters, which will then be applied here to the matrix.
        
        if (!is.null(effect)) { 
          
          if (effect$mode == "global") {
            
            mat_final[, col] <- mat_final[, col - 1] * effect$global_value # -> global mode : the same multiplier is applied to all mutations
            
          } else {
            
            for (i in seq_len(nrow(mat_final))) {
              
              node_id <- rownames(mat_final)[i]
              if (!node_id %in% names(effect$per_mutation)) {
                mult <- 1
              } else {
                mult <- effect$per_mutation[[node_id]]
              }
              mat_final[i, col] <- mat_final[i, col - 1] * mult # -> per-mutation mode : each mutation has its own multiplier
            }
          }
        }
      }
    }
    
    mat_final <- fix_zero_reappearance(mat_final)
    mat_final <- safe_hierarchy_rebalance(mat_final, parents)
    
    mat_final[mat_final < 0] <- 0
    mat_final[!is.finite(mat_final)] <- 0
    
    return(mat_final)
  })
  
  # Helper function to get depth-first order of nodes based on parent-child relationships by recursively traversing the tree structure.
  get_depth_first_order <- function(node_ids, parent_index) {
    
    children_map <- split(seq_along(parent_index), parent_index) # list where names are parent indices and values are vectors of child indices
    
    traverse <- function(idx) {
      
      current <- node_ids[idx] # start with the current node ID
      child_idx <- children_map[[as.character(idx)]] # get child indices of the current node
      
      if (is.null(child_idx)) { # if no children, return current node ID -> leaf node
        return(current)
      }
      
      # depth-first order = biological order
      c(
        current,
        unlist(lapply(child_idx, traverse)) # children are traversed recursively in the same way
      )
    }
    
    # roots = nodes without parent
    roots <- which(parent_index == 0)
    unlist(lapply(roots, traverse))
  }
  
  # Helper function to reindex node IDs after full deletion of a mutation, to ensure that node IDs remain consecutive and consistent.
  # This is important for the fishplot construction and to avoid issues with missing node IDs.
  reindex_nodes <- function(df) {
    
    unique_nodes <- df %>%
      distinct(node_id) %>%
      arrange(as.numeric(node_id)) %>%
      pull(node_id)
    
    new_ids <- as.character(seq_along(unique_nodes))
    id_map <- setNames(new_ids, unique_nodes)
    
    df <- df %>%
      mutate(
        node_id = id_map[node_id],
        parent_id = case_when(
          parent_id == "root" ~ "root",
          parent_id %in% names(id_map) ~ id_map[parent_id],
          TRUE ~ "root" 
        )
      )
    
    return(df)
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
  
  # Initialize available drugs (can be extended later with dynamic addition)
  rv$available_drugs <- c(
    "Midostaurin",
    "Venetoclax",
    "Cytarabine",
    "Azacitidine"
  )
  
  # Drug observer to update drug selection choices
  observe({
    if (!is.null(rv$available_drugs) && length(rv$available_drugs) > 0) {
      updateSelectInput(
        session,
        "selected_drug",
        choices = rv$available_drugs,
        selected = rv$available_drugs[1]
      )
    }
  })
  
  # load file server
  observeEvent(input$load_file, {
    
    req(input$load_file)
    
    my_json <- jsonlite::fromJSON(
      input$load_file$datapath,
      simplifyVector = FALSE
    )
    
    # Clone extraction
    base_clones_df <- extract_all_samples(my_json) %>%
      mutate(
        node_id = as.character(node_id),
        parent_id = ifelse(is.na(parent_id), "root", as.character(parent_id)),
        size_percent = as.numeric(size_percent) * 100,
        mutation = ifelse(is.na(mutation), "none", as.character(mutation)),
        sample_id = as.character(sample_id)
      )
    
    #Metadata extraction 
    base_metadata_df <- extract_metadata(my_json)
    
    # Store in reactive values
    rv$metadata_df <- base_metadata_df
    rv$clones_df <- base_clones_df
    rv$objects <- build_all_objects(base_clones_df)
    
    rv$drug_effects <- list()   
    
    for (patient in my_json) {
      
      if (!is.null(patient$drug_effects)) {
        
        rv$drug_effects <- c(
          rv$drug_effects,
          patient$drug_effects
        )
      }
    }
    
    updateSelectInput(
      session,
      "patient",
      choices = sort(unique(get_patient_id(base_clones_df$sample_id))),
      selected = unique(get_patient_id(base_clones_df$sample_id))[1]
    )
    
    showNotification("File loaded successfully", type = "message")
  })
  
  # Update timepoint list
  observeEvent(input$patient, {
    
    req(rv$clones_df)
    
    tp <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient) %>%
      pull(sample_id) %>%
      unique() %>%
      sort()
    
    updateSelectInput(
      session,
      "selected_timepoint",
      choices = tp,
      selected = tp[1]
    )
  })
  
  # Adding new timepoint
  observeEvent(input$new_timepoint, {
    
    req(rv$clones_df) # Security : only allow adding timepoints if data is loaded
    
    existing <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient) %>%
      pull(sample_id) %>%
      unique()
    
    nums <- as.numeric(str_extract(existing, "\\d+$"))
    nums <- nums[!is.na(nums)]
    new_idx <- ifelse(length(nums) == 0, 1, max(nums) + 1) # Find the highest existing index and add 1, or start at 1 if none exist
    new_sample <- paste0(input$patient, "-", sprintf("%03d", new_idx))
    
    base_rows <- rv$clones_df %>%
      filter(sample_id == existing[1]) %>%
      mutate(
        sample_id = new_sample,
        size_percent = 0
      ) # Copied from the first existing timepoint of the patient, but with 0 percentages
    
    rv$clones_df <- bind_rows(rv$clones_df, base_rows)
    rv$objects <- build_all_objects(rv$clones_df)
    rand_val <- generate_random_global_effect()
    
    #Initialize drug effects for new drug events, so that 
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    
    if (!is.null(mat) && ncol(mat) >= 3) {
      
      drug_cols <- seq(2, ncol(mat), by = 2)
      node_ids <- rownames(mat)
      
      for (drug in rv$available_drugs) {
        
        for (col in drug_cols) {
          
          key <- paste(input$patient, drug, col, sep = "_")
          
          if (is.null(rv$drug_effects[[key]])) {
            
            rand_val <- generate_random_global_effect()
            
            rv$drug_effects[[key]] <- list( #global mode + per mutation mode are ready to use. Initialized with random values, but can be updated if the user edits drug effects.
              mode = "global",
              global_value = rand_val,
              per_mutation = setNames(
                rep(rand_val, length(node_ids)),
                node_ids
              )
            )
          }
        }
      }
    }
    
    updateSelectInput(
      session,
      "selected_timepoint",
      choices = c(existing, new_sample),
      selected = new_sample
    )
  })
  
  # Edit all observer
  observeEvent(input$edit_all, {
    
    req(rv$clones_df, input$selected_timepoint, input$patient)
    
    patient_id <- input$patient
    
    hier <- rv$objects$all_patient_hierarchies[[patient_id]]
    parents <- hier$parent_index
    node_ids <- hier$final_nodes
    
    ordered_nodes <- get_depth_first_order(node_ids, parents)
    
    # Build a dataframe with all nodes for the patient, their mutations, and their current percentages at the selected timepoint (or 0 if not present), ordered by the depth-first order of the tree structure.
    df_structure <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == patient_id) %>%
      distinct(node_id, mutation) %>%
      mutate(node_id = as.character(node_id))
    
    df_tp_values <- rv$clones_df %>%
      filter(sample_id == input$selected_timepoint) %>%
      select(node_id, size_percent) %>%
      mutate(node_id = as.character(node_id))
    
    df_tp <- df_structure %>%
      left_join(df_tp_values, by = "node_id") %>%
      mutate(
        size_percent = ifelse(is.na(size_percent), 0, size_percent),
        node_id = factor(node_id, levels = ordered_nodes)
      ) %>%
      arrange(node_id)
    
    
    node_info <- get_node_labels(input$patient, rv$clones_df)
    
    showModal(
      modalDialog(
        title = paste("Edit all percentages –", input$selected_timepoint),
        size = "l",
        easyClose = FALSE,
        fade = TRUE,
        
        # top section with mutation tree
        h4("Mutation tree"),
        div(
          style = "
        height: 500px;
        overflow: auto;
        border: 1px solid #ddd;
        margin-bottom: 20px;
      ",
          grVizOutput("mini_tree", height = "480px")
        ),
        
        hr(),
        
        fluidRow(
          
          # left column with numeric inputs for percentages (0-100)
          column(
            7,
            h4("Edit clone percentages"),
            div(
              style = "
            max-height: 400px;
            overflow-y: auto;
            padding-right: 10px;
          ",
              tagList(
                lapply(seq_len(nrow(df_tp)), function(i) {
                  
                  current_node <- df_tp$node_id[i]
                  
                  current_label <- node_info$label[
                    match(current_node, node_info$node_id)
                  ]
                  
                  numericInput(
                    inputId = paste0("bulk_", current_node),
                    label   = current_label,
                    value   = df_tp$size_percent[i],
                    min     = 0,
                    max     = 100,
                    step    = 0.1
                  )
                })
              )
            )
          ),
          
          # Right column with mutation summary table
          column(
            5,
            h4("Mutation summary"),
            DTOutput("mini_table")
          )
        ),
        
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "confirm_bulk_edit",
            "Apply changes",
            class = "btn-success"
          )
        )
      )
    )
  })
  
  # Confirm bulk edit observer
  observeEvent(input$confirm_bulk_edit, {
    
    df_tp <- rv$clones_df %>%
      filter(sample_id == input$selected_timepoint) # retrieve the clones of the current timepoint to get the list of nodes to update, and their current percentages (to fill in the numeric inputs)
    
    for (node in df_tp$node_id) {
      
      input_id <- paste0("bulk_", node)
      
      if (!is.null(input[[input_id]])) {
        rv$clones_df <- rv$clones_df %>%
          mutate(
            size_percent = ifelse(
              sample_id == input$selected_timepoint &
                node_id == node,
              input[[input_id]],
              size_percent
            )
          )
      }
    }
    
    # After editing the percentages, the normalization function is applied to ensure that the total percentage of all clones at the current timepoint does not exceed 100%.
    rv$clones_df <- normalize_timepoint_percentages(rv$clones_df)
    
    rv$objects <- build_all_objects(rv$clones_df)
    removeModal()
  })
  
  # New mutation observer
  observeEvent(input$add_mutation, {
    
    req(rv$clones_df)
    
    df_patient <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient)
    
    df_tree <- df_patient %>%
      distinct(node_id, parent_id, mutation) %>%
      arrange(as.numeric(node_id))
    
    showModal(
      modalDialog(
        title = paste("Add new mutation – Patient", input$patient),
        size = "l",
        
        fluidRow(
          column(
            6,
            h4("Current mutation tree"),
            grVizOutput("mini_tree", height = "400px")
          ),
          column(
            6,
            h4("Mutation summary"),
            DTOutput("mini_table")
          )
        ),
        
        hr(),
        h4("New mutation parameters"),
        
        selectInput(
          "new_parent",
          "Select parent node",
          choices = df_tree$node_id
        ),
        
        textInput(
          "new_mutation_name",
          "Mutation name",
          placeholder = "e.g. TP53 R175H"
        ),
        
        numericInput(
          "new_mutation_percent",
          "Percentage (current timepoint)",
          value = 0,
          min = 0,
          max = 100,
          step = 0.1
        ),
        
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "confirm_add_mutation",
            "Add mutation",
            class = "btn-success"
          )
        )
      )
    )
  })
  
  # Confirm add mutation observer
  observeEvent(input$confirm_add_mutation, {
    
    req(rv$clones_df)
    removeModal()
    
    # new node ID generation (incremental based on the maximum existing node ID to ensure uniqueness)
    new_id <- max(as.numeric(rv$clones_df$node_id)) + 1
    new_id <- as.character(new_id)
    
    # New row creation
    if (input$new_mutation_name == "") {
      showNotification("Mutation name cannot be empty", type = "error")
      return()
    }
    
    new_row <- rv$clones_df %>%
      filter(sample_id == input$selected_timepoint) %>%
      slice(1) %>%
      mutate(
        node_id = new_id,
        parent_id = input$new_parent,
        mutation = input$new_mutation_name,
        size_percent = input$new_mutation_percent
      )
    
    # Error handling
    if (input$new_parent == new_id) {
      showNotification("Parent cannot be itself", type = "error")
      return()
    }
    
    patient_samples <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient) %>%
      pull(sample_id) %>%
      unique()
    
    # When a new mutation is added, at the selected timepoint it will have the user-defined percentage, but at all other timepoints it will start with 0%. 
    all_new_rows <- lapply(patient_samples, function(s) {
      new_row %>%
        mutate(
          sample_id = s,
          size_percent = ifelse(
            s == input$selected_timepoint,
            input$new_mutation_percent,
            0
          )
        )
    }) %>%
      bind_rows()
    
    rv$clones_df <- bind_rows(rv$clones_df, all_new_rows)
    # Normalization + objects rebuilding to ensure the new mutation is properly integrated. 
    rv$clones_df <- normalize_timepoint_percentages(rv$clones_df)
    rv$objects <- build_all_objects(rv$clones_df)
  })
  
  # when a button is clicked th UI is updated to display the corresponding view.
  observeEvent(input$show_fish, { rv$view <- "fish" })
  observeEvent(input$show_tree, { rv$view <- "tree" })
  observeEvent(input$show_data, { rv$view <- "data" })
  
  observeEvent(input$sidebar_menu, {
    
    if (input$sidebar_menu == "viz") {
      rv$view <- "fish"
    }
    
    if (input$sidebar_menu == "format") {
      rv$view <- "format"
    }
  })
  
  # Global display area + parameters 
  output$main_view <- renderUI({
    
    if (rv$view == "fish") {
      
      tagList(
        div(
          style = "margin-bottom:10px;",
          actionButton(
            "rename_events",
            "Rename events",
            icon = icon("edit"),
            class = "btn-primary"
          )
        ),
        plotOutput("fishplot", height = "700px")
      )
    } else if (rv$view == "tree") {
      
      tagList(
        div(
          style = "margin-bottom:10px;",
          actionButton("zoom_in", "+", class = "btn-primary"),
          actionButton("zoom_out", "-", class = "btn-primary"),
          actionButton("zoom_reset", "Reset", class = "btn-default")
        ),
        div(
          style = "overflow:auto; border:1px solid #ddd;",
          grVizOutput("mutation_tree", height = "800px")
        )
      )
      
    } else if (rv$view == "format") {
      box(width = 12, includeMarkdown("data_format.md"))
      
    } else if (rv$view == "metadata") {
      tagList(
        h3(paste("Metadata – Patient", input$patient)),
        DTOutput("metadata_table"),
        br(),
        div(
          style = "display:flex; gap:10px;",
          actionButton("add_metadata_field", "Add metadata field", class = "btn-primary"),
          actionButton("delete_metadata_field", "Delete metadata field", class = "btn-danger")
        )
      )
      
    } else if (rv$view == "matrix") {
      tagList(
        h3(paste("Fishplot matrix - Patient", input$patient)),
        br(),
        p(strong("Matrix rules applied:")),
        tags$ul(
          tags$li("Drug effect applied"),
          tags$li("Zero reappearance correction"),
          tags$li("Parent-child constraints enforced"),
          tags$li("Root constraint enforced")
        ),
        br(),
        DTOutput("matrix_table")
      )
      
    } else {
      DTOutput("table")
      
    }
  })
  
  # Mini mutation tree displayed in the bulk edit modal to help better understand the tree structure while editing percentages.
  output$mini_tree <- renderGrViz({
    req(rv$clones_df, rv$objects)
    plot_mutation_tree_shiny(
      input$patient,
      rv$clones_df,
      rv$objects$clone_palette
    )
  })
  
  # Mini table displayed in the bulk edit modal to summarize the mutations and their current percentages at the selected timepoint.
  output$table <- renderDT({
    
    req(rv$clones_df, rv$objects, input$selected_timepoint)
    
    patient_id <- input$patient
    
    mat <- rv$objects$concat_by_patient_hierarchy[[patient_id]]
    parents <- rv$objects$all_patient_hierarchies[[patient_id]]$parent_index
    node_ids <- rownames(mat)
    
    # Biological ordering
    ordered_nodes <- get_depth_first_order(node_ids, parents)
    
    df_display <- rv$clones_df %>%
      filter(sample_id == input$selected_timepoint) %>%
      mutate(node_id = as.character(node_id)) %>%
      filter(node_id %in% ordered_nodes) %>%
      mutate(node_id = factor(node_id, levels = ordered_nodes)) %>%
      arrange(node_id)
    
    root_row <- data.frame(
      node_id = ordered_nodes[1],
      parent_id = NA,
      size_percent = NA,
      mutation = "Root",
      sample_id = input$selected_timepoint
    )
    
    df_display$size_percent <- round(df_display$size_percent, 2)
    
    datatable(df_display, rownames = FALSE)
  })
  
  # Matrix table to display the reactive matrix with all rules applied, for better understanding of how the final fishplot is constructed.
  output$fishplot <- renderPlot({
    
    req(input$patient)
    
    mat_corrected <- fishplot_matrix_reactive()
    req(mat_corrected)
    
    custom_labels <- NULL
    
    if (!is.null(rv$event_labels[[input$patient]])) {
      custom_labels <- rv$event_labels[[input$patient]]
    }
    
    corrected_list <- list()
    corrected_list[[input$patient]] <- mat_corrected
    
    plot_fishplot(
      patient = input$patient,
      concat_by_patient_hierarchy = corrected_list,
      all_patient_hierarchies = rv$objects$all_patient_hierarchies,
      mutation_lookup = rv$objects$mutation_lookup,
      clone_palette = rv$objects$clone_palette,
      input = input,
      rv = rv,
      mini = FALSE,
      event_labels = custom_labels
    )
  })
  
  # Mutation tree with zoom functionality
  output$mutation_tree <- renderGrViz({
    
    req(rv$clones_df, rv$objects)
    
    zoom <- rv$zoom_level 
    
    plot_mutation_tree_shiny(
      input$patient,
      rv$clones_df,
      rv$objects$clone_palette
    ) %>%
      htmlwidgets::onRender(sprintf("
      function(el) {

        var zoom = %f;

        var svg = el.querySelector('svg');

        if(svg){
          svg.style.transformOrigin = '0 0';
          svg.style.transform = 'scale(' + zoom + ')';
        }
      }
    ", zoom))
  })
  
  # Save file server
  output$save_file <- downloadHandler(
    filename = function() {
      paste0("modified_fishplot_data_", Sys.Date(), ".json")
    },
    
    content = function(file) {
      
      req(rv$clones_df)
      
      df <- rv$clones_df %>%
        mutate(
          patient_id = get_patient_id(sample_id),
          size_percent = size_percent / 100
        )
      
      meta_df <- rv$metadata_df
      
      if (is.null(meta_df)) {
        meta_df <- data.frame(sample_id = unique(df$sample_id))
      }
      
      json_structure <- df %>%
        group_split(patient_id) %>%
        lapply(function(patient_df) {
          
          patient_id_val <- unique(patient_df$patient_id)
          
          patient_drug_effects <- rv$drug_effects[
            startsWith(
              names(rv$drug_effects),
              paste0(patient_id_val, "_")
            )
          ]
          
          if (is.null(patient_drug_effects)) {
            patient_drug_effects <- list()
          }
          
          samples_list <- patient_df %>%
            group_split(sample_id) %>%
            lapply(function(sample_df) {
              
              sample_id_val <- unique(sample_df$sample_id)
              
              meta_row <- meta_df %>%
                filter(sample_id == sample_id_val) %>%
                select(-sample_id)
              
              list(
                sample_id = sample_id_val,
                
                metadata = if (nrow(meta_row) > 0) {
                  as.list(meta_row[1, ])
                } else {
                  list()
                },
                
                clones = sample_df %>%
                  arrange(as.numeric(node_id)) %>%
                  select(node_id, parent_id, mutation, size_percent) %>%
                  mutate(node_id = as.numeric(node_id)) %>%
                  split(seq_len(nrow(.))) %>%
                  lapply(as.list)
              )
            })
          
          list(
            patient_id = patient_id_val,
            samples = samples_list,
            drug_effects = patient_drug_effects  
          )
        })
      
      jsonlite::write_json(
        json_structure,
        file,
        pretty = TRUE,
        auto_unbox = TRUE,
        null = "null"
      )
    }
  )
  
  #Observer for deleting timepoint
  observeEvent(input$delete_timepoint, {
    
    req(rv$clones_df, input$selected_timepoint)
    
    showModal(
      modalDialog(
        title = "Confirm deletion",
        paste(
          "Are you sure you want to delete timepoint",
          input$selected_timepoint,
          "?"
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "confirm_delete_timepoint",
            "Delete",
            class = "btn-danger"
          )
        )
      )
    )
  })
  
  
  # observer for confirming timepoint deletion
  observeEvent(input$confirm_delete_timepoint, {
    
    req(rv$clones_df, input$selected_timepoint, input$patient)
    
    tp_to_delete <- input$selected_timepoint
    patient <- input$patient
    
    rv$clones_df <- rv$clones_df %>%
      filter(sample_id != tp_to_delete)
    
    rv$clones_df <- normalize_timepoint_percentages(rv$clones_df)
    rv$objects <- build_all_objects(rv$clones_df)
    
    remaining_timepoints <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == patient) %>%
      pull(sample_id) %>%
      unique() %>%
      sort()
    
    updateSelectInput(
      session,
      "selected_timepoint",
      choices = remaining_timepoints,
      selected = if (length(remaining_timepoints) > 0)
        remaining_timepoints[1]
      else
        NULL
    )
    
    removeModal()
    
    showNotification("Timepoint deleted", type = "warning")
  })
  
  # Deletion helper text
  output$deletion_preview <- renderText({
    
    req(input$mutation_to_delete, input$selected_timepoint)
    
    node <- input$mutation_to_delete
    patient <- input$patient
    current_tp <- input$selected_timepoint
    
    first_appearance <- rv$clones_df %>%
      filter(
        node_id == node,
        get_patient_id(sample_id) == patient,
        size_percent > 0
      ) %>%
      pull(sample_id) %>%
      sort() %>%
      .[1]
    
    if (!is.na(first_appearance) && first_appearance == current_tp) {
      
      paste(
        "This mutation first appears at the current timepoint.",
        "\n→ It will be completely removed from the dataset.",
        "\n→ It will disappear from all timepoints.",
        "\n→ The evolutionary structure will be updated."
      )
      
    } else {
      
      paste(
        "This mutation existed before the current timepoint.",
        "\n→ It cannot biologically disappear.",
        "\n→ It will be reduced to near-zero from this timepoint onward.",
        "\n→ Previous timepoints will remain unchanged."
      )
    }
  })
  
  
  #Observer for mutation deletion
  observeEvent(input$delete_mutation, {
    
    req(rv$clones_df, input$patient)
    
    df_patient <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient)
    
    df_structure <- df_patient %>%
      distinct(node_id, mutation, parent_id)
    
    nodes_with_children <- df_structure$parent_id
    
    # Only mutations that are active at the current timepoint and do not have children can be deleted, to ensure biological consistency of the tree structure.
    current_active_nodes <- df_patient %>%
      filter(sample_id == input$selected_timepoint,
             size_percent > 0) %>%
      pull(node_id)
    
    deletable_clones <- df_structure %>%
      filter(
        parent_id != "root",
        !node_id %in% nodes_with_children,
        node_id %in% current_active_nodes
      ) %>%
      arrange(as.numeric(node_id))
    
    if (nrow(deletable_clones) == 0) {
      showNotification(
        "No deletable mutations available.\nOnly active leaf mutations (without children) can be removed.",
        type = "warning"
      )
      return()
    }
    
    showModal(
      modalDialog(
        title = "Delete mutation",
        size = "l",
        
        h4("Deletion rules"),
        tags$ul(
          tags$li("Only leaf mutations can be deleted."),
          tags$li("Root and internal nodes are protected."),
          tags$li("Historical mutations cannot fully disappear."),
          tags$li("New mutations (first appearing now) can be completely removed.")
        ),
        
        hr(),
        
        selectInput(
          "mutation_to_delete",
          "Select mutation",
          choices = setNames(
            deletable_clones$node_id,
            paste0("Node ", deletable_clones$node_id,
                   " | ", deletable_clones$mutation)
          )
        ),
        
        hr(),
        
        h4("What will happen?"),
        verbatimTextOutput("deletion_preview"),
        
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "confirm_delete_mutation",
            "Confirm deletion",
            class = "btn-danger"
          )
        )
      )
    )
  })
  
  # Observer for confirming mutation deletion with different scenarios based on whether the mutation first appears at the current timepoint or existed before.
  observeEvent(input$confirm_delete_mutation, {
    
    req(rv$clones_df, input$mutation_to_delete, input$selected_timepoint)
    
    node_to_delete <- input$mutation_to_delete
    patient <- input$patient
    current_tp <- input$selected_timepoint
    
    patient_timepoints <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == patient) %>%
      pull(sample_id) %>%
      unique() %>%
      sort()
    
    first_appearance <- rv$clones_df %>%
      filter(
        node_id == node_to_delete,
        get_patient_id(sample_id) == patient,
        size_percent > 0
      ) %>%
      pull(sample_id) %>%
      sort() %>%
      .[1]
    
    if (!is.na(first_appearance) && first_appearance == current_tp) {
      
      rv$clones_df <- rv$clones_df %>%
        filter(node_id != node_to_delete)
      
      # Clean the drug effect associated
      rv$drug_effects <- rv$drug_effects[
        !grepl(paste0("_", node_to_delete, "$"),
               names(rv$drug_effects))
      ]
      
      # Reindexing nodes
      rv$clones_df <- reindex_nodes(rv$clones_df)
      
      showNotification(
        paste(
          "Mutation completely deleted."
        ),
        type = "message",
        duration = 6
      )
      
    } else {
      
      later_timepoints <- patient_timepoints[
        patient_timepoints >= current_tp
      ]
      
      rv$clones_df <- rv$clones_df %>%
        mutate(
          size_percent = ifelse(
            node_id == node_to_delete &
              sample_id %in% later_timepoints,
            1e-6,
            size_percent
          )
        )
      
      showNotification(
        paste(
          "Mutation biologically suppressed.",
          "It has been reduced to near-zero from this timepoint onward."
        ),
        type = "warning",
        duration = 8
      )
    }
    
    rv$clones_df <- normalize_timepoint_percentages(rv$clones_df)
    rv$objects <- build_all_objects(rv$clones_df)
    
    removeModal()
  })
  
  # Dynamic UI for drug effect inputs, which updates based on the selected drug event and impact mode (global vs per-mutation).
  output$drug_numeric_inputs <- renderUI({
    
    req(rv$clones_df, input$patient, input$selected_drug_col, input$drug_mode)
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    node_ids <- rownames(mat)
    
    key <- paste(
      input$patient,
      input$selected_drug,
      input$selected_drug_col,
      sep = "_"
    )
    
    effect <- rv$drug_effects[[key]]
    
    if (is.null(effect)) {
      
      effect <- list(
        mode = "global",
        global_value = 1,
        per_mutation = setNames(rep(1, length(node_ids)), node_ids)
      )
      
      rv$drug_effects[[key]] <- effect
    }
    
    if (input$drug_mode == "global") { # If global mode is selected, display a single numeric input to control the overall drug impact on all mutations simultaneously.
      numericInput(
        "drug_global",
        "Global drug impact",
        value = effect$global_value,
        min = 0,
        max = 2,
        step = 0.05
      )
      
    } else {
      
      # If per-mutation mode is selected, display a numeric input for each mutation to allow specific control of drug impact on each clone. 
      node_info <- get_node_labels(input$patient, rv$clones_df)
      
      tagList(
        lapply(node_ids, function(id) {
          
          current_label <- node_info$label[
            match(id, node_info$node_id)
          ]
          
          numericInput( # drug impact per mutation with values between 0 and 2 (0 = complete suppression, 1 = no effect, >1 = expansion), and a step of 0.05 for finer control
            inputId = paste0("drug_", id),
            label = current_label,
            value = effect$per_mutation[[id]],
            min = 0,
            max = 2,
            step = 0.05
          )
        })
      )
    }
  })
  
  # Observer for editing drug effect
  observeEvent(input$edit_drug_effect, {
    
    req(rv$objects, rv$clones_df, input$patient, input$selected_drug)
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    
    # If a sample is a monotimepoint -> drug event cannot be created
    if (is.null(mat) || ncol(mat) < 3) {
      
      showModal(
        modalDialog(
          title = "Drug editing not available",
          icon = icon("exclamation-triangle", class = "text-warning"),
          
          p("This patient has only one real timepoint."),
          br(),
          p("Drug impact can only be edited when at least two real timepoints exist."),
          
          tags$ul(
            tags$li("Drug events are inserted between real timepoints."),
            tags$li("With only one timepoint, no drug event can be created.")
          ),
          
          br(),
          strong("To enable drug editing:"),
          tags$ol(
            tags$li("Click 'Create new timepoint'."),
            tags$li("Then return to 'Edit drug impact'.")
          ),
          
          footer = modalButton("OK"),
          easyClose = TRUE
        )
      )
      
      return()
    }
    
    drug_cols <- seq(2, ncol(mat), by = 2)
    
    if (length(drug_cols) == 0) {
      showNotification("No drug events detected.", type = "warning")
      return()
    }
    
    drug_labels <- paste("Drug event", seq_along(drug_cols))
    
    selected_col <- drug_cols[1]
    
    node_ids <- rownames(mat)
    node_info <- get_node_labels(input$patient, rv$clones_df)
    
    key <- paste(
      input$patient,
      input$selected_drug,
      selected_col,
      sep = "_"
    )
    
    if (is.null(rv$drug_effects[[key]])) {
      
      rv$drug_effects[[key]] <- list(
        mode = "global",
        global_value = 1,
        per_mutation = setNames(rep(1, length(node_ids)), node_ids)
      )
      
    }
    
    showModal(
      modalDialog(
        title = paste("Drug impact –", input$selected_drug),
        size = "l",
        easyClose = FALSE,
        fade = TRUE,
        
        selectInput(
          "selected_drug_col",
          "Select drug event",
          choices = setNames(drug_cols, drug_labels),
          selected = selected_col
        ),
        
        helpText(
          "Drug impact rules:",
          "• 1 = no effect",
          "• <1 = clone shrinkage",
          "• >1 = clone expansion",
          "• Values are automatically rebalanced",
          "to respect clonal hierarchy constraints."
        ),
        
        h4("Mutation tree"),
        div(
          style = "
          height: 500px;
          overflow: auto;
          border: 1px solid #ddd;
          margin-bottom: 20px;
        ",
          grVizOutput("mini_tree", height = "480px")
        ),
        
        hr(),
        
        fluidRow(
          
          column(
            7,
            h4("Drug effect per clone"),
            div(
              style = "
              max-height: 400px;
              overflow-y: auto;
              padding-right: 10px;
            ",
              radioButtons(
                "drug_mode",
                "Impact mode",
                choices = c(
                  "Same impact for all mutations" = "global",
                  "Specific impact per mutation" = "per_mutation"
                ),
                selected = rv$drug_effects[[key]]$mode,
                inline = FALSE
              ),
              uiOutput("drug_numeric_inputs")
            )
          ),
          
          column(
            5,
            h4("Mutation summary"),
            DTOutput("mini_table")
          )
        ),
        
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "confirm_drug_edit",
            "Apply changes",
            class = "btn-success"
          )
        )
      )
    )
    
  })
  
  # Observer for confirming drug effect edits, which updates the reactive values storing drug effects based on user input and the selected impact mode.
  observeEvent(input$confirm_drug_edit, {
    
    req(input$selected_drug_col)
    
    key <- paste(
      input$patient,
      input$selected_drug,
      input$selected_drug_col,
      sep = "_"
    )
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    node_ids <- rownames(mat)
    
    if (is.null(rv$drug_effects[[key]])) {
      rv$drug_effects[[key]] <- list(
        mode = input$drug_mode,
        global_value = 1,
        per_mutation = setNames(rep(1, length(node_ids)), node_ids)
      )
    }
    
    rv$drug_effects[[key]]$mode <- input$drug_mode
    
    if (input$drug_mode == "global") {
      
      rv$drug_effects[[key]]$global_value <- input$drug_global
      
      rv$drug_effects[[key]]$per_mutation <-
        setNames(rep(input$drug_global, length(node_ids)), node_ids)
      
    } else {
      
      for (id in node_ids) {
        rv$drug_effects[[key]]$per_mutation[[id]] <-
          input[[paste0("drug_", id)]]
      }
    }
    
    removeModal()
  })
  
  # Minitable in the drug effect editing modal to summarize the mutations and their current percentages at the selected tiemepoint.
  output$mini_table <- renderDT({
    
    df_summary <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient) %>%
      select(node_id, mutation, parent_id) %>%
      distinct()
    
    datatable(
      df_summary,
      rownames = FALSE,
      options = list(
        pageLength = nrow(df_summary),
        scrollY = "400px",
        scrollCollapse = TRUE,
        paging = FALSE
      )
    )
  })
  
  # Zoom Observers
  observeEvent(input$zoom_in, {
    rv$zoom_level <- min(rv$zoom_level * 1.2, 5)
  })
  
  observeEvent(input$zoom_out, {
    rv$zoom_level <- max(rv$zoom_level / 1.2, 0.3)
  })
  
  observeEvent(input$zoom_reset, {
    rv$zoom_level <- 1
  })
  
  # Metadata table with dynamic editing capabilities, allowing users to modify metadata fields directly in the table, as well as add or delete metadata fields through modal dialogs.
  output$metadata_table <- renderDT({
    
    req(rv$metadata_df, input$patient)
    
    df_patient <- rv$metadata_df %>%
      filter(get_patient_id(sample_id) == input$patient)
    
    datatable(
      df_patient,
      editable = list(target = "cell", disable = list(columns =0)),
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
    
  })
  
  # Observer for handling cell edits in the metadata table, which updates the reactive metadata dataframe based on user modifications.
  observeEvent(input$metadata_table_cell_edit, {
    
    info <- input$metadata_table_cell_edit
    req(info)
    
    df_patient <- rv$metadata_df %>%
      filter(get_patient_id(sample_id) == input$patient)
    
    real_col <- colnames(df_patient)[info$col+1] # +1 because DT indexes start at 0
    
    real_sample_id <- df_patient$sample_id[info$row]
    
    rv$metadata_df[
      rv$metadata_df$sample_id == real_sample_id,
      real_col
    ] <- info$value
  })
  
  # Add metadata field observer
  observeEvent(input$add_metadata_field, {
    
    showModal(
      modalDialog(
        title = "Add metadata field",
        textInput("new_meta_name", "Field name"),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_add_meta", "Add", class = "btn-success")
        )
      )
    )
  })
  
  # Confirm add metadata field observer
  observeEvent(input$confirm_add_meta, {
    
    req(input$new_meta_name)
    
    new_col <- make.names(input$new_meta_name)
    
    if (new_col %in% colnames(rv$metadata_df)) {
      showNotification("Field already exists",
                       type = "error")
      return()
    }
    
    rv$metadata_df[[new_col]] <- NA
    
    removeModal()
    
    showNotification(
      paste0(
        "Metadata field '", input$new_meta_name,
        "' added.\nClick inside the empty cells of this new column to enter values."
      ),
      type = "message",
      duration = 6
    )
  })
  
  # Delete metadata field observer
  observeEvent(input$delete_metadata_field, {
    
    showModal(
      modalDialog(
        title = "Delete metadata field",
        selectInput(
          "meta_field_to_delete",
          "Select field",
          choices = setdiff(colnames(rv$metadata_df), "sample_id")
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_delete_meta",
                       "Delete",
                       class = "btn-danger")
        )
      )
    )
  })
  
  # Confirm delete metadata field observer
  observeEvent(input$confirm_delete_meta, {
    
    req(input$meta_field_to_delete)
    
    rv$metadata_df[[input$meta_field_to_delete]] <- NULL
    
    removeModal()
    showNotification("Metadata field deleted",
                     type = "warning")
  })
  
  # Add drug observer to allow users to add new drugs to the dataset, which can then be selected for editing drug effects in the fishplot matrix.
  observeEvent(input$add_drug, {
    
    showModal(
      modalDialog(
        title = "Add new drug",
        textInput("new_drug_name", "Drug name"),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_add_drug", "Add")
        )
      )
    )
  })
  
  # Confirm add drug observer
  observeEvent(input$confirm_add_drug, {
    
    req(input$new_drug_name)
    
    rv$available_drugs <- unique(
      c(rv$available_drugs, input$new_drug_name)
    )
    
    removeModal()
  })
  
  # Drug mode observer to update the drug effect reactive values when the user switches between global and per-mutation impact modes in the drug effect editing modal. 
  # This ensures that the correct input values are displayed and stored based on the selected mode.
  observeEvent(input$drug_mode, {
    
    req(input$selected_drug_col)
    
    key <- paste(
      input$patient,
      input$selected_drug,
      input$selected_drug_col,
      sep = "_"
    )
    
    if (!is.null(rv$drug_effects[[key]])) {
      rv$drug_effects[[key]]$mode <- input$drug_mode
    }
    
  })
  observeEvent(
    {
      input$selected_drug
      input$patient
    },
    {
      
      req(rv$objects, input$patient, input$selected_drug)
      
      mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
      if (is.null(mat) || ncol(mat) < 3) return()
      
      drug_cols <- seq(2, ncol(mat), by = 2)
      node_ids <- rownames(mat)
      
      for (col in drug_cols) {
        
        key <- paste(
          input$patient,
          input$selected_drug,
          col,
          sep = "_"
        )
        
        if (is.null(rv$drug_effects[[key]])) {
          
          rand_val <- generate_random_global_effect()
          
          rv$drug_effects[[key]] <- list(
            mode = "global",
            global_value = rand_val,
            per_mutation = setNames(
              rep(rand_val, length(node_ids)),
              node_ids
            )
          )
        }
      }
    }
  )
  
  # Cummulative matrix observer
  observeEvent(input$show_matrix, {
    rv$view <- "matrix"
  })
  
  output$matrix_table <- renderDT({
    
    mat <- fishplot_matrix_reactive()
    req(mat)
    
    df <- as.data.frame(mat)
    df <- tibble::rownames_to_column(df, "Clone")
    df[, -1] <- round(df[, -1], 4)
    
    datatable(
      df,
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        pageLength = nrow(df)
      )
    )
  })
  
  # Show intro modal automatically at app launch
  observeEvent(TRUE, {
    
    showModal(
      modalDialog(
        includeHTML("intro_text.html"),
        size = "l",
        easyClose = TRUE,
        footer = tagList(
          actionButton(
            inputId = "start_intro_tour",
            label = "Start guided tour",
            icon = icon("info-circle"),
            class = "btn-primary"
          ),
          modalButton("Close")
        )
      )
    )
    
  }, once = TRUE)
  
  #Intro tour observer 
  observeEvent(input$start_intro_tour, {
    
    removeModal()
    
    introjs(
      session,
      options = list(
        nextLabel = "Continue",
        prevLabel = "Previous",
        doneLabel = "Alright. Let's go",
        showProgress = TRUE,
        showBullets = FALSE,
        scrollToElement = TRUE,
        disableInteraction = TRUE
      )
    )
  })
  observeEvent(input$Metadata, {
    rv$view <- "metadata"
  })
  
  observeEvent(input$rename_events, {
    
    req(input$patient, rv$objects)
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    req(mat)
    
    current_labels <- colnames(mat)
    
    showModal(
      modalDialog(
        title = paste("Rename events – Patient", input$patient),
        size = "l",
        easyClose = FALSE,
        
        tagList(
          lapply(seq_along(current_labels), function(i) {
            
            textInput(
              inputId = paste0("event_label_", i),
              label = paste("Event", current_labels[i]),
              value = ifelse(
                !is.null(rv$event_labels[[input$patient]]) &&
                  length(rv$event_labels[[input$patient]]) >= i,
                rv$event_labels[[input$patient]][i],
                current_labels[i]
              )
            )
          })
        ),
        
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "confirm_event_rename",
            "Apply",
            class = "btn-success"
          )
        )
      )
    )
  })
  
  observeEvent(input$confirm_event_rename, {
    
    req(input$patient, rv$objects)
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    req(mat)
    
    new_labels <- sapply(seq_len(ncol(mat)), function(i) {
      input[[paste0("event_label_", i)]]
    })
    
    if (any(new_labels == "") || any(duplicated(new_labels))) {
      showNotification(
        "Event names must be unique and non-empty.",
        type = "error"
      )
      return()
    }
    
    rv$event_labels[[input$patient]] <- new_labels
    
    removeModal()
    
    showNotification(
      "Event labels updated.",
      type = "message"
    )
  })
}