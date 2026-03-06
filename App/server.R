#server

server <- function(input, output, session) {
  
  #Reactive Matrix
  
  fishplot_matrix_reactive <- reactive({
    
    req(input$patient, rv$objects)
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    if (is.null(mat)) return(NULL)
    
    parents <- rv$objects$all_patient_hierarchies[[input$patient]]$parent_index
    mat_final <- mat
    
    # Fix monotime
    if (ncol(mat_final) == 1) {
      
      original_name <- colnames(mat_final)[1]
      
      mat_final <- cbind(mat_final, mat_final)
      
      colnames(mat_final) <- c(
        original_name,
        paste0(original_name, "_dup")
      )
    }
    
    # Apply drug effects
    if (ncol(mat_final) >= 3) {
      
      for (col in seq(2, ncol(mat_final) - 1, by = 2)) {
        
        key <- paste(input$patient,
                     input$selected_drug,
                     col,
                     sep = "_")
        
        effect <- rv$drug_effects[[key]]
        
        if (!is.null(effect)) {
          
          if (effect$mode == "global") {
            
            mat_final[, col] <- mat_final[, col - 1] * effect$global_value
            
          } else {
            
            for (i in seq_len(nrow(mat_final))) {
              
              node_id <- rownames(mat_final)[i]
              mult <- effect$per_mutation[[node_id]]
              if (is.null(mult)) mult <- 1
              
              mat_final[i, col] <- mat_final[i, col - 1] * mult
            }
          }
        }
      }
    }
    
    mat_final <- fix_zero_reappearance(mat_final)
    mat_final <- enforce_fishplot_constraints(mat_final, parents)
    
    mat_final[!is.finite(mat_final)] <- 0
    mat_final[mat_final < 0] <- 0
    
    return(mat_final)
  })
  
  get_depth_first_order <- function(node_ids, parent_index) {
    
    children_map <- split(seq_along(parent_index), parent_index)
    
    traverse <- function(idx) {
      
      current <- node_ids[idx]
      child_idx <- children_map[[as.character(idx)]]
      
      if (is.null(child_idx)) {
        return(current)
      }
      
      # depth-first order = biological order
      c(
        current,
        unlist(lapply(child_idx, traverse))
      )
    }
    
    # roots = nodes without parent
    roots <- which(parent_index == 0)
    
    unlist(lapply(roots, traverse))
  }
  rv <- reactiveValues(
    clones_df = NULL,
    metadata_df = NULL,
    objects = NULL,
    view = "fish",
    clicked_node = NULL,
    drug_effects = list(),
    zoom_level = 1
    
  )
  
  rv$available_drugs <- c(
    "Midostaurin",
    "Venetoclax",
    "Cytarabine",
    "Azacitidine"
  )
  
  observe({
    updateSelectInput(
      session,
      "selected_drug",
      choices = rv$available_drugs,
      selected = rv$available_drugs[1]
    )
  })
  
  
  # load file server
  observeEvent(input$load_file, {
    
    req(input$load_file)
    
    my_json <- jsonlite::fromJSON(
      input$load_file$datapath,
      simplifyVector = FALSE
    )
    
    base_clones_df <- extract_all_samples(my_json) %>%
      mutate(
        node_id = as.character(node_id),
        parent_id = ifelse(is.na(parent_id), "root", as.character(parent_id)),
        size_percent = as.numeric(size_percent) * 100,
        mutation = ifelse(is.na(mutation), "none", as.character(mutation)),
        sample_id = as.character(sample_id)
      )
    
    base_metadata_df <- extract_metadata(my_json)
    
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
    
    req(rv$clones_df)
    
    existing <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient) %>%
      pull(sample_id) %>%
      unique()
    
    nums <- as.numeric(str_extract(existing, "\\d+$"))
    nums <- nums[!is.na(nums)]
    new_idx <- ifelse(length(nums) == 0, 1, max(nums) + 1)
    new_sample <- paste0(input$patient, "-", sprintf("%03d", new_idx))
    
    base_rows <- rv$clones_df %>%
      filter(sample_id == existing[1]) %>%
      mutate(
        sample_id = new_sample,
        size_percent = 0
      )
    
    rv$clones_df <- bind_rows(rv$clones_df, base_rows)
    rv$objects <- build_all_objects(rv$clones_df)
    
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
        title = paste(
          "Edit all percentages – Patient",
          input$patient,
          "| Timepoint:",
          input$selected_timepoint
        ),
        size = "l",
        easyClose = FALSE,
        fade = TRUE,
        
        div(
          style = "background-color:#f8f9fa; padding:15px; border-radius:8px; margin-bottom:20px;",
          
          h4("Current editing context"),
          
          p(
            strong("IMPORTANT :You are currently editing timepoint: "),
            strong(input$selected_timepoint)
          ),
          
          br(),
          
          p(
            strong("To change timepoint:"),
            br(),
            "Use the timepoint selector in the left sidebar ",
            "and select another timepoint before clicking 'Edit all percentages'."
          ),
          
          hr(),
          
          h4("How to edit clone percentages"),
          
          tags$ol(
            tags$li(
              strong("Modify clone percentages in the numeric fields."),
              " Values must be between 0 and 100."
            ),
            tags$li(
              strong("Parent clones must be ≥ the sum of their children."),
              " If not, automatic rebalancing will occur."
            ),
            tags$li(
              strong("Click 'Apply changes' to validate.")
            )
          ),
          
          hr(),
          
          p(
            strong("Important biological rules applied automatically:"),
            br(),
            "• Total clone size cannot exceed 100%.",
            br(),
            "• Parent-child hierarchy constraints are enforced.",
            br(),
            "• Zero reappearance correction is applied.",
            br(),
            "• Root clone is fixed at 100%."
          )
        ),
        
        # top section with tree
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
          
          # left column
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
          
          # Right column
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
  
  
  observeEvent(input$confirm_bulk_edit, {
    
    df_tp <- rv$clones_df %>%
      filter(sample_id == input$selected_timepoint)
    
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
        
        div(
          style = "background-color:#f8f9fa; padding:15px; border-radius:8px; margin-bottom:20px;",
          
          h4("Current editing context"),
          
          p(
            strong("IMPORTANT : You are currently editing timepoint: "),
            strong(input$selected_timepoint)
          ),
          
          br(),
          
          p(
            strong("To change timepoint:"),
            br(),
            "Use the timepoint selector in the left sidebar ",
            "and choose another timepoint before clicking 'Create mutation'."
          ),
          
          hr(),
          
          h4("How to create a new mutation"),
          
          tags$ol(
            tags$li(
              strong("Select a parent node: "),
              "Choose the clone from which the new mutation evolves."
            ),
            tags$li(
              strong("Enter the mutation name: "),
              "Example: TP53 R175H."
            ),
            tags$li(
              strong("Set the percentage for the current timepoint: "),
              "This represents the clone size at ",
              strong(input$selected_timepoint), "."
            )
          ),
          
          hr(),
          
          p(
            strong("Important biological rule:"),
            br(),
            "• The mutation will appear at the current timepoint.",
            br(),
            "• It will be automatically set to 0% at all other timepoints.",
            br(),
            "• Parent-child hierarchy constraints will be enforced automatically."
          )
        ),
        
        
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
  
  observeEvent(input$confirm_add_mutation, {
    
    req(rv$clones_df)
    removeModal()
    
    # new node ID generation
    # IDs existants uniquement pour le patient courant
    existing_ids <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient) %>%
      pull(node_id) %>%
      unique() %>%
      as.numeric()
    
    # Sécurité si vide
    if (length(existing_ids) == 0) {
      new_id <- 0
    } else {
      new_id <- max(existing_ids) + 1
    }
    
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
    
    if (input$new_parent == new_id) {
      showNotification("Parent cannot be itself", type = "error")
      return()
    }
    
    patient_samples <- rv$clones_df %>%
      filter(get_patient_id(sample_id) == input$patient) %>%
      pull(sample_id) %>%
      unique()
    
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
    rv$objects <- build_all_objects(rv$clones_df)
  })
  
  
  
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
  
  output$main_view <- renderUI({
    
    if (rv$view == "fish") {
      introBox(
        plotOutput("fishplot", height = "700px"),
        data.step = 11,
        data.intro = get_step(11)
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
          introBox(
            grVizOutput("mutation_tree", height = "800px"),
            data.step = 12,
            data.intro = get_step(12)
          )
        )
      )
      
      
    } else if (rv$view == "format") {
      box(width = 12, includeMarkdown(file.path("..","Texts_and_Extras","data_format.md")))
      
    }
    
    else if (rv$view == "metadata") {
      tagList(
        h3(paste("Metadata – Patient", input$patient)),
        introBox(
          DTOutput("metadata_table"),
          data.step = 15,
          data.intro = get_step(15)
        ),
        
        br(),
        
        div(
          style = "display:flex; gap:10px;",
          actionButton("add_metadata_field", "Add metadata field", class = "btn-primary"),
          actionButton("delete_metadata_field", "Delete metadata field", class = "btn-danger")
        )
      )
    }
    
    else if (rv$view =="matrix") {
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
        
        introBox(
          DTOutput("matrix_table"),
          data.step = 13,
          data.intro = get_step(13)
        )
      )
    }
    else {
      DTOutput("table")
    }
  })
  
  output$mini_tree <- renderGrViz({
    req(rv$clones_df, rv$objects)
    plot_mutation_tree_shiny(
      input$patient,
      rv$clones_df,
      rv$objects$clone_palette
    )
  })
  
  output$table <- renderDT({
    
    req(rv$clones_df, rv$objects, input$selected_timepoint)
    
    patient_id <- input$patient
    
    mat <- rv$objects$concat_by_patient_hierarchy[[patient_id]]
    parents <- rv$objects$all_patient_hierarchies[[patient_id]]$parent_index
    node_ids <- rownames(mat)
    
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
  
  
  output$fishplot <- renderPlot({
    
    req(input$patient, rv$objects)
    
    plot_fishplot(
      patient = input$patient,
      concat_by_patient_hierarchy = rv$objects$concat_by_patient_hierarchy,
      all_patient_hierarchies = rv$objects$all_patient_hierarchies,
      mutation_lookup = rv$objects$mutation_lookup,
      clone_palette = global_clone_palette,
      input = input,
      rv = rv,
      mini = FALSE
    )
    
  })
  
  
  output$mutation_tree <- renderGrViz({
    
    req(rv$clones_df, rv$objects)
    
    zoom <- rv$zoom_level   # ← IMPORTANT : crée la dépendance
    
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

        el.querySelectorAll('.node').forEach(function(node) {
          node.onclick = function() {
            Shiny.setInputValue(
              'mutation_clicked',
              node.querySelector('title').textContent,
              {priority: 'event'}
            );
          };
        });

      }
    ", zoom))
  })
  
  
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
            drug_effects = patient_drug_effects   # 🔹 AJOUT ICI
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
  
  observeEvent(input$mutation_clicked, {
    showModal(
      modalDialog(
        title = paste("Set percentage for", input$mutation_clicked),
        numericInput(
          "edit_percent",
          "Percentage",
          value = 10,
          min = 0,
          max = 100,
          step = 0.1
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_edit", "Apply", class = "btn-success")
        )
      )
    )
  })
  
  observeEvent(input$confirm_edit, {
    
    req(rv$clones_df)
    removeModal()
    
    rv$clones_df <- rv$clones_df %>%
      mutate(
        size_percent = ifelse(
          sample_id == input$selected_timepoint &
            node_id == input$mutation_clicked,
          input$edit_percent,
          size_percent
        )
      )
    
    rv$objects <- build_all_objects(rv$clones_df)
  })
  
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
    
    active_nodes <- df_patient %>%
      group_by(node_id) %>%
      summarise(total = sum(size_percent, na.rm = TRUE)) %>%
      filter(total > 0) %>%
      pull(node_id)
    
    deletable_clones <- df_structure %>%
      filter(
        parent_id != "root",
        !node_id %in% nodes_with_children,
        node_id %in% active_nodes
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
        
        p(
          strong("IMPORTANT : You are currently editing timepoint: "),
          strong(input$selected_timepoint)
        ),
        
        br(),
        
        p(
          strong("To change timepoint:"),
          br(),
          "Use the timepoint selector in the left sidebar ",
          "and choose another sample before clicking 'Create mutation'."
        ),
        
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
    
    rv$objects <- build_all_objects(rv$clones_df)
    
    removeModal()
  })
  
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
    
    if (input$drug_mode == "global") {
      
      numericInput(
        "drug_global",
        "Global drug impact",
        value = effect$global_value,
        min = 0,
        max = 2,
        step = 0.05
      )
      
    } else {
      
      # --- MODE PER MUTATION ---
      node_info <- get_node_labels(input$patient, rv$clones_df)
      
      tagList(
        lapply(node_ids, function(id) {
          
          current_label <- node_info$label[
            match(id, node_info$node_id)
          ]
          
          numericInput(
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
  
  observeEvent(input$edit_drug_effect, {
    
    req(rv$objects, rv$clones_df, input$patient, input$selected_drug)
    
    mat <- rv$objects$concat_by_patient_hierarchy[[input$patient]]
    
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
  
  
  observeEvent(input$zoom_in, {
    rv$zoom_level <- min(rv$zoom_level * 1.2, 5)
  })
  
  observeEvent(input$zoom_out, {
    rv$zoom_level <- max(rv$zoom_level / 1.2, 0.3)
  })
  
  
  observeEvent(input$zoom_reset, {
    rv$zoom_level <- 1
  })
  
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
    
    showNotification("Metadata field added",
                     type = "message")
  })
  
  
  observeEvent(input$delete_metadata_field, {
    
    # Vérification de l'existence du metadata
    if (is.null(rv$metadata_df) || ncol(rv$metadata_df) <= 1) {
      showNotification(
        "No metadata fields available for deletion.",
        type = "warning"
      )
      return()
    }
    
    selectable_fields <- setdiff(colnames(rv$metadata_df), "sample_id")
    
    if (length(selectable_fields) == 0) {
      showNotification(
        "No removable metadata fields found.",
        type = "warning"
      )
      return()
    }
    
    showModal(
      modalDialog(
        title = "Delete metadata field",
        selectInput(
          "meta_field_to_delete",
          "Select field",
          choices = selectable_fields
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "confirm_delete_meta",
            "Delete",
            class = "btn-danger"
          )
        ),
        easyClose = TRUE
      )
    )
  })
  
  observeEvent(input$confirm_delete_meta, {
    
    req(input$meta_field_to_delete)
    
    rv$metadata_df[[input$meta_field_to_delete]] <- NULL
    
    removeModal()
    showNotification("Metadata field deleted",
                     type = "warning")
  })
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
  
  observeEvent(input$confirm_add_drug, {
    
    req(input$new_drug_name)
    
    new_drug <- trimws(input$new_drug_name)
    
    if (new_drug == "") {
      showNotification("Drug name cannot be empty", type = "error")
      return()
    }
    
    if (new_drug %in% rv$available_drugs) {
      showNotification("Drug already exists", type = "warning")
      return()
    }
    
    rv$available_drugs <- c(rv$available_drugs, new_drug)
    
    updateSelectInput(
      session,
      "selected_drug",
      choices = rv$available_drugs,
      selected = new_drug
    )
    
    removeModal()
    
    showNotification(
      paste("Drug", new_drug, "added successfully"),
      type = "message"
    )
    
  })
  
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
  observeEvent(input$show_matrix, {
    rv$view <- "matrix"
  })
  
  output$matrix_table <- renderDT({
    
    mat <- fishplot_matrix_reactive()
    req(mat)
    
    df <- as.data.frame(mat)
    
    n_cols <- ncol(df)
    
    new_names <- character(n_cols)
    real_tp_index <- 1
    drug_index <- 1
    
    for (i in seq_len(n_cols)) {
      
      if (i %% 2 == 1) {
        new_names[i] <- paste0("t", real_tp_index)
        real_tp_index <- real_tp_index + 1
      } else {
        new_names[i] <- paste0("drug_event_", drug_index)
        drug_index <- drug_index + 1
      }
    }
    
    colnames(df) <- new_names
    
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
        includeHTML(file.path("..","Texts_and_Extras","intro_text.html")),
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
}