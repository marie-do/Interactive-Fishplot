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