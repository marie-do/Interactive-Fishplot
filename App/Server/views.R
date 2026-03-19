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
    box(width = 12, includeMarkdown("../Texts_and_Extras/data_format.md"))
    
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

# Minitable in the drug effect editing modal to summarize the mutations and their current percentages at the selected tiemepoint.
output$mini_table <- renderDT({
  
  df_summary <- rv$clones_df %>%
    filter(get_patient_id(sample_id) == input$patient) %>%
    select(node_id, mutation, parent_id) %>%
    distinct() %>%
    filter(mutation != "root")
  
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
      includeHTML("../Texts_and_Extras/intro_text.html"),
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