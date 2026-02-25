#ui

intro_steps <- read_delim(
  "fishplot_tour.csv",
  delim = ";",
  col_types = cols()
)

get_step <- function(n) {
  step_text <- intro_steps$text[intro_steps$step == n]
  if (length(step_text) == 0) return("")
  step_text[[1]]
}


ui <- dashboardPage(
  skin = "blue",
  
  
  dashboardHeader(title = "Fishplot Interactive Editor"),
  
  dashboardSidebar(
    width = 300,
    
    introBox(
      fileInput(
        "load_file",
        "Load JSON data file",
        accept = ".json"
      ),
      data.step = 2,
      data.intro = get_step(2)
    ),
    
    introBox(
      downloadButton(
        "save_file",
        "Save modified data",
        class = "btn-success"
      ),
      data.step = 13,
      data.intro = get_step(13)
    ),
    
    hr(),
    
    introBox(
      selectInput(
        "patient",
        "Select patient",
        choices = NULL
      ),
      data.step = 3,
      data.intro = get_step(3)
    ),
    
    introBox(
      selectInput(
        "selected_timepoint",
        "Select timepoint",
        choices = NULL
      ),
      data.step = 4,
      data.intro = get_step(4)
    ),
    
    introBox(
      selectInput(
        "selected_drug",
        "Select drug",
        choices = NULL
      ),
      data.step = 5,
      data.intro = get_step(5)
    ),
    
    introBox(
      actionButton(
        "edit_all",
        "Edit all percentages",
        icon = icon("edit"),
        class = "btn-info"
      ),
      data.step = 6,
      data.intro = get_step(6)
    ),
    
    introBox(
      tagList(
        actionButton("new_timepoint","Create new timepoint",icon=icon("plus"),class="btn-warning"),
        actionButton("add_mutation","Create new mutation",icon=icon("plus"),class="btn-warning"),
        actionButton("edit_drug_effect","Edit drug impact",icon=icon("flask"),class="btn-info")
      ),
      data.step = 7,
      data.intro = get_step(7)
    ),
    
    introBox(
      tagList(
        actionButton("delete_timepoint","Delete timepoint",icon=icon("trash"),class="btn-danger"),
        actionButton("delete_mutation","Delete mutation",icon=icon("trash"),class="btn-danger")
      ),
      data.step = 8,
      data.intro = get_step(8)
    ),
    
    introBox(
      actionButton(
        "add_drug",
        "Add new drug",
        icon = icon("plus"),
        class = "btn-primary"
      ),
      data.step = 9,
      data.intro = get_step(9)
    ),
    
    br(), br(),
    
    sidebarMenu(
      id = "sidebar_menu",
      menuItem("Visualization", tabName = "viz", icon = icon("chart-area")),
      menuItem("Data format", tabName = "format", icon = icon("info-circle")),
      menuItem("Metadata", tabName = "metadata", icon = icon("database"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    introjsUI(),
    introBox(
      tags$div(id = "tour_start"),
      data.step = 1,
      data.intro = get_step(1)
    ),
    tags$head(
      tags$style(HTML("
      .modal-lg {
        width: 90% !important;
        max-width: 1400px !important;
      }
    "))
    ),
    fluidRow(
      introBox(
        column(
          12,
          actionButton("show_fish", "FISHPLOT", class = "btn-success"),
          actionButton("show_tree", "MUTATION TREE", class = "btn-success"),
          actionButton("show_data", "DATA INFO", class = "btn-success"),
          actionButton("show_matrix", "DATA MATRIX", class = "btn-success")
        ),
        data.step = 10,
        data.intro = get_step(10)
      )
    ),
    br(),
    uiOutput("main_view")
  )
)