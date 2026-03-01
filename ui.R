# ui

help_steps <- read_delim(
  "help_tour.csv",
  delim = ";",
  col_types = cols(.default = "c"),
  trim_ws = TRUE
)

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = tags$span(
      "Fishplot Interactive",
      'data-step' = 1,
      'data-intro' = get_step(1),
      style = "display:inline-block;"
    ),
    
    dropdownMenu(
      type = "notifications", 
      headerText = strong("HELP"), 
      icon = icon("question"), 
      badgeStatus = NULL,
      notificationItem(
        text = (help_steps$text[1]),
        icon = icon("upload")
      ),
      notificationItem(
        text = help_steps$text[2],
        icon = icon("user")
      ),
      notificationItem(
        text = help_steps$text[3],
        icon = icon("edit")
      ),
      notificationItem(
        text = help_steps$text[4],
        icon = icon("plus")
      ),
      notificationItem(
        text = help_steps$text[5],
        icon = icon("edit")
      ),
      notificationItem(
        text = help_steps$text[6],
        icon = icon("trash")
      ),
      notificationItem(
        text = help_steps$text[7],
        icon = icon("chart-area")
      ),
      notificationItem(
        text = help_steps$text[8],
        icon = icon("save")
      )
    ),
    tags$li(
      a(
        strong("ABOUT THE APP"),
        height = 40,
        href = "https://github.com/marie-do/Interactive-Fishplot",
        title = "",
        target = "_blank"
      ),
      class = "dropdown"
    )
  ),
  
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
        class = "btn-success",
        style = "width: 100%;"
      ),
      data.step = 16,
      data.intro = get_step(16)
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
        actionButton("new_timepoint","Create new timepoint",
                     icon=icon("plus"), class="btn-warning"),
        actionButton("add_mutation","Create new mutation",
                     icon=icon("plus"), class="btn-warning")
      ),
      data.step = 7,
      data.intro = get_step(7)
    ),
    
    introBox(
      actionButton(
        "edit_drug_effect",
        "Edit drug impact",
        icon = icon("flask"),
        class = "btn-info"
      ),
      data.step = 8,
      data.intro = get_step(8)
    ),
    
    introBox(
      tagList(
        actionButton("delete_timepoint","Delete timepoint",
                     icon=icon("trash"), class="btn-danger"),
        actionButton("delete_mutation","Delete mutation",
                     icon=icon("trash"), class="btn-danger")
      ),
      data.step = 9,
      data.intro = get_step(9)
    ),
    
    introBox(
      actionButton(
        "add_drug",
        "Add new drug",
        icon = icon("plus"),
        class = "btn-primary"
      ),
      data.step = 10,
      data.intro = get_step(10)
    ),
    
    br(), br(),
    
    sidebarMenu(
      id = "sidebar_menu",
      
      menuItem(
        "Visualization",
        tabName = "viz",
        icon = icon("chart-area")
      ),
      
      menuItem(
        "Data format",
        tabName = "format",
        icon = icon("info-circle")
      ),
      
      introBox(
        actionButton(
          "Metadata",
          label = "Metadata",
          icon = icon("database"),
          class = "btn btn-secondary"
        ),
        `data-step` = 15,
        `data-intro` = get_step(15)
      )
    )
    
  ),
  
  dashboardBody(
    
    useShinyjs(),
    introjsUI(),
    
    tags$head(
      tags$style(HTML("
      
      .navbar-nav > .notifications-menu > .dropdown-menu > li > .menu > li > a {
        white-space: normal !important;
        overflow: visible !important;
        text-overflow: unset !important;
        height: auto !important;
      }

      .navbar-nav > .notifications-menu > .dropdown-menu > li > .menu > li {
        height: auto !important;
      }

      .navbar-nav > .notifications-menu > .dropdown-menu {
        max-height: 600px !important;
        overflow-y: auto !important;
      }
      table.dataTable tbody td {
        cursor: pointer;
      }
      table.dataTable tbody td:not(:first-child):hover {
        background-color: #e8f4ff !important;
      }
    "))
    ),
    
    fluidRow(
      column(
        width = 12,
        
        actionButton(
          inputId = "show_fish",
          label   = "FISHPLOT",
          class   = "btn-success",
          `data-step`  = 11,
          `data-intro` = get_step(11)
        ),
        
        actionButton(
          inputId = "show_tree",
          label   = "MUTATION TREE",
          class   = "btn-success",
          `data-step`  = 12,
          `data-intro` = get_step(12)
        ),
        
        actionButton(
          inputId = "show_data",
          label   = "DATA INFO",
          class   = "btn-success",
          `data-step`  = 13,
          `data-intro` = get_step(13)
        ),
        
        actionButton(
          inputId = "show_matrix",
          label   = "DATA MATRIX",
          class   = "btn-success",
          `data-step`  = 14,
          `data-intro` = get_step(14)
        )
      )
    ),
    
    br(),
    
    uiOutput("main_view")
  )
)