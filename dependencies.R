# List of required packages for the app + fishplot code
# If fishplot can be installed from GitHub, try the step 1 first:

# Step 1: Try to install fishplot from GitHub 

fishplot_available <- requireNamespace("fishplot", quietly = TRUE)

if (!fishplot_available) {
  
  message("fishplot not found. Attempting GitHub installation...")
  
  # Ensure devtools is installed
  if (!requireNamespace("devtools", quietly = TRUE)) {
    tryCatch({
      install.packages("devtools")
    }, error = function(e) {
      message("devtools installation failed:")
      message(e$message)
    })
  }
  
  # Try installing fishplot
  tryCatch({
    
    devtools::install_github(
      "chrisamiller/fishplot",
      upgrade = "never",
      dependencies = TRUE,
      quiet = TRUE
    )
    
    message("fishplot GitHub installation completed.")
    
  }, error = function(e) {
    
    message("GitHub installation of fishplot failed:")
    message(e$message)
    
  })
}

# Re-check AFTER installation attempt
fishplot_available <- requireNamespace("fishplot", quietly = TRUE)

# Load ONLY if truly available
if (fishplot_available) {
  library(fishplot)
  message("fishplot loaded successfully.")
} else {
  message("fishplot not available. Continuing without it.")
}

required_packages <- c(
  "jsonlite",
  "dplyr",
  "purrr",
  "tidyr",
  "tibble",
  "shinyjs",
  "rintrojs",
  "readr",
  "shiny",
  "shinydashboard",
  "stringr",
  "DT",
  "DiagrammeR",
  "markdown"
)

# install missing packages

new.packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(new.packages)) {
  install.packages(new.packages)
}

rm(new.packages)

library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(shinyjs)
library(rintrojs)
library(readr)
library(shiny)
library(shinydashboard)
library(stringr)
library(DT)
library(DiagrammeR)
library(dplyr)
library(markdown)
