# Production-safe dependency loader.
# Important: never install packages at app startup on shinyapps.io.

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
  "markdown",
  "methods"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running/deploying the app."
    )
  )
}

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
library(markdown)
library(methods)
