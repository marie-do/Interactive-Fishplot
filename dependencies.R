#  List of required packages for the app + fishplot code

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
  "DiagrammeR"
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
