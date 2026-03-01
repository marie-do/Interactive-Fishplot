library(methods)

fishplot_source <- file.path("fishplot_manual", "fishplot-master")

if(!dir.exists(fishplot_source)) {
  cat("Folder not found. Contents of C:/fishplot_manual/:\n")
  print(list.files("C:/fishplot_manual/"))
  stop("Adjust the path")
}

r_files <- list.files(file.path(fishplot_source, "R"),
                      pattern = "\\.R$",
                      full.names = TRUE)

cat("R files found:\n")
print(r_files)

object_file <- r_files[grepl("object.R", r_files)]
if(length(object_file) > 0) {
  cat("Loading object.R...\n")
  source(object_file)
}

if(exists("initFishClass")) {
  cat("initFishClass found, initializing...\n")
  initFishClass()
} else {
  cat("initFishClass not found. Available functions:\n")
  print(ls(.GlobalEnv)[grepl("init|fish|Fish", ls(.GlobalEnv))])
}

other_files <- r_files[!grepl("object.R", r_files)]
for(file in other_files) {
  cat("Loading:", basename(file), "\n")
  source(file)
}

required_packages <- c("ape", "gtools")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

cat("Setup complete! You can now use fishplot.\n")