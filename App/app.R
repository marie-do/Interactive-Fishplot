app_source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

app_dir <- if (!is.null(app_source_file)) {
	dirname(normalizePath(app_source_file))
} else if (file.exists("ui.R")) {
	getwd()
} else if (file.exists(file.path("app", "ui.R"))) {
	normalizePath("app")
} else {
	stop("Unable to determine the app directory.")
}

setwd(app_dir)

source("../global.R")
source("ui.R")
source("server/server.R")

app <- shinyApp(ui = ui, server = server)

# Start the app when this file is sourced in an interactive console.
# When called via shiny::runApp("App"), return the app object without auto-running.
called_from_runapp <- any(vapply(sys.calls(), function(cl) {
	fn <- tryCatch(as.character(cl[[1]]), error = function(e) "")
	any(grepl("runApp", fn, fixed = TRUE))
}, logical(1)))

if (interactive() && !called_from_runapp) {
	print(app)
}

app
