source("../global.R")
source("ui.R")
source("Server/server.R")

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
