app_directory <- normalizePath(".", winslash = "/", mustWork = TRUE)
source(file.path(app_directory, "UI.R"), local = TRUE)
source(file.path(app_directory, "server.R"), local = TRUE)
shiny::shinyApp(ui = ui, server = server)
