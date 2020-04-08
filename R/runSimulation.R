#' @export
runSimulation <- function() {
  appDir <- system.file("shiny-examples", "project_2", package = "project_2")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `project_2`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
