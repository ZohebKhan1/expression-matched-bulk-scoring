#' Load bulk module scoring functions
#'
#' Sources the three scoring function files into the global environment.
#' Run from the repository root with `source("load_functions.R")`.
#'
#' @param path Directory containing `module_score.R`, `permutation.R`, and
#'   `plot_module_score.R`. Defaults to `"functions"`.
load_bulk_module_score_functions <- function(path = "functions") {
  source(file.path(path, "module_score.R"), local = FALSE)
  source(file.path(path, "permutation.R"), local = FALSE)
  source(file.path(path, "plot_module_score.R"), local = FALSE)
  invisible(TRUE)
}

load_bulk_module_score_functions()
