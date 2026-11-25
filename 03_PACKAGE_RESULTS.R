# Package a completed full experiment without fitting any models.
local({
  args <- commandArgs(trailingOnly = TRUE)
  root <- if (length(args)) args[1L] else file.path("results", "full")
  archive <- if (length(args) > 1L) args[2L] else "WRR_TVINE_FULL_RESULTS.zip"
  engine <- new.env(parent = .GlobalEnv)
  source("R/wrr_tvine.R", local = engine, encoding = "UTF-8")
  engine$wrr_package_results(root = root, log_file = "FULL_RUN_CONSOLE.log", archive = archive)
})
