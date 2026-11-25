# Independent base R audit of the completed experiment.
# Accepts an extracted results directory or a results ZIP.
local({
  args <- commandArgs(trailingOnly = TRUE)
  input <- if (length(args)) args[1L] else if (dir.exists(file.path("results", "full"))) {
    file.path("results", "full")
  } else "WRR_TVINE_FULL_RESULTS.zip"
  output <- if (length(args) > 1L) args[2L] else file.path("analysis_output", "verification")
  run <- function() {
    root <- input
    if (!dir.exists(root)) {
      if (!file.exists(input) || !grepl("[.]zip$", input, ignore.case = TRUE))
        stop("Supply a completed results directory or a results ZIP: ", input)
      entries <- gsub("\\\\", "/", unzip(input, list = TRUE)$Name)
      if (any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", entries)))
        stop("The archive contains unsafe paths.")
      extracted <- tempfile(pattern = "wrr_verify_", tmpdir = tempdir(check = TRUE))
      if (!dir.create(extracted)) stop("Could not create temporary extraction directory")
      on.exit(unlink(extracted, recursive = TRUE), add = TRUE)
      unzip(input, exdir = extracted)
      root <- extracted
    }
    if (!file.exists(file.path(root, "RUN_STATUS.txt"))) {
      statuses <- list.files(root, pattern = "^RUN_STATUS[.]txt$", recursive = TRUE, full.names = TRUE)
      if (length(statuses) != 1L) stop("Could not identify one completed results directory under ", root)
      root <- dirname(statuses[1L])
    }
    verification <- new.env(parent = .GlobalEnv)
    source("R/verify_results.R", local = verification, encoding = "UTF-8")
    verification$wrr_verify_results(root, file.path(output, "core"))
    verification$wrr_supplementary(root, file.path(output, "supplementary"))
  }
  run()
})
