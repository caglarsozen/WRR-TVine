# Short execution check using the same engine as the full experiment.
# Run from the repository root in RStudio or with Rscript.
local({
  args <- commandArgs(trailingOnly = TRUE)
  prices <- if (length(args)) args[1L] else "data/prices.csv"
  source("tests/test_invariants.R", local = new.env(parent = .GlobalEnv), encoding = "UTF-8")
  engine <- new.env(parent = .GlobalEnv)
  source("R/wrr_tvine.R", local = engine, encoding = "UTF-8")
  engine$wrr_launch("smoke", prices_csv = prices)
})
