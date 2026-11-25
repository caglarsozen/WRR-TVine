# Complete empirical specification, with reusable checkpoints.
# Run from the repository root in RStudio or with Rscript.
local({
  args <- commandArgs(trailingOnly = TRUE)
  prices <- if (length(args)) args[1L] else "data/prices.csv"
  engine <- new.env(parent = .GlobalEnv)
  source("R/wrr_tvine.R", local = engine, encoding = "UTF-8")
  engine$wrr_launch("full", prices_csv = prices)
})
