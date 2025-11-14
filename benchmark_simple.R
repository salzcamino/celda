#!/usr/bin/env Rscript
# Simple Benchmarking Script using reportCeldaCGRun
# This demonstrates the integrated workflow with parallelization

library(celda)
library(SingleCellExperiment)

cat("Simple celda Benchmarking with reportCeldaCGRun\n")
cat("================================================\n\n")

set.seed(12345)

# Detect OS and set parallel cores accordingly
# mclapply doesn't work on Windows, so we skip parallel tests on Windows
isWindows <- .Platform$OS.type == "windows"
if (isWindows) {
  cat("NOTE: Windows detected. Parallel tests (nCores > 1) will be skipped.\n")
  cat("      mclapply() is not supported on Windows.\n\n")
  parallelCores <- 1
} else {
  parallelCores <- 4
}

# Generate simulated data
cat("Generating simulated dataset...\n")
sim <- simulateContamination(
  C = 10,
  G = 200,
  N = 500,
  beta = 0.1,
  seed = 12345
)

sce <- SingleCellExperiment(assays = list(counts = sim$observedCounts))

# Select features (required before running celdaGridSearch)
# Use minCount=1, minCell=1 to ensure all features have at least some counts
sce <- selectFeatures(sce, minCount = 1, minCell = 1)

cat("Dataset:", nrow(sce), "genes x", ncol(sce), "cells\n\n")

# Test 1: Serial execution
cat("Test 1: Serial execution (nCores = 1)\n")
time_serial <- system.time({
  model_serial <- celdaGridSearch(
    x = sce,
    paramsTest = list(L = c(5, 10), K = c(8, 10)),
    paramsFixed = list(sampleLabel = sim$z),
    maxIter = 50,
    nchains = 1,
    cores = 1,
    verbose = FALSE,
    perplexity = FALSE
  )

  best_serial <- selectBestModel(model_serial)

  # Test recursive splits with serial
  split_module_serial <- recursiveSplitModule(
    x = best_serial,
    maxL = 15,
    nCores = 1,
    verbose = FALSE
  )

  split_cell_serial <- recursiveSplitCell(
    x = best_serial,
    maxK = 15,
    nCores = 1,
    verbose = FALSE
  )
})
cat("Total time:", round(time_serial["elapsed"], 2), "seconds\n\n")

# Test 2: Parallel execution
if (!isWindows) {
  cat("Test 2: Parallel execution (nCores = 4)\n")
  time_parallel <- system.time({
    model_parallel <- celdaGridSearch(
      x = sce,
      paramsTest = list(L = c(5, 10), K = c(8, 10)),
      paramsFixed = list(sampleLabel = sim$z),
      maxIter = 50,
      nchains = 1,
      cores = parallelCores,
      verbose = FALSE,
      perplexity = FALSE
    )

    best_parallel <- selectBestModel(model_parallel)

    # Test recursive splits with parallelization
    split_module_parallel <- recursiveSplitModule(
      x = best_parallel,
      maxL = 15,
      nCores = parallelCores,
      verbose = FALSE
    )

    split_cell_parallel <- recursiveSplitCell(
      x = best_parallel,
      maxK = 15,
      nCores = parallelCores,
      verbose = FALSE
    )
  })
  cat("Total time:", round(time_parallel["elapsed"], 2), "seconds\n\n")
} else {
  cat("Test 2: Parallel execution SKIPPED (not supported on Windows)\n\n")
  time_parallel <- time_serial
}

# Summary
cat("================================================\n")
cat("RESULTS:\n")
cat("  Serial time:", round(time_serial["elapsed"], 2), "seconds\n")
if (!isWindows) {
  cat("  Parallel time:", round(time_parallel["elapsed"], 2), "seconds\n")
  cat("  Speedup:", round(time_serial["elapsed"] / time_parallel["elapsed"], 2), "x\n")
} else {
  cat("  Parallel tests: SKIPPED (Windows does not support mclapply)\n")
}
cat("================================================\n\n")

# Optional: Generate HTML report with parallelization
cat("To generate a full HTML report with parallelization, run:\n\n")
cat("reportCeldaCGRun(\n")
cat("  sce = sce,\n")
cat("  L = 10,\n")
cat("  K = 10,\n")
cat("  maxL = 15,\n")
cat("  maxK = 15,\n")
cat("  nCores = 4,  # Use parallelization!\n")
cat("  outputFile = 'celda_report.html',\n")
cat("  outputDir = '.'\n")
cat(")\n")
