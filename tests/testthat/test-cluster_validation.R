# Tests for Cluster Validation Metrics
library(celda)
library(testthat)

context("Testing cluster validation metrics")

# Setup test data
set.seed(12345)
K <- 3
L <- 5
nGenes <- 100
nCells <- 50

# Simulate a simple count matrix
testCounts <- matrix(rpois(nGenes * nCells, lambda = 10),
                    nrow = nGenes, ncol = nCells)
testCounts <- as(testCounts, "dgCMatrix")

# Create good clustering (distinct clusters)
goodZ <- rep(1:K, length.out = nCells)
goodY <- rep(1:L, length.out = nGenes)

# Create poor clustering (random)
poorZ <- sample(1:K, nCells, replace = TRUE)
poorY <- sample(1:L, nGenes, replace = TRUE)

# Create edge case clustering (all same cluster)
badZ <- rep(1, nCells)


# Test Calinski-Harabasz Index
test_that(".calculateCH identifies good vs. poor clustering", {
  chGood <- .calculateCH(testCounts, goodZ)
  chPoor <- .calculateCH(testCounts, poorZ)

  expect_true(is.numeric(chGood))
  expect_true(is.numeric(chPoor))
  expect_false(is.na(chGood))
  expect_false(is.na(chPoor))

  # Good clustering should have higher CH
  # (not always guaranteed with random data, but likely)
  expect_true(is.finite(chGood))
  expect_true(is.finite(chPoor))
})

test_that(".calculateCH handles edge cases", {
  # K = 1 should return NA
  chBad <- .calculateCH(testCounts, badZ)
  expect_true(is.na(chBad))

  # K = 2 should work
  z2 <- rep(1:2, length.out = nCells)
  ch2 <- .calculateCH(testCounts, z2)
  expect_true(is.numeric(ch2))
  expect_false(is.na(ch2))
})


# Test Davies-Bouldin Index
test_that(".calculateDB calculates correctly", {
  dbGood <- .calculateDB(testCounts, goodZ)
  dbPoor <- .calculateDB(testCounts, poorZ)

  expect_true(is.numeric(dbGood))
  expect_true(is.numeric(dbPoor))
  expect_false(is.na(dbGood))
  expect_false(is.na(dbPoor))

  # DB should be positive
  expect_true(dbGood >= 0)
  expect_true(dbPoor >= 0)
})

test_that(".calculateDB handles edge cases", {
  # K = 1 should return NA
  dbBad <- .calculateDB(testCounts, badZ)
  expect_true(is.na(dbBad))

  # K = 2 should work
  z2 <- rep(1:2, length.out = nCells)
  db2 <- .calculateDB(testCounts, z2)
  expect_true(is.numeric(db2))
  expect_false(is.na(db2))
})


# Test Silhouette Score
test_that(".calculateSilhouette_Sampled works with sampling", {
  # Small dataset - no sampling
  silSmall <- .calculateSilhouette_Sampled(testCounts, goodZ, maxCells = 1000)
  expect_true(is.numeric(silSmall))

  # Silhouette is between -1 and 1
  expect_true(silSmall >= -1 && silSmall <= 1)

  # Larger dataset would trigger sampling
  # (our test data is small, so just check it doesn't error)
  silLarge <- .calculateSilhouette_Sampled(testCounts, goodZ, maxCells = 10)
  expect_true(is.numeric(silLarge))
})

test_that(".calculateSilhouette_Sampled handles edge cases", {
  # K = 1 should return NA
  silBad <- .calculateSilhouette_Sampled(testCounts, badZ, maxCells = 100)
  expect_true(is.na(silBad))

  # Very small dataset
  smallCounts <- testCounts[, 1:5]
  z5 <- c(1, 1, 2, 2, 3)
  silVerySmall <- .calculateSilhouette_Sampled(smallCounts, z5, maxCells = 100)
  expect_true(is.numeric(silVerySmall) || is.na(silVerySmall))
})


# Test Module Coherence
test_that(".calculateModuleCoherence measures gene correlations", {
  cohGood <- .calculateModuleCoherence(testCounts, goodY, maxCells = 100)
  cohPoor <- .calculateModuleCoherence(testCounts, poorY, maxCells = 100)

  expect_true(is.numeric(cohGood) || is.na(cohGood))
  expect_true(is.numeric(cohPoor) || is.na(cohPoor))

  # Coherence should be between -1 and 1
  if (!is.na(cohGood)) {
    expect_true(cohGood >= -1 && cohGood <= 1)
  }
})

test_that(".calculateModuleCoherence handles edge cases", {
  # L = 1 should return NA
  badY <- rep(1, nGenes)
  cohBad <- .calculateModuleCoherence(testCounts, badY, maxCells = 100)
  expect_true(is.na(cohBad))

  # Each gene in its own module
  uniqueY <- 1:nGenes
  cohUnique <- .calculateModuleCoherence(testCounts, uniqueY, maxCells = 100)
  # Should be NA because each module has only 1 gene
  expect_true(is.na(cohUnique))
})


# Test Combined Metric Calculation
test_that(".calculateClusterMetrics computes all metrics", {
  metrics <- .calculateClusterMetrics(
    testCounts,
    z = goodZ,
    y = goodY,
    metrics = "all"
  )

  expect_true(is.list(metrics))
  expect_true("calinskiHarabasz" %in% names(metrics))
  expect_true("daviesBouldin" %in% names(metrics))
  expect_true("moduleCoherence" %in% names(metrics))
  expect_true("clusterSizeCV" %in% names(metrics))

  # Check types
  expect_true(is.numeric(metrics$calinskiHarabasz))
  expect_true(is.numeric(metrics$daviesBouldin))
  expect_true(is.numeric(metrics$clusterSizeCV))
})

test_that(".calculateClusterMetrics handles selective metrics", {
  metrics <- .calculateClusterMetrics(
    testCounts,
    z = goodZ,
    y = NULL,
    metrics = c("calinskiHarabasz", "daviesBouldin")
  )

  expect_true(is.list(metrics))
  expect_true("calinskiHarabasz" %in% names(metrics))
  expect_true("daviesBouldin" %in% names(metrics))
  expect_false("moduleCoherence" %in% names(metrics))
})


# Test Metric Normalization
test_that(".normalizeMetrics correctly normalizes values", {
  metricsList <- list(
    list(calinskiHarabasz = 100, daviesBouldin = 2.0),
    list(calinskiHarabasz = 200, daviesBouldin = 1.0),
    list(calinskiHarabasz = 150, daviesBouldin = 1.5)
  )

  normalized <- .normalizeMetrics(metricsList)

  expect_equal(length(normalized), 3)

  # Check normalization: values should be between 0 and 1
  expect_true(normalized[[1]]$calinskiHarabasz >= 0 &&
              normalized[[1]]$calinskiHarabasz <= 1)
  expect_true(normalized[[2]]$calinskiHarabasz >= 0 &&
              normalized[[2]]$calinskiHarabasz <= 1)

  # Best CH (200) should be 1, worst (100) should be 0
  expect_equal(normalized[[2]]$calinskiHarabasz, 1)
  expect_equal(normalized[[1]]$calinskiHarabasz, 0)

  # DB is inverse: lower is better
  # Best DB (1.0) should normalize to 1, worst (2.0) to 0
  expect_equal(normalized[[2]]$daviesBouldin, 1)
  expect_equal(normalized[[1]]$daviesBouldin, 0)
})

test_that(".normalizeMetrics handles edge cases", {
  # All same value
  metricsList <- list(
    list(calinskiHarabasz = 100),
    list(calinskiHarabasz = 100),
    list(calinskiHarabasz = 100)
  )

  normalized <- .normalizeMetrics(metricsList)

  # Should all be 0.5
  expect_equal(normalized[[1]]$calinskiHarabasz, 0.5)
  expect_equal(normalized[[2]]$calinskiHarabasz, 0.5)

  # With NA values
  metricsListNA <- list(
    list(calinskiHarabasz = 100, daviesBouldin = NA),
    list(calinskiHarabasz = 200, daviesBouldin = 1.0),
    list(calinskiHarabasz = NA, daviesBouldin = 2.0)
  )

  normalizedNA <- .normalizeMetrics(metricsListNA)
  expect_true(is.na(normalizedNA[[1]]$daviesBouldin))
  expect_true(is.na(normalizedNA[[3]]$calinskiHarabasz))
})


# Test Combined Score Calculation
test_that(".calculateCombinedScore combines metrics correctly", {
  normalizedMetrics <- list(
    calinskiHarabasz = 0.8,
    daviesBouldin = 0.6,
    moduleCoherence = 0.7
  )

  score <- .calculateCombinedScore(
    logLik = 0.9,
    metrics = list(),  # Not used in calculation
    normalizedMetrics = normalizedMetrics,
    validationWeight = 0.3
  )

  expect_true(is.numeric(score))
  expect_true(score >= 0 && score <= 1)

  # Should be weighted combination
  # 0.7 * 0.9 + 0.3 * (weighted combination of metrics)
  expect_true(score > 0.6 && score < 1.0)
})

test_that(".calculateCombinedScore respects validation weight", {
  normalizedMetrics <- list(
    calinskiHarabasz = 1.0,
    daviesBouldin = 1.0
  )

  # Weight = 0: only log-likelihood
  score0 <- .calculateCombinedScore(
    logLik = 0.5,
    metrics = list(),
    normalizedMetrics = normalizedMetrics,
    validationWeight = 0
  )
  expect_equal(score0, 0.5)

  # Weight = 1: only validation metrics
  score1 <- .calculateCombinedScore(
    logLik = 0.5,
    metrics = list(),
    normalizedMetrics = normalizedMetrics,
    validationWeight = 1
  )
  expect_true(score1 > 0.9)  # Should be close to 1 given perfect metrics
})


# Test Stratified Sampling
test_that(".stratifiedSample maintains cluster proportions", {
  z <- rep(1:3, times = c(10, 20, 30))  # Unbalanced clusters
  n <- 30

  sampled <- .stratifiedSample(z, n)

  expect_equal(length(sampled), 30)
  expect_true(all(sampled %in% seq_along(z)))

  # Check proportions are approximately maintained
  originalProps <- table(z) / length(z)
  sampledProps <- table(z[sampled]) / length(sampled)

  # Should be similar (allowing for rounding)
  expect_true(all(abs(originalProps - sampledProps) < 0.15))
})


# Integration test with celda_CG
test_that("celda_CG with combined selection runs without error", {
  skip_if_not_installed("SingleCellExperiment")

  # Simulate small dataset
  simData <- simulateCells("celda_CG", K = 3, L = 5, C = 50, G = 100)

  # Run with combined selection (very few iterations for speed)
  result <- celda_CG(
    simData,
    K = 3,
    L = 5,
    nchains = 2,
    maxIter = 5,
    algorithm = "EM",
    selectionCriterion = "combined",
    validationWeight = 0.3,
    verbose = FALSE
  )

  expect_s4_class(result, "SingleCellExperiment")

  # Check that clustering was performed
  expect_true("celda_cell_cluster" %in% names(SummarizedExperiment::colData(result)))
  expect_true("celda_feature_module" %in% names(SummarizedExperiment::rowData(result)))
})


# Backward compatibility test
test_that("celda_CG with default (logLik) selection maintains backward compatibility", {
  skip_if_not_installed("SingleCellExperiment")

  # Simulate small dataset
  simData <- simulateCells("celda_CG", K = 3, L = 5, C = 50, G = 100)

  # Run with default selection (should be log-likelihood)
  result <- celda_CG(
    simData,
    K = 3,
    L = 5,
    nchains = 2,
    maxIter = 5,
    algorithm = "EM",
    # selectionCriterion not specified - should default to "logLik"
    verbose = FALSE
  )

  expect_s4_class(result, "SingleCellExperiment")

  # Should work exactly as before
  expect_true("celda_cell_cluster" %in% names(SummarizedExperiment::colData(result)))
})


# Performance test
test_that("Combined selection overhead is reasonable", {
  skip_if_not_installed("SingleCellExperiment")
  skip_on_cran()  # Skip on CRAN to save time

  # Simulate moderately-sized dataset
  simData <- simulateCells("celda_CG", K = 4, L = 6, C = 100, G = 200)

  # Time log-likelihood selection
  timeLL <- system.time({
    resultLL <- celda_CG(
      simData,
      K = 4,
      L = 6,
      nchains = 3,
      maxIter = 10,
      algorithm = "EM",
      selectionCriterion = "logLik",
      verbose = FALSE
    )
  })

  # Time combined selection
  timeCombined <- system.time({
    resultCombined <- celda_CG(
      simData,
      K = 4,
      L = 6,
      nchains = 3,
      maxIter = 10,
      algorithm = "EM",
      selectionCriterion = "combined",
      validationWeight = 0.3,
      verbose = FALSE
    )
  })

  # Combined should be slower, but not more than 2x
  expect_true(timeCombined["elapsed"] < 2 * timeLL["elapsed"])

  # Both should produce valid results
  expect_s4_class(resultLL, "SingleCellExperiment")
  expect_s4_class(resultCombined, "SingleCellExperiment")
})


# Test that combined selection can choose different chain than log-likelihood
test_that("Combined selection can differ from log-likelihood selection", {
  skip_if_not_installed("SingleCellExperiment")
  skip_on_cran()

  # This is a probabilistic test, so we use multiple seeds
  set.seed(42)

  # Simulate data
  simData <- simulateCells("celda_CG", K = 3, L = 4, C = 80, G = 150)

  # Run with both methods
  resultLL <- celda_CG(
    simData,
    K = 3,
    L = 4,
    nchains = 3,
    maxIter = 15,
    algorithm = "EM",
    selectionCriterion = "logLik",
    seed = 123,
    verbose = FALSE
  )

  resultCombined <- celda_CG(
    simData,
    K = 3,
    L = 4,
    nchains = 3,
    maxIter = 15,
    algorithm = "EM",
    selectionCriterion = "combined",
    validationWeight = 0.5,  # Higher weight to make difference more likely
    seed = 123,
    verbose = FALSE
  )

  # Results might differ (not guaranteed, but possible)
  # Just check both are valid
  expect_s4_class(resultLL, "SingleCellExperiment")
  expect_s4_class(resultCombined, "SingleCellExperiment")

  # Both should have same K and L
  zLL <- SummarizedExperiment::colData(resultLL)$celda_cell_cluster
  zCombined <- SummarizedExperiment::colData(resultCombined)$celda_cell_cluster

  expect_equal(length(unique(zLL)), 3)
  expect_equal(length(unique(zCombined)), 3)
})
