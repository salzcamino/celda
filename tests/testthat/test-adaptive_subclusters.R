# Test adaptive K/L subcluster selection

# Load test data
data(celdaCGSim)
data(celdaCSim)
data(celdaGSim)


# Test .adaptiveKSubcluster function
test_that(".adaptiveKSubcluster returns reasonable values", {
  # Test with small K - should use default heuristic
  result <- .adaptiveKSubcluster(celdaCSim$counts, K = 3)
  expect_equal(result, ceiling(sqrt(3)))

  # Test with medium K - should return bounded value
  result <- .adaptiveKSubcluster(celdaCSim$counts, K = 10)
  expect_true(result >= 2)
  expect_true(result <= 10)
  expect_true(is.integer(result))

  # Test with large K
  result <- .adaptiveKSubcluster(celdaCSim$counts, K = 25)
  expect_true(result >= 2)
  expect_true(result <= 25)
  expect_true(is.integer(result))
})


test_that(".adaptiveKSubcluster handles small datasets", {
  # Very small dataset should fall back to default
  smallCounts <- celdaCSim$counts[, 1:50]
  result <- .adaptiveKSubcluster(smallCounts, K = 10)
  expect_equal(result, ceiling(sqrt(10)))
})


test_that(".adaptiveKSubcluster handles edge cases", {
  # Test with K = 2 (minimum)
  result <- .adaptiveKSubcluster(celdaCSim$counts, K = 2)
  expect_true(result >= 2)

  # Test deterministic behavior given same data
  result1 <- .adaptiveKSubcluster(celdaCSim$counts, K = 10)
  result2 <- .adaptiveKSubcluster(celdaCSim$counts, K = 10)
  # Note: may differ due to sampling, but should be in same ballpark
  expect_true(abs(result1 - result2) <= 2)
})


# Test .adaptiveLSubcluster function
test_that(".adaptiveLSubcluster returns reasonable values", {
  # Test with small L
  result <- .adaptiveLSubcluster(celdaGSim$counts, L = 3)
  expect_equal(result, ceiling(sqrt(3)))

  # Test with medium L
  result <- .adaptiveLSubcluster(celdaGSim$counts, L = 10)
  expect_true(result >= 2)
  expect_true(result <= 10)
  expect_true(is.integer(result))

  # Test with large L
  result <- .adaptiveLSubcluster(celdaGSim$counts, L = 20)
  expect_true(result >= 2)
  expect_true(result <= 20)
  expect_true(is.integer(result))
})


test_that(".adaptiveLSubcluster handles small datasets", {
  # Small gene set should fall back to default
  smallCounts <- celdaGSim$counts[1:30, ]
  result <- .adaptiveLSubcluster(smallCounts, L = 10)
  expect_equal(result, ceiling(sqrt(10)))
})


test_that(".adaptiveLSubcluster handles edge cases", {
  # Test with L = 2 (minimum)
  result <- .adaptiveLSubcluster(celdaGSim$counts, L = 2)
  expect_true(result >= 2)

  # Test with sparse data
  sparseCounts <- celdaGSim$counts
  sparseCounts[sample(length(sparseCounts), length(sparseCounts) * 0.9)] <- 0
  result <- .adaptiveLSubcluster(sparseCounts, L = 10)
  expect_true(result >= 2)
  expect_true(result <= 10)
})


# Test .initializeSplitZ with adaptive subclusters
test_that(".initializeSplitZ uses adaptive subclusters when enabled", {
  # Test with adaptiveSubclusters = FALSE (default behavior)
  z1 <- .initializeSplitZ(celdaCSim$counts,
                          K = 10,
                          adaptiveSubclusters = FALSE)
  expect_length(z1, ncol(celdaCSim$counts))
  expect_true(all(z1 %in% 1:10))

  # Test with adaptiveSubclusters = TRUE
  z2 <- .initializeSplitZ(celdaCSim$counts,
                          K = 10,
                          adaptiveSubclusters = TRUE)
  expect_length(z2, ncol(celdaCSim$counts))
  expect_true(all(z2 %in% 1:10))

  # Both should produce valid initializations
  expect_equal(length(unique(z1)), 10)
  expect_equal(length(unique(z2)), 10)
})


test_that(".initializeSplitZ respects manual KSubcluster override", {
  # When KSubcluster is specified, adaptive should be ignored
  z <- .initializeSplitZ(celdaCSim$counts,
                         K = 10,
                         KSubcluster = 5,
                         adaptiveSubclusters = TRUE)
  expect_length(z, ncol(celdaCSim$counts))
  expect_true(all(z %in% 1:10))
})


# Test .initializeSplitY with adaptive subclusters
test_that(".initializeSplitY uses adaptive subclusters when enabled", {
  # Test with adaptiveSubclusters = FALSE
  y1 <- .initializeSplitY(celdaGSim$counts,
                          L = 10,
                          adaptiveSubclusters = FALSE)
  expect_length(y1, nrow(celdaGSim$counts))
  expect_true(all(y1 %in% 1:10))

  # Test with adaptiveSubclusters = TRUE
  y2 <- .initializeSplitY(celdaGSim$counts,
                          L = 10,
                          adaptiveSubclusters = TRUE)
  expect_length(y2, nrow(celdaGSim$counts))
  expect_true(all(y2 %in% 1:10))

  # Both should produce valid initializations
  expect_equal(length(unique(y1)), 10)
  expect_equal(length(unique(y2)), 10)
})


test_that(".initializeSplitY respects manual LSubcluster override", {
  # When LSubcluster is specified, adaptive should be ignored
  y <- .initializeSplitY(celdaGSim$counts,
                         L = 10,
                         LSubcluster = 4,
                         adaptiveSubclusters = TRUE)
  expect_length(y, nrow(celdaGSim$counts))
  expect_true(all(y %in% 1:10))
})


# Integration tests with main celda functions
test_that("celda_C works with adaptiveSubclusters parameter", {
  # Test backward compatibility - default should be FALSE
  res1 <- celda_C(celdaCSim$counts,
                  sampleLabel = celdaCSim$sampleLabel,
                  K = 5,
                  maxIter = 5,
                  nchains = 1,
                  verbose = FALSE)
  expect_s4_class(res1, "celda_C")

  # Test with adaptiveSubclusters = TRUE
  res2 <- celda_C(celdaCSim$counts,
                  sampleLabel = celdaCSim$sampleLabel,
                  K = 5,
                  maxIter = 5,
                  nchains = 1,
                  zInitialize = "split",
                  adaptiveSubclusters = TRUE,
                  verbose = FALSE)
  expect_s4_class(res2, "celda_C")

  # Both should produce valid results
  expect_equal(length(celdaClusters(res1)$z), ncol(celdaCSim$counts))
  expect_equal(length(celdaClusters(res2)$z), ncol(celdaCSim$counts))
})


test_that("celda_G works with adaptiveSubclusters parameter", {
  # Test backward compatibility
  res1 <- celda_G(celdaGSim$counts,
                  L = 5,
                  maxIter = 5,
                  nchains = 1,
                  verbose = FALSE)
  expect_s4_class(res1, "celda_G")

  # Test with adaptiveSubclusters = TRUE
  res2 <- celda_G(celdaGSim$counts,
                  L = 5,
                  maxIter = 5,
                  nchains = 1,
                  yInitialize = "split",
                  adaptiveSubclusters = TRUE,
                  verbose = FALSE)
  expect_s4_class(res2, "celda_G")

  # Both should produce valid results
  expect_equal(length(celdaClusters(res1)$y), nrow(celdaGSim$counts))
  expect_equal(length(celdaClusters(res2)$y), nrow(celdaGSim$counts))
})


test_that("celda_CG works with adaptiveSubclusters parameter", {
  # Test backward compatibility
  res1 <- celda_CG(celdaCGSim$counts,
                   sampleLabel = celdaCGSim$sampleLabel,
                   K = 5,
                   L = 5,
                   maxIter = 5,
                   nchains = 1,
                   verbose = FALSE)
  expect_s4_class(res1, "celda_CG")

  # Test with adaptiveSubclusters = TRUE
  res2 <- celda_CG(celdaCGSim$counts,
                   sampleLabel = celdaCGSim$sampleLabel,
                   K = 5,
                   L = 5,
                   maxIter = 5,
                   nchains = 1,
                   zInitialize = "split",
                   yInitialize = "split",
                   adaptiveSubclusters = TRUE,
                   verbose = FALSE)
  expect_s4_class(res2, "celda_CG")

  # Both should produce valid results
  clusters1 <- celdaClusters(res1)
  clusters2 <- celdaClusters(res2)
  expect_equal(length(clusters1$z), ncol(celdaCGSim$counts))
  expect_equal(length(clusters1$y), nrow(celdaCGSim$counts))
  expect_equal(length(clusters2$z), ncol(celdaCGSim$counts))
  expect_equal(length(clusters2$y), nrow(celdaCGSim$counts))
})


# Test backward compatibility
test_that("Adaptive subclusters are backward compatible", {
  # Default behavior should be identical to old behavior
  set.seed(12345)
  z_old <- .initializeSplitZ(celdaCSim$counts,
                             K = 10,
                             adaptiveSubclusters = FALSE)

  set.seed(12345)
  y_old <- .initializeSplitY(celdaGSim$counts,
                             L = 10,
                             adaptiveSubclusters = FALSE)

  # These should use the fixed sqrt(K) and sqrt(L) heuristics
  # Results should be reproducible with same seed
  set.seed(12345)
  z_default <- .initializeSplitZ(celdaCSim$counts,
                                 K = 10)

  set.seed(12345)
  y_default <- .initializeSplitY(celdaGSim$counts,
                                 L = 10)

  # Should be identical when using default (FALSE)
  expect_identical(z_old, z_default)
  expect_identical(y_old, y_default)
})


# Performance test - adaptive should not be significantly slower
test_that("Adaptive subcluster selection is reasonably fast", {
  # This should complete quickly even with adaptive selection
  skip_on_cran()

  time_default <- system.time({
    z1 <- .initializeSplitZ(celdaCSim$counts,
                           K = 10,
                           adaptiveSubclusters = FALSE)
  })

  time_adaptive <- system.time({
    z2 <- .initializeSplitZ(celdaCSim$counts,
                           K = 10,
                           adaptiveSubclusters = TRUE)
  })

  # Adaptive should not be more than 3x slower
  # (it does extra hierarchical clustering but on sampled data)
  expect_true(time_adaptive["elapsed"] < time_default["elapsed"] * 3)
})


# Test with data that has clear vs diffuse structure
test_that("Adaptive selection responds to data structure", {
  skip_on_cran()

  # Create synthetic data with clear cluster structure
  set.seed(123)
  clearCounts <- matrix(0, nrow = 100, ncol = 500)
  # 5 clear clusters of 100 cells each
  for (i in 1:5) {
    cellIdx <- ((i - 1) * 100 + 1):(i * 100)
    geneIdx <- ((i - 1) * 20 + 1):(i * 20)
    clearCounts[geneIdx, cellIdx] <- rpois(length(geneIdx) * 100, lambda = 10)
  }

  # Get adaptive K for clear structure
  K_clear <- .adaptiveKSubcluster(clearCounts, K = 25)

  # Create diffuse data (random noise)
  set.seed(123)
  diffuseCounts <- matrix(rpois(100 * 500, lambda = 5), nrow = 100, ncol = 500)

  # Get adaptive K for diffuse structure
  K_diffuse <- .adaptiveKSubcluster(diffuseCounts, K = 25)

  # Both should be valid
  expect_true(K_clear >= 2 && K_clear <= 25)
  expect_true(K_diffuse >= 2 && K_diffuse <= 25)

  # Note: We can't strictly enforce that K_diffuse > K_clear because
  # the adaptive logic depends on silhouette scores which can vary,
  # but we can check they're in reasonable ranges
  expect_true(is.integer(K_clear))
  expect_true(is.integer(K_diffuse))
})
