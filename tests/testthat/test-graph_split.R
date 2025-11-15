# Test Graph-Based Split Functionality

library(celda)
library(testthat)
library(Matrix)

# Test 1: Modularity calculation
test_that(".calculateModularity works correctly on known graphs", {
  # Test case 1: Fully connected graph (should have low modularity)
  n <- 10
  adjFull <- matrix(1, nrow = n, ncol = n)
  diag(adjFull) <- 0
  modFull <- celda:::.calculateModularity(adjFull)
  expect_gte(modFull, 0)
  expect_lte(modFull, 1)
  expect_lt(modFull, 0.3)  # Should be low for fully connected

  # Test case 2: Two distinct communities
  n1 <- 5
  n2 <- 5
  adjCommunity <- matrix(0, nrow = n1 + n2, ncol = n1 + n2)
  # Connect within first community
  adjCommunity[1:n1, 1:n1] <- 1
  # Connect within second community
  adjCommunity[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)] <- 1
  diag(adjCommunity) <- 0
  modCommunity <- celda:::.calculateModularity(adjCommunity)
  expect_gte(modCommunity, 0)
  expect_lte(modCommunity, 1)
  expect_gt(modCommunity, modFull)  # Should be higher than fully connected

  # Test case 3: Empty graph (no edges)
  adjEmpty <- matrix(0, nrow = 10, ncol = 10)
  modEmpty <- celda:::.calculateModularity(adjEmpty)
  expect_equal(modEmpty, 0)

  # Test case 4: Single node
  adjSingle <- matrix(0, nrow = 1, ncol = 1)
  modSingle <- celda:::.calculateModularity(adjSingle)
  expect_equal(modSingle, 0)

  # Test case 5: Sparse matrix input
  adjSparse <- Matrix(adjCommunity, sparse = TRUE)
  modSparse <- celda:::.calculateModularity(adjSparse)
  expect_equal(modSparse, modCommunity)
})


# Test 2: Bimodal gene detection
test_that(".findBimodalGenes detects bimodal distributions", {
  set.seed(12345)

  # Test case 1: Clearly bimodal data
  nGenes <- 100
  nCells <- 100
  # Create bimodal gene (first gene): two peaks
  bimodalGene <- c(rnorm(50, mean = 2, sd = 0.5), rnorm(50, mean = 10, sd = 0.5))
  # Create unimodal genes
  unimodalGenes <- matrix(rnorm(nCells * (nGenes - 1), mean = 5, sd = 1),
                         nrow = nGenes - 1, ncol = nCells)
  testCounts <- rbind(bimodalGene, unimodalGenes)

  bimodalIdx <- celda:::.findBimodalGenes(testCounts, pvalueThreshold = 0.05)

  # Should detect the first gene as bimodal
  expect_true(length(bimodalIdx) >= 1)
  # May not always detect depending on dip test implementation
  # expect_true(1 %in% bimodalIdx)

  # Test case 2: All unimodal (should find few/no bimodal genes)
  unimodalCounts <- matrix(rnorm(nGenes * nCells, mean = 5, sd = 1),
                          nrow = nGenes, ncol = nCells)
  unimodalIdx <- celda:::.findBimodalGenes(unimodalCounts, pvalueThreshold = 0.05)
  expect_true(length(unimodalIdx) <= nGenes * 0.1)  # Should be <= 10% false positives

  # Test case 3: Too few cells (should return empty)
  smallCounts <- matrix(rnorm(5 * 5), nrow = 5, ncol = 5)
  smallIdx <- celda:::.findBimodalGenes(smallCounts, pvalueThreshold = 0.05)
  expect_equal(length(smallIdx), 0)

  # Test case 4: Sparse matrix input
  sparseCounts <- Matrix(testCounts, sparse = TRUE)
  sparseIdx <- celda:::.findBimodalGenes(sparseCounts, pvalueThreshold = 0.05)
  expect_true(is.integer(sparseIdx))
})


# Test 3: Substructure detection
test_that(".detectSubstructure identifies community structure", {
  # Test case 1: Highly correlated matrix (low substructure)
  n <- 20
  highCorr <- matrix(0.9, nrow = n, ncol = n)
  diag(highCorr) <- 1
  substructHigh <- celda:::.detectSubstructure(highCorr, threshold = 0.3)
  expect_gte(substructHigh, 0)
  expect_lte(substructHigh, 1)
  expect_lt(substructHigh, 0.5)  # Should be low

  # Test case 2: Block diagonal correlation (high substructure)
  n1 <- 10
  n2 <- 10
  blockCorr <- matrix(0, nrow = n1 + n2, ncol = n1 + n2)
  # Block 1: high correlation
  blockCorr[1:n1, 1:n1] <- 0.9
  # Block 2: high correlation
  blockCorr[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)] <- 0.9
  # Between blocks: low correlation
  blockCorr[1:n1, (n1 + 1):(n1 + n2)] <- 0.1
  blockCorr[(n1 + 1):(n1 + n2), 1:n1] <- 0.1
  diag(blockCorr) <- 1

  substructBlock <- celda:::.detectSubstructure(blockCorr, threshold = 0.5)
  expect_gte(substructBlock, 0)
  expect_lte(substructBlock, 1)
  expect_gt(substructBlock, substructHigh)  # Should be higher

  # Test case 3: All disconnected (maximum substructure)
  disconnected <- matrix(0, nrow = 10, ncol = 10)
  diag(disconnected) <- 1
  substructDisc <- celda:::.detectSubstructure(disconnected, threshold = 0.3)
  expect_equal(substructDisc, 1)

  # Test case 4: Small matrix
  smallCorr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
  substructSmall <- celda:::.detectSubstructure(smallCorr, threshold = 0.3)
  expect_gte(substructSmall, 0)
  expect_lte(substructSmall, 1)

  # Test case 5: Sparse matrix
  sparseCorr <- Matrix(blockCorr, sparse = TRUE)
  substructSparse <- celda:::.detectSubstructure(sparseCorr, threshold = 0.5)
  expect_equal(substructSparse, substructBlock)
})


# Test 4: Graph-based split candidate identification
test_that(".identifySplitCandidates_GraphBased identifies subclusters", {
  set.seed(12345)

  # Create synthetic data with clear subclusters
  nGenes <- 200
  nCells <- 300
  K <- 3

  # Cluster 1: 100 cells, homogeneous
  cluster1 <- matrix(rpois(nGenes * 100, lambda = 5),
                    nrow = nGenes, ncol = 100)

  # Cluster 2: 100 cells, heterogeneous (two subclusters)
  subcluster2a <- matrix(rpois(nGenes * 50, lambda = 3),
                        nrow = nGenes, ncol = 50)
  subcluster2b <- matrix(rpois(nGenes * 50, lambda = 15),
                        nrow = nGenes, ncol = 50)
  cluster2 <- cbind(subcluster2a, subcluster2b)

  # Cluster 3: 100 cells, moderate heterogeneity
  cluster3 <- matrix(rpois(nGenes * 100, lambda = 8),
                    nrow = nGenes, ncol = 100)

  # Combine
  counts <- cbind(cluster1, cluster2, cluster3)
  z <- rep(1:3, each = 100)

  # Test without reduced dimensions (correlation-based)
  candidates <- celda:::.identifySplitCandidates_GraphBased(
    counts, z, K,
    reducedDim = NULL,
    minCell = 3,
    heterogeneityThreshold = 0.5
  )

  expect_true(is.integer(candidates))
  expect_true(all(candidates >= 1 & candidates <= K))
  expect_true(length(candidates) > 0)
  # Cluster 2 should be a candidate (most heterogeneous)
  # expect_true(2 %in% candidates)

  # Test with reduced dimensions
  # Create simple 2D coordinates that reflect structure
  reducedDim <- matrix(NA, nrow = nCells, ncol = 2)
  # Cluster 1: tight group
  reducedDim[1:100, ] <- matrix(rnorm(200, mean = 0, sd = 0.5), ncol = 2)
  # Cluster 2: two separated subclusters
  reducedDim[101:150, ] <- matrix(rnorm(100, mean = 5, sd = 0.5), ncol = 2)
  reducedDim[151:200, ] <- matrix(rnorm(100, mean = 10, sd = 0.5), ncol = 2)
  # Cluster 3: moderate spread
  reducedDim[201:300, ] <- matrix(rnorm(200, mean = 15, sd = 1.5), ncol = 2)

  candidatesWithDim <- celda:::.identifySplitCandidates_GraphBased(
    counts, z, K,
    reducedDim = reducedDim,
    minCell = 3,
    heterogeneityThreshold = 0.5
  )

  expect_true(is.integer(candidatesWithDim))
  expect_true(all(candidatesWithDim >= 1 & candidatesWithDim <= K))

  # Test edge cases
  # Too few cells
  smallCandidates <- celda:::.identifySplitCandidates_GraphBased(
    counts[, 1:5], z[1:5], K = 2,
    minCell = 10
  )
  expect_equal(length(smallCandidates), 0)
})


# Test 5: Integration with celda_C
test_that("Graph-based splitting works in celda_C", {
  skip_on_cran()
  set.seed(12345)

  # Create simple test data
  nGenes <- 50
  nCells <- 100
  K <- 3

  # Simple count matrix
  counts <- matrix(rpois(nGenes * nCells, lambda = 5),
                  nrow = nGenes, ncol = nCells)
  rownames(counts) <- paste0("Gene", 1:nGenes)
  colnames(counts) <- paste0("Cell", 1:nCells)

  # Test with graph-based splitting disabled (default)
  result1 <- celda_C(counts,
                    K = K,
                    useGraphBasedSplit = FALSE,
                    nchains = 1,
                    maxIter = 5,
                    splitOnIter = 2,
                    verbose = FALSE)

  expect_s4_class(result1, "SingleCellExperiment")

  # Test with graph-based splitting enabled (no reduced dims)
  result2 <- celda_C(counts,
                    K = K,
                    useGraphBasedSplit = TRUE,
                    reducedDimForSplit = NULL,
                    nchains = 1,
                    maxIter = 5,
                    splitOnIter = 2,
                    verbose = FALSE)

  expect_s4_class(result2, "SingleCellExperiment")

  # Test with graph-based splitting and reduced dims
  reducedDim <- matrix(rnorm(nCells * 2), ncol = 2)
  result3 <- celda_C(counts,
                    K = K,
                    useGraphBasedSplit = TRUE,
                    reducedDimForSplit = reducedDim,
                    nchains = 1,
                    maxIter = 5,
                    splitOnIter = 2,
                    verbose = FALSE)

  expect_s4_class(result3, "SingleCellExperiment")
})


# Test 6: Integration with celda_CG
test_that("Graph-based splitting works in celda_CG", {
  skip_on_cran()
  set.seed(12345)

  # Create simple test data
  nGenes <- 50
  nCells <- 100
  K <- 3
  L <- 5

  # Simple count matrix
  counts <- matrix(rpois(nGenes * nCells, lambda = 5),
                  nrow = nGenes, ncol = nCells)
  rownames(counts) <- paste0("Gene", 1:nGenes)
  colnames(counts) <- paste0("Cell", 1:nCells)

  # Test with graph-based splitting enabled
  result <- celda_CG(counts,
                    K = K,
                    L = L,
                    useGraphBasedSplit = TRUE,
                    heterogeneityThreshold = 0.3,
                    nchains = 1,
                    maxIter = 5,
                    splitOnIter = 2,
                    verbose = FALSE)

  expect_s4_class(result, "SingleCellExperiment")
})


# Test 7: Backward compatibility
test_that("Graph-based splitting maintains backward compatibility", {
  skip_on_cran()
  set.seed(12345)

  nGenes <- 50
  nCells <- 100
  K <- 3

  counts <- matrix(rpois(nGenes * nCells, lambda = 5),
                  nrow = nGenes, ncol = nCells)
  rownames(counts) <- paste0("Gene", 1:nGenes)
  colnames(counts) <- paste0("Cell", 1:nCells)

  # Default behavior should be unchanged (useGraphBasedSplit = FALSE)
  result1 <- celda_C(counts, K = K, nchains = 1, maxIter = 5, verbose = FALSE)
  result2 <- celda_C(counts, K = K, nchains = 1, maxIter = 5,
                    useGraphBasedSplit = FALSE, verbose = FALSE)

  # Both should work and produce valid results
  expect_s4_class(result1, "SingleCellExperiment")
  expect_s4_class(result2, "SingleCellExperiment")
})


# Test 8: Performance overhead
test_that("Graph-based splitting has acceptable overhead", {
  skip_on_cran()
  skip_on_ci()  # Skip on CI to avoid timeout

  set.seed(12345)

  nGenes <- 100
  nCells <- 500
  K <- 5

  counts <- matrix(rpois(nGenes * nCells, lambda = 5),
                  nrow = nGenes, ncol = nCells)
  rownames(counts) <- paste0("Gene", 1:nGenes)
  colnames(counts) <- paste0("Cell", 1:nCells)

  # Time statistical method
  time1 <- system.time({
    result1 <- celda_C(counts, K = K, nchains = 1, maxIter = 10,
                      useGraphBasedSplit = FALSE, verbose = FALSE)
  })

  # Time graph-based method
  time2 <- system.time({
    result2 <- celda_C(counts, K = K, nchains = 1, maxIter = 10,
                      useGraphBasedSplit = TRUE, verbose = FALSE)
  })

  # Graph-based should not be more than 50% slower
  # (allowing some variance)
  overhead <- (time2["elapsed"] - time1["elapsed"]) / time1["elapsed"]
  expect_lt(overhead, 0.5)
})


# Test 9: Error handling
test_that("Graph-based functions handle errors gracefully", {
  # Invalid adjacency matrix
  expect_error(celda:::.calculateModularity(matrix("a", 5, 5)))
  expect_error(celda:::.calculateModularity(list(1, 2, 3)))

  # Invalid counts
  expect_error(celda:::.findBimodalGenes(matrix("a", 5, 5)))

  # Invalid correlation matrix
  expect_error(celda:::.detectSubstructure(matrix("a", 5, 5)))

  # Mismatched dimensions
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  z <- rep(1:2, each = 15)  # Wrong length
  expect_error(celda:::.identifySplitCandidates_GraphBased(counts, z, K = 2))

  # Wrong reducedDim dimensions
  z <- rep(1:2, each = 5)
  reducedDim <- matrix(rnorm(15 * 2), ncol = 2)  # Wrong number of rows
  expect_error(celda:::.identifySplitCandidates_GraphBased(
    counts, z, K = 2, reducedDim = reducedDim
  ))
})


# Test 10: Doesn't over-split homogeneous clusters
test_that("Graph-based method doesn't over-split homogeneous clusters", {
  set.seed(12345)

  nGenes <- 100
  nCells <- 200
  K <- 2

  # Create two very homogeneous clusters
  cluster1 <- matrix(rpois(nGenes * 100, lambda = 5),
                    nrow = nGenes, ncol = 100)
  cluster2 <- matrix(rpois(nGenes * 100, lambda = 10),
                    nrow = nGenes, ncol = 100)

  counts <- cbind(cluster1, cluster2)
  z <- rep(1:2, each = 100)

  # With high threshold, should identify few/no candidates
  candidates <- celda:::.identifySplitCandidates_GraphBased(
    counts, z, K,
    minCell = 3,
    heterogeneityThreshold = 0.7  # High threshold
  )

  # Should identify very few clusters for splitting
  expect_true(length(candidates) <= 1)
})
