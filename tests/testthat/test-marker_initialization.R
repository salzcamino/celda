# Test marker-guided initialization
library(celda)
library(SingleCellExperiment)
library(Matrix)
context("Testing marker-guided initialization")

# Create a simple test dataset
set.seed(12345)
nGenes <- 100
nCells <- 200
counts <- matrix(rpois(nGenes * nCells, lambda = 5),
    nrow = nGenes,
    ncol = nCells)
rownames(counts) <- paste0("Gene_", seq_len(nGenes))
colnames(counts) <- paste0("Cell_", seq_len(nCells))

# Add some marker genes with distinctive expression patterns
# Create 3 cell types with specific markers
cellType1 <- seq_len(70)
cellType2 <- 71:140
cellType3 <- 141:200

# Marker genes for type 1 (genes 1-5)
counts[1:5, cellType1] <- counts[1:5, cellType1] + 20
# Marker genes for type 2 (genes 6-10)
counts[6:10, cellType2] <- counts[6:10, cellType2] + 20
# Marker genes for type 3 (genes 11-15)
counts[11:15, cellType3] <- counts[11:15, cellType3] + 20

markerList <- list(
    "Type1" = paste0("Gene_", 1:5),
    "Type2" = paste0("Gene_", 6:10),
    "Type3" = paste0("Gene_", 11:15)
)

# Test .calculateMarkerScores
test_that("calculateMarkerScores calculates correct scores", {
    scores <- celda:::.calculateMarkerScores(counts, markerList)

    expect_equal(nrow(scores), ncol(counts))
    expect_equal(ncol(scores), length(markerList))
    expect_equal(colnames(scores), names(markerList))

    # Check that scores sum to 1 for each cell (or 0 if no markers)
    rowSums <- rowSums(scores)
    expect_true(all(abs(rowSums - 1) < 1e-10 | rowSums == 0))

    # Check that Type1 markers have highest scores for type 1 cells
    type1Cells <- apply(scores[cellType1, , drop = FALSE], 1, which.max)
    expect_true(mean(type1Cells == 1) > 0.8)  # At least 80% correct

    # Check that Type2 markers have highest scores for type 2 cells
    type2Cells <- apply(scores[cellType2, , drop = FALSE], 1, which.max)
    expect_true(mean(type2Cells == 2) > 0.8)

    # Check that Type3 markers have highest scores for type 3 cells
    type3Cells <- apply(scores[cellType3, , drop = FALSE], 1, which.max)
    expect_true(mean(type3Cells == 3) > 0.8)
})


test_that("calculateMarkerScores handles missing markers", {
    # Create marker list with some non-existent genes
    markerListPartial <- list(
        "Type1" = c("Gene_1", "NonExistent1", "Gene_2"),
        "Type2" = c("Gene_6", "Gene_7")
    )

    expect_warning(
        scores <- celda:::.calculateMarkerScores(counts, markerListPartial),
        "not found"
    )

    expect_equal(nrow(scores), ncol(counts))
    expect_equal(ncol(scores), length(markerListPartial))
})


test_that("calculateMarkerScores handles all missing markers", {
    markerListMissing <- list(
        "Type1" = c("NonExistent1", "NonExistent2")
    )

    expect_warning(
        scores <- celda:::.calculateMarkerScores(counts, markerListMissing),
        "No markers from set"
    )

    expect_equal(nrow(scores), ncol(counts))
    # All scores should be 0
    expect_true(all(scores == 0))
})


test_that("calculateMarkerScores validates input", {
    # Test with unnamed list
    expect_error(
        celda:::.calculateMarkerScores(counts, list(c("Gene_1", "Gene_2"))),
        "must be a named list"
    )

    # Test with counts without rownames
    countsNoNames <- counts
    rownames(countsNoNames) <- NULL
    expect_error(
        celda:::.calculateMarkerScores(countsNoNames, markerList),
        "must have row names"
    )
})


# Test .initializeSplitZ_MarkerGuided
test_that("initializeSplitZ_MarkerGuided produces K clusters", {
    K <- 3
    z <- celda:::.initializeSplitZ_MarkerGuided(counts, K, markerList)

    expect_equal(length(z), ncol(counts))
    expect_equal(length(unique(z)), K)
    expect_true(all(z >= 1 & z <= K))
})


test_that("initializeSplitZ_MarkerGuided clusters match marker expression", {
    K <- 3
    z <- celda:::.initializeSplitZ_MarkerGuided(counts, K, markerList)

    # Calculate purity of each marker-defined group
    # For cells in cellType1, most should be in the same cluster
    type1Clusters <- z[cellType1]
    type1Purity <- max(table(type1Clusters)) / length(type1Clusters)
    expect_true(type1Purity > 0.7)  # At least 70% purity

    # For cells in cellType2
    type2Clusters <- z[cellType2]
    type2Purity <- max(table(type2Clusters)) / length(type2Clusters)
    expect_true(type2Purity > 0.7)

    # For cells in cellType3
    type3Clusters <- z[cellType3]
    type3Purity <- max(table(type3Clusters)) / length(type3Clusters)
    expect_true(type3Purity > 0.7)
})


test_that("initializeSplitZ_MarkerGuided handles K > number of marker sets", {
    K <- 6
    z <- celda:::.initializeSplitZ_MarkerGuided(counts, K, markerList)

    expect_equal(length(z), ncol(counts))
    expect_true(length(unique(z)) >= 3)  # Should have at least 3 clusters
    expect_true(all(z >= 1))
})


test_that("initializeSplitZ_MarkerGuided handles K < number of marker sets", {
    # Add more marker sets
    markerListLarge <- list(
        "Type1" = paste0("Gene_", 1:2),
        "Type2" = paste0("Gene_", 6:7),
        "Type3" = paste0("Gene_", 11:12),
        "Type4" = paste0("Gene_", 16:17),
        "Type5" = paste0("Gene_", 21:22)
    )

    K <- 3
    z <- celda:::.initializeSplitZ_MarkerGuided(counts, K, markerListLarge)

    expect_equal(length(z), ncol(counts))
    expect_equal(length(unique(z)), K)
})


test_that("initializeSplitZ_MarkerGuided falls back when no markers expressed", {
    # Create counts with very low expression
    lowCounts <- matrix(0, nrow = 10, ncol = 50)
    rownames(lowCounts) <- paste0("Gene_", seq_len(nrow(lowCounts)))

    markerListLow <- list("Type1" = paste0("Gene_", 1:3))

    expect_warning(
        z <- celda:::.initializeSplitZ_MarkerGuided(lowCounts, K = 3, markerListLow),
        "Falling back to random"
    )

    expect_equal(length(z), ncol(lowCounts))
})


# Test .initializeSplitZ_PriorClustering
test_that("initializeSplitZ_PriorClustering uses exact match when K matches", {
    priorZ <- sample(1:3, ncol(counts), replace = TRUE)
    K <- 3

    z <- celda:::.initializeSplitZ_PriorClustering(counts, K, priorZ)

    expect_equal(length(z), ncol(counts))
    expect_equal(length(unique(z)), K)
    # Should be identical after renumbering
    expect_equal(as.integer(as.factor(priorZ)), z)
})


test_that("initializeSplitZ_PriorClustering splits when prior has fewer clusters", {
    priorZ <- sample(1:2, ncol(counts), replace = TRUE)
    K <- 5

    z <- celda:::.initializeSplitZ_PriorClustering(counts, K, priorZ)

    expect_equal(length(z), ncol(counts))
    expect_true(length(unique(z)) >= 2)  # Should have at least original 2
    expect_true(length(unique(z)) <= K)  # But not necessarily exactly K
})


test_that("initializeSplitZ_PriorClustering merges when prior has more clusters", {
    priorZ <- sample(1:10, ncol(counts), replace = TRUE)
    K <- 3

    z <- celda:::.initializeSplitZ_PriorClustering(counts, K, priorZ)

    expect_equal(length(z), ncol(counts))
    expect_equal(length(unique(z)), K)
})


test_that("initializeSplitZ_PriorClustering validates input", {
    priorZ <- sample(1:3, ncol(counts) - 10, replace = TRUE)  # Wrong length
    K <- 3

    expect_error(
        celda:::.initializeSplitZ_PriorClustering(counts, K, priorZ),
        "must have length equal"
    )
})


# Test integration with celda_C
test_that("celda_C works with marker-guided initialization", {
    # Create SCE object
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts)
    )
    sce <- selectFeatures(sce, minCount = 1)

    # Run celda_C with marker genes
    result <- celda_C(sce,
        K = 3,
        zInitialize = "split",
        markerGenes = markerList,
        maxIter = 5,
        nchains = 1,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE)

    expect_s4_class(result, "SingleCellExperiment")

    # Extract cluster assignments
    z <- celdaClusters(result)$z
    expect_equal(length(z), ncol(counts))
    expect_equal(length(unique(z)), 3)
})


test_that("celda_C works with prior clustering initialization", {
    # Create SCE object
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts)
    )
    sce <- selectFeatures(sce, minCount = 1)

    # Create prior clustering
    priorZ <- sample(1:3, ncol(counts), replace = TRUE)

    # Run celda_C with prior clustering
    result <- celda_C(sce,
        K = 3,
        zInitialize = "split",
        priorClustering = priorZ,
        maxIter = 5,
        nchains = 1,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE)

    expect_s4_class(result, "SingleCellExperiment")

    # Extract cluster assignments
    z <- celdaClusters(result)$z
    expect_equal(length(z), ncol(counts))
    expect_equal(length(unique(z)), 3)
})


# Test backward compatibility
test_that("celda_C backward compatibility without marker genes", {
    # Create SCE object
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts)
    )
    sce <- selectFeatures(sce, minCount = 1)

    # Run celda_C without marker genes (should use existing split method)
    result <- celda_C(sce,
        K = 3,
        zInitialize = "split",
        maxIter = 5,
        nchains = 1,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE)

    expect_s4_class(result, "SingleCellExperiment")

    # Extract cluster assignments
    z <- celdaClusters(result)$z
    expect_equal(length(z), ncol(counts))
    expect_equal(length(unique(z)), 3)
})


test_that("celda_C with random initialization still works", {
    # Create SCE object
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts)
    )
    sce <- selectFeatures(sce, minCount = 1)

    # Run celda_C with random initialization
    result <- celda_C(sce,
        K = 3,
        zInitialize = "random",
        maxIter = 5,
        nchains = 1,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE)

    expect_s4_class(result, "SingleCellExperiment")

    # Extract cluster assignments
    z <- celdaClusters(result)$z
    expect_equal(length(z), ncol(counts))
    expect_equal(length(unique(z)), 3)
})
