# Test consensus clustering functionality
library(celda)
library(testthat)

context("Consensus Clustering")

# Create test data
set.seed(12345)
testCounts <- simulateCells(model = "celda_CG", K = 3, L = 5, S = 1, G = 50, C = 50)

test_that(".buildCooccurrenceMatrix works correctly", {
    # Create mock chain results with perfect agreement
    chain1 <- list(z = rep(1:3, each = 10), y = rep(1:5, each = 10))
    chain2 <- list(z = rep(1:3, each = 10), y = rep(1:5, each = 10))
    chainResults <- list(chain1, chain2)

    # Build co-occurrence matrix for cells
    cooccur_z <- celda:::.buildCooccurrenceMatrix(chainResults, type = "z")

    # Check dimensions
    expect_equal(nrow(cooccur_z), 30)
    expect_equal(ncol(cooccur_z), 30)

    # Check diagonal is 1
    expect_true(all(diag(cooccur_z) == 1))

    # Check symmetry
    expect_equal(cooccur_z, t(cooccur_z))

    # Check cells in same cluster have cooccurrence = 1
    expect_equal(cooccur_z[1, 2], 1)  # Both in cluster 1
    expect_equal(cooccur_z[11, 12], 1)  # Both in cluster 2

    # Check cells in different clusters have cooccurrence = 0
    expect_equal(cooccur_z[1, 11], 0)  # Different clusters
})

test_that(".buildCooccurrenceMatrix handles partial agreement", {
    # Create chains with partial agreement
    chain1 <- list(z = c(rep(1, 10), rep(2, 10), rep(3, 10)))
    chain2 <- list(z = c(rep(1, 10), rep(3, 10), rep(2, 10)))  # Swapped clusters 2 and 3
    chainResults <- list(chain1, chain2)

    cooccur <- celda:::.buildCooccurrenceMatrix(chainResults, type = "z")

    # First 10 cells should have perfect agreement
    expect_equal(cooccur[1, 2], 1)

    # Cells 11 and 21 disagree (different clusters in both chains)
    expect_equal(cooccur[11, 21], 0)

    # Cells 11 and 12 agree in chain1 but not chain2
    expect_equal(cooccur[11, 12], 0.5)
})

test_that(".consensusClustering with cooccurrence method works", {
    # Perfect agreement case
    chain1 <- list(z = rep(1:3, each = 10), y = rep(1:2, each = 25))
    chain2 <- list(z = rep(1:3, each = 10), y = rep(1:2, each = 25))
    chainResults <- list(chain1, chain2)

    consensus <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "cooccurrence",
        type = "z"
    )

    # Check output structure
    expect_true(is.list(consensus))
    expect_true(all(c("assignments", "confidence", "lowConfidenceIndices") %in% names(consensus)))

    # Check assignments
    expect_equal(length(consensus$assignments), 30)
    expect_equal(length(unique(consensus$assignments)), 3)  # 3 clusters

    # Check confidence
    expect_equal(length(consensus$confidence), 30)
    expect_true(all(consensus$confidence >= 0 & consensus$confidence <= 1))

    # With perfect agreement, confidence should be high
    expect_true(all(consensus$confidence > 0.9))

    # No low-confidence cells with perfect agreement
    expect_equal(length(consensus$lowConfidenceIndices), 0)
})

test_that(".consensusClustering with median method works", {
    # Create chains with some disagreement
    chain1 <- list(z = c(rep(1, 10), rep(2, 10), rep(3, 10)))
    chain2 <- list(z = c(rep(1, 10), rep(2, 10), rep(3, 10)))
    chain3 <- list(z = c(rep(1, 10), rep(3, 10), rep(2, 10)))  # Different for last 20
    chainResults <- list(chain1, chain2, chain3)

    consensus <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "median",
        type = "z",
        minAgreement = 0.6
    )

    # First 10 cells should have full agreement (confidence = 1)
    expect_equal(consensus$confidence[1:10], rep(1, 10))

    # Last 20 cells have 2/3 agreement (confidence = 0.67)
    expect_true(all(consensus$confidence[11:30] < 1))
    expect_true(all(consensus$confidence[11:30] >= 0.6))

    # Some cells might be low confidence depending on threshold
    if (length(consensus$lowConfidenceIndices) > 0) {
        expect_true(all(consensus$lowConfidenceIndices > 10))
    }
})

test_that(".consensusClustering_CG works for both z and y", {
    # Create mock chains
    chain1 <- list(
        z = rep(1:3, each = 10),
        y = rep(1:5, each = 10),
        finalLogLik = 1000
    )
    chain2 <- list(
        z = rep(1:3, each = 10),
        y = rep(1:5, each = 10),
        finalLogLik = 1000
    )
    chainResults <- list(chain1, chain2)

    consensus <- celda:::.consensusClustering_CG(
        allChainResults = chainResults,
        method = "cooccurrence",
        minAgreement = 0.7
    )

    # Check output structure
    expect_true(all(c("z", "y", "zConfidence", "yConfidence",
                      "lowConfidenceCells", "lowConfidenceGenes") %in% names(consensus)))

    # Check z (cells)
    expect_equal(length(consensus$z), 30)
    expect_equal(length(consensus$zConfidence), 30)

    # Check y (genes)
    expect_equal(length(consensus$y), 50)
    expect_equal(length(consensus$yConfidence), 50)
})

test_that("Consensus is robust to label switching", {
    # Chains with permuted labels (label switching problem)
    chain1 <- list(z = c(rep(1, 10), rep(2, 10), rep(3, 10)))
    chain2 <- list(z = c(rep(2, 10), rep(3, 10), rep(1, 10)))  # Permuted labels
    chain3 <- list(z = c(rep(3, 10), rep(1, 10), rep(2, 10)))  # Different permutation
    chainResults <- list(chain1, chain2, chain3)

    consensus <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "cooccurrence",
        type = "z"
    )

    # Consensus should still correctly identify the 3 groups
    # Elements 1-10, 11-20, 21-30 should cluster together
    assign1 <- consensus$assignments[1:10]
    assign2 <- consensus$assignments[11:20]
    assign3 <- consensus$assignments[21:30]

    # Each group should have only one unique assignment
    expect_equal(length(unique(assign1)), 1)
    expect_equal(length(unique(assign2)), 1)
    expect_equal(length(unique(assign3)), 1)

    # And these should be different from each other
    expect_true(unique(assign1) != unique(assign2))
    expect_true(unique(assign2) != unique(assign3))
    expect_true(unique(assign1) != unique(assign3))
})

test_that("Confidence scores reflect chain agreement", {
    # High agreement for some cells, low for others
    chain1 <- list(z = c(rep(1, 10), rep(2, 5), rep(3, 5)))
    chain2 <- list(z = c(rep(1, 10), rep(2, 5), rep(3, 5)))
    chain3 <- list(z = c(rep(1, 10), rep(3, 5), rep(2, 5)))  # Disagree on last 10
    chainResults <- list(chain1, chain2, chain3)

    consensus <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "median",
        type = "z"
    )

    # First 10 cells have perfect agreement - high confidence
    expect_true(all(consensus$confidence[1:10] == 1))

    # Last 10 cells have disagreement - lower confidence
    expect_true(all(consensus$confidence[11:20] < 1))
})

test_that("Edge cases are handled correctly", {
    # Single chain - should just return that chain's assignments
    chain1 <- list(z = rep(1:3, each = 10))
    chainResults <- list(chain1)

    consensus <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "median",
        type = "z"
    )

    # Should have same assignments as input
    expect_equal(length(unique(consensus$assignments)), 3)
    expect_equal(length(consensus$assignments), 30)

    # Confidence should be perfect (only one chain to agree with)
    expect_true(all(consensus$confidence == 1))
})

test_that("Empty chain list throws error", {
    expect_error(
        celda:::.consensusClustering(
            allChainResults = list(),
            method = "median",
            type = "z"
        ),
        "must contain at least one chain"
    )
})

test_that("Chains with different lengths throw error", {
    chain1 <- list(z = 1:10)
    chain2 <- list(z = 1:15)  # Different length
    chainResults <- list(chain1, chain2)

    expect_error(
        celda:::.buildCooccurrenceMatrix(chainResults, type = "z"),
        "same number of elements"
    )
})

test_that("Consensus with many chains is stable", {
    # Create 10 chains with mostly agreement
    set.seed(123)
    nChains <- 10
    chainResults <- lapply(seq_len(nChains), function(i) {
        # Add small random noise to cluster assignments
        baseZ <- rep(1:3, each = 10)
        if (i > 1) {
            # Randomly swap a few assignments
            swapIdx <- sample(30, 3)
            baseZ[swapIdx] <- sample(1:3, 3, replace = TRUE)
        }
        list(z = baseZ)
    })

    consensus <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "cooccurrence",
        type = "z"
    )

    # Should still recover 3 main clusters
    expect_true(length(unique(consensus$assignments)) <= 4)  # Allow slight over-clustering

    # Median method should be similar
    consensusMedian <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "median",
        type = "z"
    )

    expect_true(length(unique(consensusMedian$assignments)) <= 4)
})

test_that("Low confidence threshold works correctly", {
    chain1 <- list(z = c(rep(1, 10), rep(2, 10)))
    chain2 <- list(z = c(rep(1, 10), rep(2, 10)))
    chain3 <- list(z = c(rep(2, 10), rep(1, 10)))  # Opposite
    chainResults <- list(chain1, chain2, chain3)

    # With high threshold, more cells should be flagged
    consensus_high <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "median",
        type = "z",
        minAgreement = 0.9
    )

    # With low threshold, fewer cells should be flagged
    consensus_low <- celda:::.consensusClustering(
        allChainResults = chainResults,
        method = "median",
        type = "z",
        minAgreement = 0.5
    )

    # High threshold should flag more cells
    expect_true(
        length(consensus_high$lowConfidenceIndices) >=
        length(consensus_low$lowConfidenceIndices)
    )
})
