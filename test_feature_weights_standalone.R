#!/usr/bin/env Rscript
# Standalone test script for feature weighting
# This can be run independently to test the basic logic

# Simulate the .calculateGeneWeights function logic
calculateGeneWeights_test <- function(counts, z, y, method = "variance") {
    nG <- nrow(counts)
    K <- max(z)

    weights <- rep(1, nG)

    if (method == "variance") {
        for (i in seq_len(nG)) {
            geneExpression <- as.numeric(counts[i, ])

            if (stats::var(geneExpression) == 0) {
                weights[i] <- 0.1
                next
            }

            clusterMeans <- vapply(seq_len(K), function(k) {
                cellsInCluster <- which(z == k)
                if (length(cellsInCluster) > 0) {
                    mean(geneExpression[cellsInCluster])
                } else {
                    0
                }
            }, numeric(1))

            betweenVar <- stats::var(clusterMeans)
            totalVar <- stats::var(geneExpression)

            if (totalVar > 0) {
                score <- betweenVar / totalVar
                weights[i] <- score
            } else {
                weights[i] <- 0.1
            }
        }

        # Apply softmax-like transformation
        maxWeight <- max(weights)
        expWeights <- exp(weights - maxWeight)
        sumExpWeights <- sum(expWeights)
        weights <- (expWeights / sumExpWeights) * nG

        # Clip weights
        weights <- pmax(0.1, pmin(10, weights))
    }

    return(weights)
}

# Test with synthetic data
set.seed(12345)
nGenes <- 100
nCells <- 50
K <- 3

# Create simple clustering
z <- rep(1:K, length.out = nCells)
counts <- matrix(0, nrow = nGenes, ncol = nCells)

# First 20 genes are informative
for (i in 1:20) {
    for (k in 1:K) {
        cells_in_cluster <- which(z == k)
        counts[i, cells_in_cluster] <- rpois(length(cells_in_cluster),
            lambda = k * 10)
    }
}

# Remaining genes are noisy
for (i in 21:nGenes) {
    counts[i, ] <- rpois(nCells, lambda = 2)
}

y <- rep(1:5, length.out = nGenes)

# Calculate weights
weights <- calculateGeneWeights_test(counts, z, y)

# Check results
cat("Weight calculation test:\n")
cat("- Number of genes:", length(weights), "\n")
cat("- All weights finite:", all(is.finite(weights)), "\n")
cat("- All weights in range [0.1, 10]:",
    all(weights >= 0.1 & weights <= 10), "\n")
cat("- Mean weight (informative genes):",
    mean(weights[1:20]), "\n")
cat("- Mean weight (noisy genes):",
    mean(weights[21:nGenes]), "\n")
cat("- Informative > Noisy:",
    mean(weights[1:20]) > mean(weights[21:nGenes]), "\n")

cat("\nTest PASSED if all checks are TRUE\n")
