#' @title Calculate Gene Weights for Adaptive Feature Weighting
#' @description Calculates importance scores for genes based on their
#'  informativeness for clustering. Genes with high between-cluster variance
#'  relative to total variance receive higher weights.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param z Integer vector. Cell cluster assignments.
#' @param y Integer vector. Gene module assignments.
#' @param method Character. Method for calculating gene weights. Currently only
#'  "variance" is supported. Default "variance".
#' @return Numeric vector of gene weights (one per gene). Weights are scaled
#'  such that mean weight equals 1, and are clipped to range [0.1, 10] to
#'  avoid extreme values.
#' @keywords internal
.calculateGeneWeights <- function(counts,
                                   z,
                                   y,
                                   method = "variance") {
    nG <- nrow(counts)
    K <- max(z)

    # Initialize weights to 1
    weights <- rep(1, nG)

    if (method == "variance") {
        # Calculate variance-based weights for each gene
        for (i in seq_len(nG)) {
            geneExpression <- as.numeric(counts[i, ])

            # Skip genes with zero or constant expression
            if (stats::var(geneExpression) == 0) {
                weights[i] <- 0.1  # Minimum weight for uninformative genes
                next
            }

            # Calculate cluster means
            clusterMeans <- vapply(seq_len(K), function(k) {
                cellsInCluster <- which(z == k)
                if (length(cellsInCluster) > 0) {
                    mean(geneExpression[cellsInCluster])
                } else {
                    0
                }
            }, numeric(1))

            # Calculate between-cluster variance
            betweenVar <- stats::var(clusterMeans)

            # Calculate total variance
            totalVar <- stats::var(geneExpression)

            # Calculate importance score (0 to 1)
            # Higher when gene expression varies between clusters
            if (totalVar > 0) {
                score <- betweenVar / totalVar
                weights[i] <- score
            } else {
                weights[i] <- 0.1
            }
        }

        # Apply softmax-like transformation to convert scores to weights
        # This emphasizes differences while keeping weights positive
        # Use log-sum-exp trick for numerical stability
        maxWeight <- max(weights)
        expWeights <- exp(weights - maxWeight)
        sumExpWeights <- sum(expWeights)
        weights <- (expWeights / sumExpWeights) * nG

        # Clip weights to avoid extreme values
        # Minimum of 0.1, maximum of 10
        weights <- pmax(0.1, pmin(10, weights))
    }

    return(weights)
}


#' @title Calculate Gene Weights with Marker Gene Boosting
#' @description Calculates gene importance scores with optional boosting for
#'  known marker genes. First calculates variance-based weights, then
#'  multiplies weights of specified marker genes by a boost factor.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param z Integer vector. Cell cluster assignments.
#' @param y Integer vector. Gene module assignments.
#' @param markerGenes Integer vector or NULL. Indices of genes to boost
#'  (1-indexed). If NULL, no boosting is performed. Default NULL.
#' @param boostFactor Numeric. Multiplicative factor to boost marker gene
#'  weights. Default 2.
#' @return Numeric vector of gene weights (one per gene). Weights are scaled
#'  such that mean weight equals 1 after boosting.
#' @keywords internal
.calculateGeneWeights_MarkerBoosted <- function(counts,
                                                  z,
                                                  y,
                                                  markerGenes = NULL,
                                                  boostFactor = 2) {
    # Calculate base variance weights
    weights <- .calculateGeneWeights(counts, z, y, method = "variance")

    # Apply marker gene boosting if specified
    if (!is.null(markerGenes) && length(markerGenes) > 0) {
        # Validate marker gene indices
        nG <- nrow(counts)
        validMarkers <- markerGenes[markerGenes >= 1 & markerGenes <= nG]

        if (length(validMarkers) > 0) {
            # Boost marker gene weights
            weights[validMarkers] <- weights[validMarkers] * boostFactor

            # Renormalize weights to have mean = 1
            weights <- weights / mean(weights) * nG / nG
        }
    }

    return(weights)
}
