# Internal Cluster Validation Metrics for Model Selection
#
# This file contains internal functions for calculating cluster validation
# metrics during model fitting. These metrics are used in combination with
# log-likelihood for better model selection.


#' @title Calculate Silhouette Scores with Sampling
#' @description Calculates silhouette scores for cell clustering with optional
#'   sampling for computational efficiency. Uses correlation-based distance.
#' @param counts Sparse or dense matrix (features x cells)
#' @param z Integer vector of cell cluster assignments
#' @param maxCells Integer. Maximum number of cells to use. If ncol(counts) >
#'   maxCells, will sample maxCells cells. Default 5000.
#' @return Numeric. Mean silhouette score across all cells
#' @keywords internal
.calculateSilhouette_Sampled <- function(counts, z, maxCells = 5000) {
    nCells <- ncol(counts)

    # Handle edge cases
    if (nCells < 3) {
        return(NA_real_)
    }

    # Number of clusters
    K <- length(unique(z))
    if (K < 2) {
        return(NA_real_)
    }

    # Sample cells if dataset is large
    if (nCells > maxCells) {
        # Stratified sampling to maintain cluster proportions
        idx <- .stratifiedSample(z, maxCells)
        counts <- counts[, idx, drop = FALSE]
        z <- z[idx]
        nCells <- length(idx)
    }

    # Calculate correlation-based distance
    # Use pearson correlation for efficiency with sparse matrices
    tryCatch({
        # For sparse matrices, convert to dense for correlation
        if (methods::is(counts, "sparseMatrix")) {
            # Only use subset of genes with highest variance for speed
            if (nrow(counts) > 1000) {
                geneVar <- Matrix::rowMeans((counts - Matrix::rowMeans(counts))^2)
                topGenes <- order(geneVar, decreasing = TRUE)[seq_len(1000)]
                counts <- as.matrix(counts[topGenes, ])
            } else {
                counts <- as.matrix(counts)
            }
        }

        # Calculate correlation distance: dist = 1 - abs(cor)
        # This is faster than Euclidean distance for high-dimensional data
        corMat <- stats::cor(counts, method = "pearson")
        distMat <- as.dist(1 - corMat)

        # Calculate silhouette using cluster package
        if (requireNamespace("cluster", quietly = TRUE)) {
            silResult <- cluster::silhouette(z, distMat)
            meanSil <- mean(silResult[, 3])
        } else {
            # Fallback: simple silhouette calculation
            meanSil <- .simpleSilhouette(distMat, z)
        }

        return(meanSil)
    }, error = function(e) {
        warning("Error calculating silhouette: ", e$message)
        return(NA_real_)
    })
}


#' @title Calculate Calinski-Harabasz Index
#' @description Calculates the Calinski-Harabasz (CH) index, also known as
#'   variance ratio criterion. Higher values indicate better clustering.
#'   CH = (BCSS / (K-1)) / (WCSS / (N-K))
#' @param counts Matrix (features x cells)
#' @param z Integer vector of cell cluster assignments
#' @return Numeric. CH index (higher is better)
#' @keywords internal
.calculateCH <- function(counts, z) {
    nCells <- ncol(counts)
    K <- length(unique(z))

    # Handle edge cases
    if (K < 2) {
        return(NA_real_)
    }
    if (nCells <= K) {
        return(NA_real_)
    }

    tryCatch({
        # Calculate overall centroid
        overallCentroid <- Matrix::rowMeans(counts)

        # Calculate between-cluster sum of squares (BCSS)
        bcss <- 0
        for (k in unique(z)) {
            clusterCells <- which(z == k)
            nk <- length(clusterCells)
            if (nk > 0) {
                clusterCentroid <- Matrix::rowMeans(counts[, clusterCells, drop = FALSE])
                bcss <- bcss + nk * sum((clusterCentroid - overallCentroid)^2)
            }
        }

        # Calculate within-cluster sum of squares (WCSS)
        wcss <- 0
        for (k in unique(z)) {
            clusterCells <- which(z == k)
            if (length(clusterCells) > 0) {
                clusterCentroid <- Matrix::rowMeans(counts[, clusterCells, drop = FALSE])
                clusterCounts <- counts[, clusterCells, drop = FALSE]
                wcss <- wcss + sum((clusterCounts - clusterCentroid)^2)
            }
        }

        # Calculate CH index
        ch <- (bcss / (K - 1)) / (wcss / (nCells - K))

        return(ch)
    }, error = function(e) {
        warning("Error calculating CH index: ", e$message)
        return(NA_real_)
    })
}


#' @title Calculate Davies-Bouldin Index
#' @description Calculates the Davies-Bouldin (DB) index. Lower values indicate
#'   better clustering (better separation and compactness).
#' @param counts Matrix (features x cells)
#' @param z Integer vector of cell cluster assignments
#' @return Numeric. DB index (lower is better)
#' @keywords internal
.calculateDB <- function(counts, z) {
    K <- length(unique(z))

    # Handle edge cases
    if (K < 2) {
        return(NA_real_)
    }

    tryCatch({
        # Calculate cluster centroids and within-cluster scatter
        centroids <- list()
        scatters <- numeric(K)
        clusterIds <- sort(unique(z))

        for (i in seq_along(clusterIds)) {
            k <- clusterIds[i]
            clusterCells <- which(z == k)
            nk <- length(clusterCells)

            if (nk > 0) {
                # Cluster centroid
                centroids[[i]] <- Matrix::rowMeans(counts[, clusterCells, drop = FALSE])

                # Within-cluster scatter (average distance to centroid)
                if (nk > 1) {
                    clusterCounts <- counts[, clusterCells, drop = FALSE]
                    distances <- sqrt(colSums((clusterCounts - centroids[[i]])^2))
                    scatters[i] <- mean(distances)
                } else {
                    scatters[i] <- 0
                }
            } else {
                centroids[[i]] <- rep(0, nrow(counts))
                scatters[i] <- 0
            }
        }

        # Calculate DB index
        dbValues <- numeric(K)
        for (i in seq_len(K)) {
            maxRatio <- 0
            for (j in seq_len(K)) {
                if (i != j) {
                    # Distance between centroids
                    centroidDist <- sqrt(sum((centroids[[i]] - centroids[[j]])^2))

                    # Avoid division by zero
                    if (centroidDist > 0) {
                        ratio <- (scatters[i] + scatters[j]) / centroidDist
                        maxRatio <- max(maxRatio, ratio)
                    }
                }
            }
            dbValues[i] <- maxRatio
        }

        db <- mean(dbValues)
        return(db)
    }, error = function(e) {
        warning("Error calculating DB index: ", e$message)
        return(NA_real_)
    })
}


#' @title Calculate Module Coherence
#' @description Calculates gene module coherence by measuring within-module
#'   gene correlations. Higher values indicate more coherent modules.
#' @param counts Matrix (features x cells)
#' @param y Integer vector of gene module assignments
#' @param maxCells Integer. Maximum number of cells to use for correlation
#'   calculation. Default 1000.
#' @return Numeric. Mean within-module correlation (higher is better)
#' @keywords internal
.calculateModuleCoherence <- function(counts, y, maxCells = 1000) {
    nCells <- ncol(counts)
    L <- length(unique(y))

    # Handle edge cases
    if (L < 2) {
        return(NA_real_)
    }

    # Sample cells if too many for efficient correlation
    if (nCells > maxCells) {
        cellIdx <- sample(nCells, maxCells)
        counts <- counts[, cellIdx, drop = FALSE]
    }

    tryCatch({
        # Calculate mean correlation within each module
        moduleCorrelations <- numeric(L)

        for (l in unique(y)) {
            moduleGenes <- which(y == l)
            nGenes <- length(moduleGenes)

            if (nGenes >= 2) {
                # Get module gene counts
                moduleCounts <- as.matrix(counts[moduleGenes, , drop = FALSE])

                # Sample genes if module is very large
                if (nGenes > 100) {
                    geneIdx <- sample(nGenes, 100)
                    moduleCounts <- moduleCounts[geneIdx, , drop = FALSE]
                    nGenes <- 100
                }

                # Calculate pairwise correlations
                corMat <- stats::cor(t(moduleCounts), method = "pearson")

                # Mean of upper triangle (excluding diagonal)
                if (nGenes > 1) {
                    moduleCorrelations[l] <- mean(corMat[upper.tri(corMat)])
                } else {
                    moduleCorrelations[l] <- NA
                }
            } else {
                moduleCorrelations[l] <- NA
            }
        }

        # Return mean coherence across all modules
        meanCoherence <- mean(moduleCorrelations, na.rm = TRUE)

        # Handle case where all modules have < 2 genes
        if (is.nan(meanCoherence)) {
            meanCoherence <- NA_real_
        }

        return(meanCoherence)
    }, error = function(e) {
        warning("Error calculating module coherence: ", e$message)
        return(NA_real_)
    })
}


#' @title Calculate All Cluster Metrics
#' @description Main dispatcher function that calculates all requested cluster
#'   validation metrics.
#' @param counts Matrix (features x cells)
#' @param z Integer vector of cell cluster assignments
#' @param y Integer vector of gene module assignments (optional, for celda_CG)
#' @param metrics Character vector. Which metrics to calculate. Options:
#'   "silhouette", "calinskiHarabasz", "daviesBouldin", "clusterSizeCV",
#'   "moduleCoherence", "all". Default "all".
#' @return Named list with calculated metrics
#' @keywords internal
.calculateClusterMetrics <- function(counts,
                                     z,
                                     y = NULL,
                                     metrics = "all") {

    # Determine which metrics to calculate
    if ("all" %in% metrics) {
        calcMetrics <- c("silhouette", "calinskiHarabasz", "daviesBouldin",
                        "clusterSizeCV")
        if (!is.null(y)) {
            calcMetrics <- c(calcMetrics, "moduleCoherence")
        }
    } else {
        calcMetrics <- metrics
    }

    results <- list()

    # Silhouette score
    if ("silhouette" %in% calcMetrics) {
        results$silhouette <- tryCatch({
            .calculateSilhouette_Sampled(counts, z, maxCells = 5000)
        }, error = function(e) {
            warning("Skipping silhouette calculation: ", e$message)
            NA_real_
        })
    }

    # Calinski-Harabasz index
    if ("calinskiHarabasz" %in% calcMetrics) {
        results$calinskiHarabasz <- tryCatch({
            .calculateCH(counts, z)
        }, error = function(e) {
            warning("Skipping CH calculation: ", e$message)
            NA_real_
        })
    }

    # Davies-Bouldin index
    if ("daviesBouldin" %in% calcMetrics) {
        results$daviesBouldin <- tryCatch({
            .calculateDB(counts, z)
        }, error = function(e) {
            warning("Skipping DB calculation: ", e$message)
            NA_real_
        })
    }

    # Cluster size coefficient of variation
    if ("clusterSizeCV" %in% calcMetrics) {
        clusterSizes <- table(z)
        results$clusterSizeCV <- stats::sd(clusterSizes) / mean(clusterSizes)
    }

    # Module coherence (only if y provided)
    if ("moduleCoherence" %in% calcMetrics && !is.null(y)) {
        results$moduleCoherence <- tryCatch({
            .calculateModuleCoherence(counts, y, maxCells = 1000)
        }, error = function(e) {
            warning("Skipping module coherence calculation: ", e$message)
            NA_real_
        })
    }

    return(results)
}


#' @title Normalize Metrics for Combined Score
#' @description Normalizes validation metrics to [0, 1] scale for combining
#'   into a single score. Handles both "higher is better" and "lower is better"
#'   metrics.
#' @param metricsList List of named lists, where each list contains metrics
#'   from one chain
#' @return List with same structure but normalized values
#' @keywords internal
.normalizeMetrics <- function(metricsList) {
    # Extract all values for each metric across chains
    allMetrics <- names(metricsList[[1]])

    normalized <- vector("list", length(metricsList))

    for (metricName in allMetrics) {
        values <- sapply(metricsList, function(x) x[[metricName]])

        # Skip if all NA
        if (all(is.na(values))) {
            for (i in seq_along(metricsList)) {
                normalized[[i]][[metricName]] <- NA_real_
            }
            next
        }

        # Remove NA for min/max calculation
        validValues <- values[!is.na(values)]
        minVal <- min(validValues)
        maxVal <- max(validValues)

        # Avoid division by zero
        if (maxVal == minVal) {
            for (i in seq_along(metricsList)) {
                normalized[[i]][[metricName]] <- 0.5
            }
        } else {
            # Normalize based on metric type
            if (metricName == "daviesBouldin" || metricName == "clusterSizeCV") {
                # Lower is better: invert the normalization
                for (i in seq_along(metricsList)) {
                    if (is.na(values[i])) {
                        normalized[[i]][[metricName]] <- NA_real_
                    } else {
                        normalized[[i]][[metricName]] <-
                            1 - (values[i] - minVal) / (maxVal - minVal)
                    }
                }
            } else {
                # Higher is better: standard normalization
                for (i in seq_along(metricsList)) {
                    if (is.na(values[i])) {
                        normalized[[i]][[metricName]] <- NA_real_
                    } else {
                        normalized[[i]][[metricName]] <-
                            (values[i] - minVal) / (maxVal - minVal)
                    }
                }
            }
        }
    }

    return(normalized)
}


#' @title Calculate Combined Model Score
#' @description Combines log-likelihood and validation metrics into a single
#'   score for model selection.
#' @param logLik Numeric. Log-likelihood value
#' @param metrics Named list of validation metrics
#' @param normalizedMetrics Named list of normalized (0-1) validation metrics
#' @param validationWeight Numeric. Weight for validation metrics (0-1).
#'   Weight for log-likelihood is (1 - validationWeight). Default 0.3.
#' @return Numeric. Combined score (higher is better)
#' @keywords internal
.calculateCombinedScore <- function(logLik,
                                    metrics,
                                    normalizedMetrics,
                                    validationWeight = 0.3) {

    # Validation weight should be between 0 and 1
    validationWeight <- max(0, min(1, validationWeight))
    llWeight <- 1 - validationWeight

    # Calculate weighted validation score
    # Weights for different metrics (should sum to 1)
    metricWeights <- list(
        silhouette = 0.15,          # Silhouette is expensive, lower weight
        calinskiHarabasz = 0.35,    # CH is reliable and cheap
        daviesBouldin = 0.25,       # DB is reliable and cheap
        clusterSizeCV = 0.05,       # Penalize unbalanced clusters
        moduleCoherence = 0.20      # Module quality (if available)
    )

    # Calculate validation score from available metrics
    validationScore <- 0
    totalWeight <- 0

    for (metricName in names(normalizedMetrics)) {
        if (!is.na(normalizedMetrics[[metricName]]) &&
            metricName %in% names(metricWeights)) {
            validationScore <- validationScore +
                normalizedMetrics[[metricName]] * metricWeights[[metricName]]
            totalWeight <- totalWeight + metricWeights[[metricName]]
        }
    }

    # Normalize validation score by total weight of available metrics
    if (totalWeight > 0) {
        validationScore <- validationScore / totalWeight
    } else {
        # No valid metrics, fall back to log-likelihood only
        validationWeight <- 0
        llWeight <- 1
    }

    # Combined score (note: logLik will be normalized separately before calling)
    combinedScore <- llWeight * logLik + validationWeight * validationScore

    return(combinedScore)
}


#' @title Stratified Sampling
#' @description Performs stratified sampling to maintain cluster proportions
#' @param z Integer vector of cluster assignments
#' @param n Integer. Number of samples to draw
#' @return Integer vector of indices
#' @keywords internal
.stratifiedSample <- function(z, n) {
    K <- length(unique(z))
    nTotal <- length(z)

    if (n >= nTotal) {
        return(seq_len(nTotal))
    }

    # Calculate samples per cluster proportional to cluster size
    clusterSizes <- table(z)
    samplesPerCluster <- round(clusterSizes / nTotal * n)

    # Adjust to exactly n samples
    diff <- sum(samplesPerCluster) - n
    if (diff > 0) {
        # Remove from largest clusters
        largest <- order(samplesPerCluster, decreasing = TRUE)[seq_len(diff)]
        for (i in largest) {
            samplesPerCluster[i] <- samplesPerCluster[i] - 1
        }
    } else if (diff < 0) {
        # Add to largest clusters
        largest <- order(samplesPerCluster, decreasing = TRUE)[seq_len(-diff)]
        for (i in largest) {
            samplesPerCluster[i] <- samplesPerCluster[i] + 1
        }
    }

    # Sample from each cluster
    idx <- integer(0)
    for (k in names(clusterSizes)) {
        clusterIdx <- which(z == as.integer(k))
        nSamples <- samplesPerCluster[k]
        if (nSamples > 0) {
            if (nSamples >= length(clusterIdx)) {
                idx <- c(idx, clusterIdx)
            } else {
                idx <- c(idx, sample(clusterIdx, nSamples))
            }
        }
    }

    return(idx)
}


#' @title Simple Silhouette Calculation
#' @description Fallback silhouette calculation if cluster package not available
#' @param distMat Distance matrix
#' @param z Cluster assignments
#' @return Mean silhouette score
#' @keywords internal
.simpleSilhouette <- function(distMat, z) {
    distMat <- as.matrix(distMat)
    n <- length(z)
    K <- length(unique(z))

    if (K < 2) {
        return(NA_real_)
    }

    silScores <- numeric(n)

    for (i in seq_len(n)) {
        # Average distance to points in same cluster (a)
        sameCluster <- which(z == z[i] & seq_len(n) != i)
        if (length(sameCluster) > 0) {
            a <- mean(distMat[i, sameCluster])
        } else {
            a <- 0
        }

        # Minimum average distance to points in other clusters (b)
        b <- Inf
        for (k in setdiff(unique(z), z[i])) {
            otherCluster <- which(z == k)
            if (length(otherCluster) > 0) {
                avgDist <- mean(distMat[i, otherCluster])
                b <- min(b, avgDist)
            }
        }

        # Silhouette score
        if (max(a, b) > 0) {
            silScores[i] <- (b - a) / max(a, b)
        } else {
            silScores[i] <- 0
        }
    }

    return(mean(silScores))
}
