#' @title Calculate marker scores for cells
#' @description For each cell, calculates expression scores for each marker
#'  gene set. Scores are based on the mean or median expression of marker genes
#'  in each cell, normalized across all marker sets.
#' @param counts Integer matrix of counts (genes x cells)
#' @param markerGenes Named list where names are cell types and values are
#'  character vectors of gene names. Example:
#'  list("T cells" = c("CD3D", "CD3E"), "B cells" = c("CD19", "MS4A1"))
#' @param method Character. Method for calculating marker scores. One of "mean"
#'  or "median". Default "mean".
#' @return Matrix of marker scores (cells x marker sets). Scores are normalized
#'  so each cell's scores sum to 1.
#' @keywords internal
.calculateMarkerScores <- function(counts,
                                   markerGenes,
                                   method = c("mean", "median")) {

    method <- match.arg(method)

    # Validate input
    if (!is.list(markerGenes) || is.null(names(markerGenes))) {
        stop("'markerGenes' must be a named list")
    }

    if (length(markerGenes) == 0) {
        stop("'markerGenes' must contain at least one marker set")
    }

    # Get row names (gene names)
    geneNames <- rownames(counts)
    if (is.null(geneNames)) {
        stop("'counts' must have row names (gene names)")
    }

    nCells <- ncol(counts)
    nMarkerSets <- length(markerGenes)
    markerSetNames <- names(markerGenes)

    # Initialize score matrix
    markerScores <- matrix(0,
        nrow = nCells,
        ncol = nMarkerSets,
        dimnames = list(colnames(counts), markerSetNames))

    # Calculate scores for each marker set
    for (i in seq_along(markerGenes)) {
        markers <- markerGenes[[i]]

        # Find markers that exist in the count matrix
        markersFound <- intersect(markers, geneNames)

        if (length(markersFound) == 0) {
            warning("No markers from set '", markerSetNames[i],
                "' found in count matrix. Setting scores to 0.")
            next
        }

        if (length(markersFound) < length(markers)) {
            missingMarkers <- setdiff(markers, geneNames)
            warning(length(missingMarkers), " marker(s) from set '",
                markerSetNames[i], "' not found in count matrix: ",
                paste(head(missingMarkers, 3), collapse = ", "),
                if (length(missingMarkers) > 3) "..." else "")
        }

        # Extract marker gene counts
        if (length(markersFound) == 1) {
            markerCounts <- counts[markersFound, , drop = FALSE]
        } else {
            markerCounts <- counts[markersFound, , drop = FALSE]
        }

        # Calculate score based on method
        if (method == "mean") {
            if (length(markersFound) == 1) {
                markerScores[, i] <- as.numeric(markerCounts)
            } else {
                markerScores[, i] <- colMeans(markerCounts)
            }
        } else {  # median
            if (length(markersFound) == 1) {
                markerScores[, i] <- as.numeric(markerCounts)
            } else {
                markerScores[, i] <- apply(markerCounts, 2, stats::median)
            }
        }
    }

    # Normalize scores so each cell sums to 1 (or equals 0 if all scores are 0)
    rowSums <- rowSums(markerScores)
    for (i in seq_len(nCells)) {
        if (rowSums[i] > 0) {
            markerScores[i, ] <- markerScores[i, ] / rowSums[i]
        }
    }

    return(markerScores)
}


#' @title Initialize cell clustering using marker genes
#' @description Initializes cell cluster assignments using expression of marker
#'  genes. First calculates marker scores for each cell, then uses hierarchical
#'  clustering to identify initial clusters.
#' @param counts Integer matrix of counts (genes x cells)
#' @param K Integer. Number of cell clusters to initialize
#' @param markerGenes Named list where names are cell types and values are
#'  character vectors of gene names
#' @param minCell Integer. Minimum number of cells allowed in a cluster. Default 3.
#' @param method Character. Method for calculating marker scores. One of "mean"
#'  or "median". Default "mean".
#' @return Integer vector of length ncol(counts) with initial cluster assignments
#' @keywords internal
.initializeSplitZ_MarkerGuided <- function(counts,
                                           K,
                                           markerGenes,
                                           minCell = 3,
                                           method = c("mean", "median")) {

    method <- match.arg(method)
    nCells <- ncol(counts)

    # Validate inputs
    if (K < 2) {
        stop("'K' must be at least 2")
    }

    if (nCells < K) {
        stop("Number of cells (", nCells, ") must be >= K (", K, ")")
    }

    # Calculate marker scores
    markerScores <- .calculateMarkerScores(counts, markerGenes, method = method)
    nMarkerSets <- ncol(markerScores)

    # Check if any cells have non-zero scores
    cellsWithScores <- rowSums(markerScores) > 0
    if (sum(cellsWithScores) == 0) {
        warning("No marker genes expressed in any cells. ",
            "Falling back to random initialization.")
        return(.initializeCluster(K, nCells, initial = NULL, fixed = NULL))
    }

    # Initialize z vector
    z <- integer(nCells)

    # If K <= number of marker sets, assign cells to marker sets
    if (K <= nMarkerSets) {
        # For each cell, assign to marker set with highest score
        z <- apply(markerScores, 1, which.max)

        # If we need fewer clusters than marker sets, merge similar ones
        if (K < nMarkerSets) {
            # Use hierarchical clustering on marker set similarity
            # Calculate pairwise similarity based on which cells are assigned
            markerSetSimilarity <- matrix(0, nrow = nMarkerSets, ncol = nMarkerSets)
            for (i in seq_len(nMarkerSets)) {
                for (j in seq_len(nMarkerSets)) {
                    if (i != j) {
                        # Similarity = correlation of marker scores across cells
                        markerSetSimilarity[i, j] <- stats::cor(
                            markerScores[, i],
                            markerScores[, j],
                            method = "pearson")
                    }
                }
            }

            # Convert similarity to distance
            markerSetDist <- stats::as.dist(1 - markerSetSimilarity)
            hc <- stats::hclust(markerSetDist, method = "average")
            markerSetClusters <- stats::cutree(hc, k = K)

            # Remap z to merged clusters
            z <- markerSetClusters[z]
        }
    } else {
        # K > number of marker sets, need to split clusters

        # Start by assigning cells to marker sets
        zInitial <- apply(markerScores, 1, which.max)
        z <- zInitial
        currentK <- nMarkerSets

        # Iteratively split largest/most heterogeneous clusters
        iterCount <- 0
        maxIter <- 50

        while (currentK < K && iterCount < maxIter) {
            # Count cells in each cluster
            clusterCounts <- tabulate(z, nbins = currentK)

            # Find clusters eligible for splitting (> minCell * 2)
            eligibleClusters <- which(clusterCounts >= minCell * 2)

            if (length(eligibleClusters) == 0) {
                warning("Cannot split to K=", K,
                    " clusters while maintaining minCell=", minCell,
                    ". Returning ", currentK, " clusters.")
                break
            }

            # Calculate heterogeneity for each eligible cluster
            # Heterogeneity = variance of marker scores within cluster
            heterogeneity <- vapply(eligibleClusters, function(cl) {
                cellsInCluster <- which(z == cl)
                if (length(cellsInCluster) <= 1) {
                    return(0)
                }
                scoresInCluster <- markerScores[cellsInCluster, , drop = FALSE]
                mean(apply(scoresInCluster, 2, stats::var))
            }, numeric(1))

            # Split the most heterogeneous cluster
            clusterToSplit <- eligibleClusters[which.max(heterogeneity)]
            cellsInCluster <- which(z == clusterToSplit)

            # Use k-means on marker scores to split into 2
            scoresInCluster <- markerScores[cellsInCluster, , drop = FALSE]

            # Handle case where all scores are identical
            if (all(scoresInCluster == scoresInCluster[1, 1])) {
                # Random split
                splitAssignment <- sample(rep(1:2,
                    length.out = length(cellsInCluster)))
            } else {
                kmResult <- stats::kmeans(scoresInCluster,
                    centers = 2,
                    nstart = 10)
                splitAssignment <- kmResult$cluster
            }

            # Update z: keep cluster assignment for group 1,
            # assign new cluster ID for group 2
            z[cellsInCluster[splitAssignment == 2]] <- currentK + 1
            currentK <- currentK + 1
            iterCount <- iterCount + 1
        }
    }

    # Ensure all cluster IDs from 1 to K are represented
    # and renumber to be consecutive
    z <- as.integer(as.factor(z))

    # If we have fewer clusters than K due to minCell constraints,
    # split largest clusters randomly
    currentK <- max(z)
    while (currentK < K) {
        clusterCounts <- tabulate(z, nbins = currentK)
        largestCluster <- which.max(clusterCounts)

        if (clusterCounts[largestCluster] < minCell * 2) {
            warning("Cannot achieve K=", K,
                " clusters while maintaining minCell=", minCell)
            break
        }

        cellsInCluster <- which(z == largestCluster)
        # Randomly split in half
        splitAssignment <- sample(rep(1:2, length.out = length(cellsInCluster)))
        z[cellsInCluster[splitAssignment == 2]] <- currentK + 1
        currentK <- currentK + 1
    }

    # Final renumbering to ensure consecutive IDs
    z <- as.integer(as.factor(z))

    return(z)
}


#' @title Initialize cell clustering using prior clustering results
#' @description Initializes cell cluster assignments by refining a prior
#'  clustering result. If the prior has the same number of clusters as K,
#'  it is used directly. If fewer, clusters are split. If more, clusters
#'  are merged based on expression similarity.
#' @param counts Integer matrix of counts (genes x cells)
#' @param K Integer. Target number of cell clusters
#' @param priorClustering Integer vector of length ncol(counts) with prior
#'  cluster assignments
#' @param minCell Integer. Minimum number of cells allowed in a cluster. Default 3.
#' @return Integer vector of length ncol(counts) with refined cluster assignments
#' @keywords internal
.initializeSplitZ_PriorClustering <- function(counts,
                                               K,
                                               priorClustering,
                                               minCell = 3) {

    nCells <- ncol(counts)

    # Validate inputs
    if (K < 2) {
        stop("'K' must be at least 2")
    }

    if (length(priorClustering) != nCells) {
        stop("'priorClustering' must have length equal to ncol(counts)")
    }

    # Convert to integer factor
    z <- as.integer(as.factor(priorClustering))
    currentK <- max(z)

    # Case 1: Prior clustering has exactly K clusters
    if (currentK == K) {
        return(z)
    }

    # Case 2: Prior has fewer clusters than K - need to split
    if (currentK < K) {
        iterCount <- 0
        maxIter <- 50

        while (currentK < K && iterCount < maxIter) {
            # Find clusters eligible for splitting
            clusterCounts <- tabulate(z, nbins = currentK)
            eligibleClusters <- which(clusterCounts >= minCell * 2)

            if (length(eligibleClusters) == 0) {
                warning("Cannot split to K=", K,
                    " clusters while maintaining minCell=", minCell,
                    ". Returning ", currentK, " clusters.")
                break
            }

            # Calculate within-cluster variance for each cluster
            # Split the cluster with highest variance
            variances <- vapply(eligibleClusters, function(cl) {
                cellsInCluster <- which(z == cl)
                if (length(cellsInCluster) <= 1) {
                    return(0)
                }
                clusterCounts <- counts[, cellsInCluster, drop = FALSE]
                # Calculate variance of total counts per cell
                stats::var(colSums(clusterCounts))
            }, numeric(1))

            clusterToSplit <- eligibleClusters[which.max(variances)]
            cellsInCluster <- which(z == clusterToSplit)

            # Use k-means on gene expression to split into 2
            clusterCounts <- counts[, cellsInCluster, drop = FALSE]

            # Use top variable genes for splitting to reduce computation
            geneVars <- apply(clusterCounts, 1, stats::var)
            topGenes <- order(geneVars, decreasing = TRUE)[seq_len(min(500, nrow(counts)))]
            clusterCountsSubset <- t(as.matrix(clusterCounts[topGenes, , drop = FALSE]))

            # Perform k-means
            kmResult <- stats::kmeans(clusterCountsSubset,
                centers = 2,
                nstart = 10)
            splitAssignment <- kmResult$cluster

            # Update z
            z[cellsInCluster[splitAssignment == 2]] <- currentK + 1
            currentK <- currentK + 1
            iterCount <- iterCount + 1
        }

        z <- as.integer(as.factor(z))
        return(z)
    }

    # Case 3: Prior has more clusters than K - need to merge
    if (currentK > K) {
        # Iteratively merge most similar clusters
        while (currentK > K) {
            # Calculate mean expression profile for each cluster
            clusterProfiles <- matrix(0, nrow = nrow(counts), ncol = currentK)

            for (cl in seq_len(currentK)) {
                cellsInCluster <- which(z == cl)
                if (length(cellsInCluster) == 1) {
                    clusterProfiles[, cl] <- counts[, cellsInCluster]
                } else {
                    clusterProfiles[, cl] <- rowMeans(counts[, cellsInCluster, drop = FALSE])
                }
            }

            # Calculate pairwise distances between cluster profiles
            # Use correlation-based distance
            clusterDist <- matrix(0, nrow = currentK, ncol = currentK)
            for (i in seq_len(currentK - 1)) {
                for (j in (i + 1):currentK) {
                    # Use 1 - correlation as distance
                    clusterDist[i, j] <- 1 - stats::cor(
                        clusterProfiles[, i],
                        clusterProfiles[, j],
                        method = "pearson")
                    clusterDist[j, i] <- clusterDist[i, j]
                }
            }

            # Find pair of clusters with minimum distance
            # Avoid merging with itself
            diag(clusterDist) <- Inf
            minDist <- which(clusterDist == min(clusterDist), arr.ind = TRUE)[1, ]
            cluster1 <- minDist[1]
            cluster2 <- minDist[2]

            # Merge cluster2 into cluster1
            z[z == cluster2] <- cluster1

            # Renumber to consecutive integers
            z <- as.integer(as.factor(z))
            currentK <- max(z)
        }

        return(z)
    }
}
