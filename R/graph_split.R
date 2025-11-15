# Graph-based split candidate identification functions
# These functions enhance subcluster identification using graph theory
# and topological analysis of cell relationships.


#' Calculate Newman-Girvan Modularity
#'
#' Calculates graph modularity score using Newman-Girvan formula.
#' Higher modularity indicates stronger community structure.
#'
#' @param adjMatrix Numeric matrix. Adjacency matrix from kNN graph (cells x cells).
#'   Values represent edge weights (e.g., 1 for connected, 0 for disconnected).
#' @return Numeric. Modularity score between 0 and 1. Higher values indicate
#'   stronger community structure. Returns 0 for degenerate cases (empty graph,
#'   all disconnected, or single community).
#' @keywords internal
#' @details
#' Newman-Girvan modularity Q is calculated as:
#' Q = (1/2m) * sum_ij [A_ij - (k_i * k_j)/(2m)] * delta(c_i, c_j)
#' where:
#' - A_ij is the adjacency matrix
#' - k_i, k_j are node degrees
#' - m is the total number of edges
#' - delta(c_i, c_j) = 1 if nodes i and j are in the same community
#'
#' For split candidate identification, we use a simple two-community detection
#' based on spectral clustering of the adjacency matrix.
.calculateModularity <- function(adjMatrix) {
    # Input validation
    if (!is.matrix(adjMatrix) && !inherits(adjMatrix, "Matrix")) {
        stop("adjMatrix must be a matrix or sparse Matrix")
    }

    n <- nrow(adjMatrix)

    # Handle edge cases
    if (n < 2) {
        return(0)
    }

    # Convert to dense matrix if sparse
    if (inherits(adjMatrix, "Matrix")) {
        adjMatrix <- as.matrix(adjMatrix)
    }

    # Calculate total edges (2m in Newman-Girvan notation)
    m2 <- sum(adjMatrix)

    # Empty graph or all disconnected
    if (m2 == 0) {
        return(0)
    }

    # Calculate node degrees
    degrees <- rowSums(adjMatrix)

    # All nodes isolated
    if (all(degrees == 0)) {
        return(0)
    }

    # Simple two-community detection using spectral clustering
    # Calculate Laplacian matrix
    D <- diag(degrees)
    L <- D - adjMatrix

    # Handle numerical issues
    L[is.na(L)] <- 0
    L[is.infinite(L)] <- 0

    # Eigendecomposition (use tryCatch for numerical stability)
    eigenResult <- tryCatch({
        eigen(L, symmetric = TRUE)
    }, error = function(e) {
        return(NULL)
    })

    if (is.null(eigenResult)) {
        return(0)
    }

    # Find second smallest eigenvalue (Fiedler vector)
    # eigenvalues are in decreasing order, so we want second-to-last
    sortedIdx <- order(eigenResult$values)
    fiedlerVector <- eigenResult$vectors[, sortedIdx[2]]

    # Assign communities based on sign of Fiedler vector
    communities <- as.integer(fiedlerVector > 0) + 1L

    # Check if we got a valid partition
    if (length(unique(communities)) < 2) {
        # All in same community - no structure
        return(0)
    }

    # Calculate modularity Q
    Q <- 0
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            if (communities[i] == communities[j]) {
                expected <- (degrees[i] * degrees[j]) / m2
                Q <- Q + (adjMatrix[i, j] - expected)
            }
        }
    }
    Q <- Q / m2

    # Ensure Q is in valid range [0, 1]
    Q <- max(0, min(1, Q))

    return(Q)
}


#' Find Bimodal Genes Using Hartigan's Dip Test
#'
#' Identifies genes with bimodal expression distributions within a cluster.
#' Bimodality suggests potential subclusters.
#'
#' @param clusterCounts Numeric matrix. Expression counts for cells in the cluster
#'   (genes x cells).
#' @param pvalueThreshold Numeric. P-value threshold for significance. Default 0.05.
#' @return Integer vector. Indices of genes with significant bimodal distributions.
#' @keywords internal
#' @details
#' Uses Hartigan's dip test to detect multimodality. The dip test measures
#' the maximum difference between the empirical distribution function and
#' the best fitting unimodal distribution.
#'
#' If the 'diptest' package is available, it uses the optimized implementation.
#' Otherwise, it falls back to a basic implementation that checks for local
#' minima in the density estimate.
#'
#' Genes must have at least 10 non-zero values to be tested.
.findBimodalGenes <- function(clusterCounts, pvalueThreshold = 0.05) {
    # Input validation
    if (!is.matrix(clusterCounts) && !inherits(clusterCounts, "Matrix")) {
        stop("clusterCounts must be a matrix or sparse Matrix")
    }

    nGenes <- nrow(clusterCounts)
    nCells <- ncol(clusterCounts)

    # Handle edge cases
    if (nGenes == 0 || nCells < 10) {
        return(integer(0))
    }

    # Convert to dense matrix if sparse
    if (inherits(clusterCounts, "Matrix")) {
        clusterCounts <- as.matrix(clusterCounts)
    }

    # Check if diptest package is available
    hasDiptest <- requireNamespace("diptest", quietly = TRUE)

    bimodalGenes <- integer(0)

    for (i in seq_len(nGenes)) {
        geneExpression <- clusterCounts[i, ]

        # Require minimum non-zero values
        nonZeroValues <- geneExpression[geneExpression > 0]
        if (length(nonZeroValues) < 10) {
            next
        }

        # Perform dip test
        if (hasDiptest) {
            # Use optimized diptest package
            dipResult <- tryCatch({
                diptest::dip.test(nonZeroValues)
            }, error = function(e) {
                return(NULL)
            })

            if (!is.null(dipResult) && dipResult$p.value < pvalueThreshold) {
                bimodalGenes <- c(bimodalGenes, i)
            }
        } else {
            # Basic implementation: check for local minima in density
            # Standardize values
            if (sd(nonZeroValues) == 0) {
                next
            }

            standardized <- scale(nonZeroValues)

            # Create density estimate
            densityResult <- tryCatch({
                density(standardized, n = 512)
            }, error = function(e) {
                return(NULL)
            })

            if (is.null(densityResult)) {
                next
            }

            # Find local minima in density
            densityValues <- densityResult$y
            n <- length(densityValues)

            # A local minimum is where density[i-1] > density[i] < density[i+1]
            localMinima <- which(
                densityValues[2:(n-1)] < densityValues[1:(n-2)] &
                densityValues[2:(n-1)] < densityValues[3:n]
            )

            # Check if there's a significant valley (local minimum)
            # between two peaks
            if (length(localMinima) > 0) {
                # Find peaks (local maxima)
                localMaxima <- which(
                    densityValues[2:(n-1)] > densityValues[1:(n-2)] &
                    densityValues[2:(n-1)] > densityValues[3:n]
                )

                # Check if we have at least 2 peaks with a valley between
                if (length(localMaxima) >= 2) {
                    # Check if valley is substantial (depth > 20% of max peak)
                    maxPeak <- max(densityValues[localMaxima + 1])
                    minValley <- min(densityValues[localMinima + 1])

                    if ((maxPeak - minValley) / maxPeak > 0.2) {
                        bimodalGenes <- c(bimodalGenes, i)
                    }
                }
            }
        }
    }

    return(bimodalGenes)
}


#' Detect Community Substructure in Correlation Matrix
#'
#' Detects community structure by analyzing correlation-based graph connectivity.
#'
#' @param corrMatrix Numeric matrix. Cell-cell correlation matrix (cells x cells).
#' @param threshold Numeric. Correlation threshold for edge creation. Default 0.3.
#'   Edges are created between cells with correlation > threshold.
#' @return Numeric. Substructure score between 0 and 1. Higher values indicate
#'   more pronounced community structure (multiple disconnected components or
#'   weak connectivity).
#' @keywords internal
#' @details
#' Builds a graph where edges connect cells with correlation above the threshold.
#' Calculates the number of disconnected components and the connectivity structure.
#'
#' Returns a score based on:
#' - Number of disconnected components (more = higher score)
#' - Relative sizes of components (balanced sizes = higher score)
#' - Overall connectivity (sparse graph = higher score)
#'
#' Score of 0 indicates fully connected graph (no substructure).
#' Score near 1 indicates clear community separation.
.detectSubstructure <- function(corrMatrix, threshold = 0.3) {
    # Input validation
    if (!is.matrix(corrMatrix) && !inherits(corrMatrix, "Matrix")) {
        stop("corrMatrix must be a matrix or sparse Matrix")
    }

    n <- nrow(corrMatrix)

    # Handle edge cases
    if (n < 2) {
        return(0)
    }

    # Convert to dense matrix if sparse
    if (inherits(corrMatrix, "Matrix")) {
        corrMatrix <- as.matrix(corrMatrix)
    }

    # Build adjacency matrix based on threshold
    # Remove self-loops (diagonal)
    diag(corrMatrix) <- 0
    adjMatrix <- (corrMatrix > threshold) * 1

    # Calculate connectivity
    totalPossibleEdges <- n * (n - 1) / 2
    actualEdges <- sum(adjMatrix) / 2  # Divide by 2 since symmetric

    # If fully connected or fully disconnected
    if (actualEdges == 0) {
        # All disconnected - maximum substructure
        return(1)
    }
    if (actualEdges == totalPossibleEdges) {
        # Fully connected - no substructure
        return(0)
    }

    # Find connected components using BFS
    visited <- rep(FALSE, n)
    components <- rep(0L, n)
    componentId <- 0L

    for (start in seq_len(n)) {
        if (!visited[start]) {
            componentId <- componentId + 1L
            # BFS to find all nodes in this component
            queue <- start
            visited[start] <- TRUE
            components[start] <- componentId

            while (length(queue) > 0) {
                current <- queue[1]
                queue <- queue[-1]

                # Find neighbors
                neighbors <- which(adjMatrix[current, ] > 0)
                for (neighbor in neighbors) {
                    if (!visited[neighbor]) {
                        visited[neighbor] <- TRUE
                        components[neighbor] <- componentId
                        queue <- c(queue, neighbor)
                    }
                }
            }
        }
    }

    nComponents <- componentId

    # Calculate substructure score
    if (nComponents == 1) {
        # Single component - score based on sparsity
        sparsity <- 1 - (actualEdges / totalPossibleEdges)
        return(sparsity * 0.5)  # Max 0.5 for single component
    } else {
        # Multiple components - higher score
        # Component size balance (more balanced = higher score)
        componentSizes <- tabulate(components, nComponents)
        sizeProportion <- componentSizes / n

        # Entropy-based balance measure (max when all equal size)
        entropy <- -sum(sizeProportion * log(sizeProportion + 1e-10))
        maxEntropy <- log(nComponents)
        balance <- entropy / maxEntropy

        # Component count score (more components = higher score, but diminishing)
        componentScore <- min(1, nComponents / 5)

        # Combined score
        score <- 0.5 + 0.3 * componentScore + 0.2 * balance
        return(min(1, score))
    }
}


#' Identify Split Candidates Using Graph-Based Methods
#'
#' Main function to identify cluster split candidates by combining statistical
#' heterogeneity with graph-based community structure analysis.
#'
#' @param counts Integer matrix. Gene expression counts (genes x cells).
#' @param z Integer vector. Current cell cluster assignments.
#' @param K Integer. Number of clusters.
#' @param reducedDim Numeric matrix. Optional reduced dimensional representation
#'   (cells x dimensions), such as UMAP or t-SNE coordinates. If NULL, uses
#'   correlation-based substructure detection. Default NULL.
#' @param minCell Integer. Minimum number of cells required to consider a cluster
#'   for splitting. Default 3.
#' @param heterogeneityThreshold Numeric. Quantile threshold for heterogeneity
#'   filtering (0-1). Default 0.3 means top 30% most heterogeneous clusters.
#' @return Integer vector. Cluster IDs that are candidates for splitting.
#' @keywords internal
#' @details
#' Combines three complementary metrics:
#'
#' 1. Statistical Heterogeneity (40% weight):
#'    - Coefficient of variation of cell total counts
#'    - Gene expression variance
#'
#' 2. Graph Structure (40% weight):
#'    - If reducedDim provided: builds kNN graph and calculates modularity
#'    - Otherwise: uses correlation-based substructure detection
#'
#' 3. Bimodal Genes (20% weight):
#'    - Proportion of genes with bimodal distributions
#'
#' Final score is weighted combination. Clusters above the heterogeneityThreshold
#' quantile are returned as split candidates.
#'
#' For large clusters (>1000 cells), samples cells to reduce computational cost.
.identifySplitCandidates_GraphBased <- function(counts,
                                                 z,
                                                 K,
                                                 reducedDim = NULL,
                                                 minCell = 3,
                                                 heterogeneityThreshold = 0.3) {
    # Input validation
    if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
        stop("counts must be a matrix or sparse Matrix")
    }
    if (!is.integer(z)) {
        z <- as.integer(z)
    }
    if (length(z) != ncol(counts)) {
        stop("Length of z must equal number of columns in counts")
    }

    # Identify candidate clusters (size >= minCell)
    zTa <- tabulate(z, K)
    zCandidate <- which(zTa >= minCell)

    if (length(zCandidate) == 0) {
        return(integer(0))
    }

    # Calculate scores for each candidate cluster
    combinedScores <- vapply(zCandidate, function(k) {
        clusterCells <- which(z == k)
        nCellsInCluster <- length(clusterCells)

        if (nCellsInCluster < minCell) {
            return(0)
        }

        # === 1. Statistical Heterogeneity (existing method) ===
        # Coefficient of variation of cell total counts
        cellTotals <- Matrix::colSums(counts[, clusterCells, drop = FALSE])
        if (length(cellTotals) < 2) {
            statScore <- 0
        } else {
            cv <- sd(cellTotals) / mean(cellTotals)

            # Gene expression variance
            if (nrow(counts) > 100) {
                # Sample genes for speed
                geneSample <- sample(nrow(counts), min(100, nrow(counts)))
                geneVar <- mean(apply(
                    as.matrix(counts[geneSample, clusterCells, drop = FALSE]),
                    1, var
                ))
            } else {
                geneVar <- mean(apply(
                    as.matrix(counts[, clusterCells, drop = FALSE]),
                    1, var
                ))
            }

            # Combined statistical heterogeneity
            statScore <- cv + log1p(geneVar)
        }

        # === 2. Graph-Based Substructure Detection ===
        graphScore <- 0

        # Sample cells if cluster is very large
        if (nCellsInCluster > 1000) {
            sampleIdx <- sample(clusterCells, 1000)
        } else {
            sampleIdx <- clusterCells
        }

        if (!is.null(reducedDim)) {
            # Use reduced dimensions for kNN graph
            if (nrow(reducedDim) == ncol(counts)) {
                clusterCoords <- reducedDim[sampleIdx, , drop = FALSE]

                # Build kNN graph (k=10 neighbors)
                k <- min(10, nrow(clusterCoords) - 1)
                if (k >= 2) {
                    # Calculate pairwise distances
                    distMatrix <- as.matrix(dist(clusterCoords))

                    # Build kNN adjacency matrix
                    adjMatrix <- matrix(0, nrow = nrow(distMatrix),
                                       ncol = ncol(distMatrix))
                    for (i in seq_len(nrow(distMatrix))) {
                        # Find k nearest neighbors (excluding self)
                        neighbors <- order(distMatrix[i, ])[2:(k + 1)]
                        adjMatrix[i, neighbors] <- 1
                    }

                    # Make symmetric
                    adjMatrix <- (adjMatrix + t(adjMatrix)) > 0
                    adjMatrix <- adjMatrix * 1

                    # Calculate modularity
                    graphScore <- .calculateModularity(adjMatrix)
                }
            }
        } else {
            # Use correlation-based substructure detection
            if (length(sampleIdx) >= 10) {
                clusterCounts <- counts[, sampleIdx, drop = FALSE]

                # Calculate cell-cell correlation
                # Transpose so cells are rows
                countsTrans <- t(as.matrix(clusterCounts))

                # Only use top variable genes for speed
                if (nrow(clusterCounts) > 500) {
                    geneVar <- apply(clusterCounts, 1, var)
                    topGenes <- order(geneVar, decreasing = TRUE)[1:500]
                    countsTrans <- countsTrans[, topGenes, drop = FALSE]
                }

                # Calculate correlation matrix
                corrMatrix <- tryCatch({
                    cor(t(countsTrans))
                }, error = function(e) {
                    return(NULL)
                })

                if (!is.null(corrMatrix)) {
                    graphScore <- .detectSubstructure(corrMatrix, threshold = 0.3)
                }
            }
        }

        # === 3. Bimodal Gene Detection ===
        bimodalScore <- 0
        if (nCellsInCluster >= 10) {
            clusterCounts <- counts[, clusterCells, drop = FALSE]
            bimodalGenes <- .findBimodalGenes(clusterCounts, pvalueThreshold = 0.05)

            # Proportion of genes that are bimodal
            bimodalScore <- length(bimodalGenes) / nrow(counts)
        }

        # === 4. Combined Score ===
        # Normalize scores to similar ranges
        statScoreNorm <- min(1, statScore / 10)  # Typical range 0-10
        graphScoreNorm <- graphScore  # Already 0-1
        bimodalScoreNorm <- min(1, bimodalScore * 5)  # Boost bimodal signal

        # Weighted combination
        combinedScore <- (
            0.4 * statScoreNorm +
            0.4 * graphScoreNorm +
            0.2 * bimodalScoreNorm
        )

        return(combinedScore)
    }, numeric(1))

    # Filter candidates by heterogeneity threshold
    if (length(combinedScores) > 0 && any(!is.na(combinedScores))) {
        threshold <- quantile(combinedScores[!is.na(combinedScores)],
                            probs = 1 - heterogeneityThreshold,
                            na.rm = TRUE)
        highHeterogeneity <- !is.na(combinedScores) & combinedScores > threshold
        return(zCandidate[highHeterogeneity])
    } else {
        return(zCandidate)
    }
}
