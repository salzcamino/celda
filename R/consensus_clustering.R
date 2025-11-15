#' @title Build co-occurrence matrix from multiple chains
#' @description Builds a symmetric co-occurrence matrix where entry (i,j) represents
#'   the proportion of chains in which elements i and j were assigned to the same
#'   cluster.
#' @param chainResults List of chain results, each containing clustering assignments
#' @param type Character. Either "z" for cells or "y" for genes.
#' @return Sparse matrix of co-occurrence probabilities
#' @keywords internal
.buildCooccurrenceMatrix <- function(chainResults, type = c("z", "y")) {
    type <- match.arg(type)

    nChains <- length(chainResults)
    if (nChains == 0) {
        stop("chainResults must contain at least one chain")
    }

    # Extract cluster assignments from each chain
    assignments <- lapply(chainResults, function(res) {
        if (type == "z") {
            return(res$z)
        } else {
            return(res$y)
        }
    })

    # Check all chains have same length
    n <- length(assignments[[1]])
    if (!all(sapply(assignments, length) == n)) {
        stop("All chains must have same number of elements")
    }

    # Build co-occurrence matrix
    # For memory efficiency with large n, we use sparse matrix
    cooccur <- matrix(0, nrow = n, ncol = n)

    for (i in seq_len(nChains)) {
        z <- assignments[[i]]
        # For each pair of elements, increment if in same cluster
        # Vectorized approach: build indicator matrix
        for (k in unique(z)) {
            idx <- which(z == k)
            if (length(idx) > 1) {
                # All pairs within this cluster co-occur
                cooccur[idx, idx] <- cooccur[idx, idx] + 1
            }
        }
    }

    # Normalize by number of chains
    cooccur <- cooccur / nChains

    # Ensure diagonal is 1 (element always clusters with itself)
    diag(cooccur) <- 1

    return(cooccur)
}


#' @title Consensus clustering from multiple chains
#' @description Combines clustering results from multiple chains using either
#'   co-occurrence based hierarchical clustering or median assignment method.
#' @param allChainResults List of chain results
#' @param method Character. Either "cooccurrence" or "median"
#' @param minAgreement Numeric. Minimum agreement threshold for flagging
#'   low-confidence elements (default 0.7)
#' @param type Character. Either "z" for cells or "y" for genes
#' @return List containing:
#'   - assignments: Consensus cluster assignments
#'   - confidence: Per-element confidence scores
#'   - lowConfidenceIndices: Indices of low-confidence elements
#' @keywords internal
.consensusClustering <- function(allChainResults,
                                  method = c("cooccurrence", "median"),
                                  minAgreement = 0.7,
                                  type = c("z", "y")) {
    method <- match.arg(method)
    type <- match.arg(type)

    nChains <- length(allChainResults)
    if (nChains == 0) {
        stop("allChainResults must contain at least one chain")
    }

    # Extract assignments
    assignments <- lapply(allChainResults, function(res) {
        if (type == "z") {
            return(res$z)
        } else {
            return(res$y)
        }
    })

    n <- length(assignments[[1]])

    if (method == "cooccurrence") {
        # Build co-occurrence matrix
        cooccur <- .buildCooccurrenceMatrix(allChainResults, type = type)

        # Convert to distance: dist = 1 - cooccurrence
        dist <- 1 - cooccur
        diag(dist) <- 0  # Distance to self is 0

        # Convert to dist object (only upper triangle needed)
        distVec <- as.dist(dist)

        # Hierarchical clustering
        hc <- stats::hclust(distVec, method = "average")

        # Determine K from median across chains
        Kvalues <- sapply(allChainResults, function(res) {
            if (type == "z") {
                return(length(unique(res$z)))
            } else {
                return(length(unique(res$y)))
            }
        })
        K <- round(stats::median(Kvalues))

        # Cut tree at K clusters
        consensusAssignments <- stats::cutree(hc, k = K)

        # Calculate confidence as mean co-occurrence within assigned cluster
        confidence <- numeric(n)
        for (i in seq_len(n)) {
            # Find all elements in same consensus cluster
            sameCluster <- which(consensusAssignments == consensusAssignments[i])
            # Confidence = mean co-occurrence with cluster members
            if (length(sameCluster) > 1) {
                confidence[i] <- mean(cooccur[i, sameCluster])
            } else {
                confidence[i] <- 1  # Singleton cluster
            }
        }

    } else if (method == "median") {
        # For each element, find mode (most common) cluster assignment
        consensusAssignments <- integer(n)
        confidence <- numeric(n)

        for (i in seq_len(n)) {
            # Get assignments for this element across all chains
            elemAssignments <- sapply(assignments, function(a) a[i])

            # Find mode
            tab <- table(elemAssignments)
            mode <- as.integer(names(tab)[which.max(tab)])

            consensusAssignments[i] <- mode

            # Confidence = proportion of chains agreeing with mode
            confidence[i] <- max(tab) / nChains
        }

        # Renumber clusters to 1:K
        uniqueClusters <- unique(consensusAssignments)
        mapping <- seq_along(uniqueClusters)
        names(mapping) <- as.character(uniqueClusters)
        consensusAssignments <- mapping[as.character(consensusAssignments)]
    }

    # Identify low-confidence elements
    lowConfidenceIndices <- which(confidence < minAgreement)

    return(list(
        assignments = as.integer(consensusAssignments),
        confidence = confidence,
        lowConfidenceIndices = lowConfidenceIndices
    ))
}


#' @title Consensus clustering for celda_CG
#' @description Applies consensus clustering to both cell (z) and gene (y)
#'   cluster assignments from multiple chains.
#' @param allChainResults List of chain results from celda_CG
#' @param method Character. Either "cooccurrence" or "median"
#' @param minAgreement Numeric. Minimum agreement threshold (default 0.7)
#' @return List containing:
#'   - z: Consensus cell cluster assignments
#'   - y: Consensus gene module assignments
#'   - zConfidence: Per-cell confidence scores
#'   - yConfidence: Per-gene confidence scores
#'   - lowConfidenceCells: Indices of low-confidence cells
#'   - lowConfidenceGenes: Indices of low-confidence genes
#' @keywords internal
.consensusClustering_CG <- function(allChainResults,
                                     method = c("cooccurrence", "median"),
                                     minAgreement = 0.7) {
    method <- match.arg(method)

    if (length(allChainResults) == 0) {
        stop("allChainResults must contain at least one chain")
    }

    # Consensus for cell clusters (z)
    zConsensus <- .consensusClustering(
        allChainResults = allChainResults,
        method = method,
        minAgreement = minAgreement,
        type = "z"
    )

    # Consensus for gene modules (y)
    yConsensus <- .consensusClustering(
        allChainResults = allChainResults,
        method = method,
        minAgreement = minAgreement,
        type = "y"
    )

    return(list(
        z = zConsensus$assignments,
        y = yConsensus$assignments,
        zConfidence = zConsensus$confidence,
        yConfidence = yConsensus$confidence,
        lowConfidenceCells = zConsensus$lowConfidenceIndices,
        lowConfidenceGenes = yConsensus$lowConfidenceIndices
    ))
}
