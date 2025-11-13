# Clustering Quality Diagnostics for Celda

#' @title Calculate cluster quality metrics
#' @description Computes silhouette scores, cluster separation, and other
#'  quality metrics to assess clustering results.
#' @param celdaMod A celda model object (celda_C, celda_G, or celda_CG)
#' @param useAssay String. Name of assay to use for calculations. Default "counts".
#' @param maxCells Integer. Maximum number of cells to use for silhouette
#'  calculation (for computational efficiency). If dataset has more cells,
#'  a random sample will be used. Default 5000.
#' @return A celdaDiagnostics object containing quality metrics
#' @export
#' @examples
#' data(celdaCMod)
#' result <- celda_C(celdaCMod$counts, K = 5, nchains = 1)
#' diag <- celdaClusterQuality(result)
#' print(diag)
#' plot(diag)
celdaClusterQuality <- function(celdaMod, useAssay = "counts", maxCells = 5000) {
  # Determine model type and extract appropriate data
  if (is(celdaMod, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(celdaMod, useAssay)

    if (is(celdaMod, "celda_C") ||
        !is.null(SummarizedExperiment::colData(celdaMod)$celda_cell_cluster)) {
      z <- as.integer(SummarizedExperiment::colData(celdaMod)$celda_cell_cluster)
      K <- length(unique(z))
      clusterType <- "cell"
    } else {
      stop("Could not determine clustering type from SingleCellExperiment object")
    }
  } else if (is(celdaMod, "celda_C")) {
    # Extract from celda_C S4 object
    stop("Please provide a SingleCellExperiment object with celda results")
  } else {
    stop("celdaMod must be a SingleCellExperiment with celda results")
  }

  # Subsample if too large
  if (ncol(counts) > maxCells) {
    idx <- sample(ncol(counts), maxCells)
    counts <- counts[, idx, drop = FALSE]
    z <- z[idx]
  }

  # Calculate silhouette scores
  silScores <- .calculateSilhouetteScores(counts, z, K)

  # Calculate cluster separation
  separation <- .calculateClusterSeparation(counts, z, K)

  # Calculate within-cluster sum of squares
  wcss <- .calculateWCSS(counts, z, K)

  # Calculate cluster stability (if available)
  stability <- NA  # Requires per-chain results

  results <- list(
    silhouette = silScores,
    separation = separation,
    wcss = wcss,
    stability = stability,
    summary = .summarizeDiagnostics(silScores, separation, wcss, stability)
  )

  class(results) <- "celdaDiagnostics"
  return(results)
}


#' @keywords internal
.calculateSilhouetteScores <- function(counts, z, K) {
  # Use normalized counts for distance calculation
  normCounts <- t(scale(t(as.matrix(counts))))
  normCounts[is.na(normCounts)] <- 0  # Handle constant features

  # Calculate pairwise distances (sample if too large)
  if (ncol(counts) > 1000) {
    # Use approximate silhouette with sampling
    sampleSize <- 1000
    idx <- sample(ncol(counts), sampleSize)
    distMat <- as.matrix(stats::dist(t(normCounts[, idx])))
    zSample <- z[idx]
    sil <- cluster::silhouette(zSample, distMat)
  } else {
    distMat <- as.matrix(stats::dist(t(normCounts)))
    sil <- cluster::silhouette(z, distMat)
  }

  return(sil)
}


#' @keywords internal
.calculateClusterSeparation <- function(counts, z, K) {
  # Calculate mean expression for each cluster
  clusterMeans <- vapply(seq(K), function(k) {
    if (sum(z == k) > 0) {
      Matrix::rowMeans(counts[, z == k, drop = FALSE])
    } else {
      rep(0, nrow(counts))
    }
  }, numeric(nrow(counts)))

  # Calculate pairwise distances between cluster centroids
  if (K > 1) {
    centroidDist <- as.matrix(stats::dist(t(clusterMeans)))

    # Return min, mean, and max separation
    minSep <- min(centroidDist[upper.tri(centroidDist)])
    meanSep <- mean(centroidDist[upper.tri(centroidDist)])
    maxSep <- max(centroidDist[upper.tri(centroidDist)])
  } else {
    minSep <- NA
    meanSep <- NA
    maxSep <- NA
  }

  return(list(min = minSep, mean = meanSep, max = maxSep))
}


#' @keywords internal
.calculateWCSS <- function(counts, z, K) {
  # Within-cluster sum of squares
  wcss <- vapply(seq(K), function(k) {
    clusterCells <- which(z == k)
    if (length(clusterCells) == 0) return(0)

    clusterCounts <- as.matrix(counts[, clusterCells, drop = FALSE])
    centroid <- Matrix::rowMeans(counts[, clusterCells, drop = FALSE])

    sum((clusterCounts - centroid)^2)
  }, numeric(1))

  return(list(total = sum(wcss), byCluster = wcss))
}


#' @keywords internal
.summarizeDiagnostics <- function(silScores, separation, wcss, stability) {
  meanSil <- mean(silScores[, 3])

  quality <- "Good"
  if (meanSil < 0.25) {
    quality <- "Poor"
  } else if (meanSil < 0.5) {
    quality <- "Fair"
  }

  return(list(
    meanSilhouette = meanSil,
    minSeparation = separation$min,
    totalWCSS = wcss$total,
    overallQuality = quality
  ))
}


#' @title Print celda diagnostics
#' @param x celdaDiagnostics object
#' @param ... Additional arguments (ignored)
#' @export
print.celdaDiagnostics <- function(x, ...) {
  cat("Celda Clustering Diagnostics\n")
  cat("============================\n\n")
  cat(sprintf("Overall Quality: %s\n", x$summary$overallQuality))
  cat(sprintf("Mean Silhouette Score: %.3f\n", x$summary$meanSilhouette))

  if (!is.na(x$summary$minSeparation)) {
    cat(sprintf("Min Cluster Separation: %.3f\n", x$summary$minSeparation))
  }

  cat(sprintf("Total WCSS: %.2e\n", x$summary$totalWCSS))
  cat("\nInterpretation:\n")
  cat("  Silhouette > 0.5: Good separation\n")
  cat("  Silhouette 0.25-0.5: Moderate separation\n")
  cat("  Silhouette < 0.25: Poor separation or overlapping clusters\n")

  invisible(x)
}


#' @title Plot celda diagnostics
#' @param x celdaDiagnostics object
#' @param ... Additional arguments (ignored)
#' @export
plot.celdaDiagnostics <- function(x, ...) {
  # Save current par settings
  oldPar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldPar))

  # Create diagnostic plots
  graphics::par(mfrow = c(2, 2))

  # 1. Silhouette plot
  cluster::plot.silhouette(x$silhouette, main = "Silhouette Plot")

  # 2. Cluster size distribution
  clusterSizes <- table(x$silhouette[, 1])
  graphics::barplot(clusterSizes, main = "Cluster Sizes",
          xlab = "Cluster", ylab = "Number of Cells",
          col = "steelblue")

  # 3. Silhouette by cluster
  silByCluster <- tapply(x$silhouette[, 3], x$silhouette[, 1], mean)
  graphics::barplot(silByCluster, main = "Mean Silhouette by Cluster",
          xlab = "Cluster", ylab = "Mean Silhouette Score",
          col = ifelse(silByCluster > 0.5, "green3",
                       ifelse(silByCluster > 0.25, "orange", "red")))
  graphics::abline(h = 0.25, lty = 2, col = "gray40")
  graphics::abline(h = 0.5, lty = 2, col = "gray40")

  # 4. WCSS by cluster
  graphics::barplot(x$wcss$byCluster, main = "Within-Cluster Sum of Squares",
          xlab = "Cluster", ylab = "WCSS",
          col = "coral")

  graphics::par(mfrow = c(1, 1))
  invisible(x)
}


#' @title Calculate convergence metrics
#' @description Analyzes log-likelihood history to assess convergence
#' @param celdaMod A celda model object with log-likelihood history
#' @return List with convergence assessment
#' @keywords internal
.assessConvergence <- function(celdaMod) {
  ll <- celdaMod@completeLogLik

  if (length(ll) < 10) {
    return(list(
      converged = NA,
      message = "Insufficient iterations to assess convergence"
    ))
  }

  # Check if log-likelihood is still increasing
  last10 <- tail(ll, 10)
  trend <- stats::coef(stats::lm(last10 ~ seq_along(last10)))[2]

  # Calculate improvement in last 10 iterations
  improvement <- (last10[10] - last10[1]) / abs(last10[1])

  converged <- (trend < 1e-6 && improvement < 0.001)

  return(list(
    converged = converged,
    trend = trend,
    improvement = improvement,
    message = if (converged) {
      "Model appears to have converged"
    } else {
      "Model may not have fully converged - consider more iterations"
    }
  ))
}
