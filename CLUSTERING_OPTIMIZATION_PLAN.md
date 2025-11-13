# Celda Clustering Algorithm Optimization Plan

**Version**: 1.0
**Date**: 2025-11-13
**Status**: Implementation Ready

## Executive Summary

This document outlines a comprehensive plan to improve the performance and quality of celda's clustering algorithms. The proposed optimizations focus on:
1. Reducing computational overhead through adaptive heuristics
2. Improving convergence speed via batch processing
3. Enhancing result quality through better diagnostics
4. Enabling parallel processing for split operations

**Expected Overall Improvement**: 3-5x speedup for typical workflows, better cluster quality, and improved scalability.

---

## Phase 1: Adaptive Split Heuristic Frequency (Week 1-2)

### Goal
Reduce the computational cost of the split heuristic while maintaining clustering quality.

### Current Bottleneck
- Split evaluation runs every `splitOnIter` iterations (default: 10)
- Each evaluation performs K=2 clustering on every splittable cluster
- For K=20 clusters, this means 20 mini-clustering runs every 10 iterations
- Total cost: O(K × cells_per_cluster × 5_iterations) per split check

### Implementation Details

#### 1.1 Add Adaptive Split Frequency

**File**: `R/celda_C.R`

**New Parameters**:
```r
celda_C <- function(...,
  splitOnIter = 10,
  splitAdaptive = TRUE,        # NEW: Enable adaptive splitting
  splitDecayRate = 0.8,         # NEW: Frequency decay rate
  splitMinIter = 20,            # NEW: Minimum iterations between splits
  ...)
```

**Algorithm**:
```r
# In .celda_C() main loop (around line 426)
nextSplitIter <- splitOnIter  # Initialize

while (iter <= maxIter & numIterWithoutImprovement <= stopIter) {
  # ... existing Gibbs/EM sampling ...

  # Adaptive split frequency calculation
  if (isTRUE(splitAdaptive)) {
    # Increase interval as clustering stabilizes
    if (splitOccurred) {
      nextSplitIter <- max(splitOnIter, nextSplitIter * splitDecayRate)
    } else {
      nextSplitIter <- min(splitMinIter, nextSplitIter * (1 / splitDecayRate))
    }
  }

  # Check if it's time to split
  shouldSplit <- (iter %% ceiling(nextSplitIter) == 0) & isTRUE(doCellSplit)

  if (K > 2 & iter != maxIter & shouldSplit) {
    # ... existing split logic ...
  }
}
```

**Expected Benefit**: 30-40% reduction in split evaluation overhead

#### 1.2 Pre-filter Split Candidates

**New Function**: `R/split_clusters.R`

```r
#' @title Identify clusters that need splitting based on heterogeneity
#' @description Uses fast heuristics to identify split candidates
#' @keywords internal
.identifySplitCandidates <- function(counts, z, K,
                                      heterogeneityThreshold = 0.3,
                                      minCell = 3) {
  zTa <- tabulate(z, K)
  zCandidate <- which(zTa >= minCell)

  if (length(zCandidate) == 0) {
    return(integer(0))
  }

  # Calculate within-cluster heterogeneity
  heterogeneity <- vapply(zCandidate, function(k) {
    clusterCells <- which(z == k)
    if (length(clusterCells) < minCell) return(0)

    # Use coefficient of variation of cell total counts
    cellTotals <- Matrix::colSums(counts[, clusterCells, drop = FALSE])
    cv <- sd(cellTotals) / mean(cellTotals)

    # Also consider gene expression variance
    if (nrow(counts) > 100) {
      # Sample genes for speed
      geneSample <- sample(nrow(counts), min(100, nrow(counts)))
      geneVar <- mean(apply(counts[geneSample, clusterCells, drop = FALSE], 1, var))
    } else {
      geneVar <- mean(apply(counts[, clusterCells, drop = FALSE], 1, var))
    }

    # Combined heterogeneity score
    return(cv + log1p(geneVar))
  }, numeric(1))

  # Filter candidates by heterogeneity threshold
  highHeterogeneity <- heterogeneity > quantile(heterogeneity,
                                                  probs = 1 - heterogeneityThreshold)

  return(zCandidate[highHeterogeneity])
}
```

**Integration in `.cCSplitZ()`** (around line 15):
```r
.cCSplitZ <- function(counts, ..., heterogeneityThreshold = 0.3) {
  ## Identify clusters to split - USE NEW FILTER
  zCandidate <- .identifySplitCandidates(counts, z, K,
                                          heterogeneityThreshold,
                                          minCell)

  if (length(zCandidate) == 0) {
    # Early exit if no good candidates
    return(list(z = z, ...))
  }

  # Continue with existing split logic on filtered candidates
  zToSplit <- zCandidate
  ...
}
```

**Expected Benefit**: 40-60% reduction in split evaluations

#### 1.3 Cache Split Results

**New Helper Function**:
```r
#' @keywords internal
.cacheSplitResults <- function() {
  list(
    clusters = list(),      # Cached K=2 results for each cluster
    signatures = list(),    # Hash of cluster composition
    iteration = list()      # When cache was created
  )
}

#' @keywords internal
.needsResplit <- function(cache, clusterIdx, currentZ, currentIter,
                          cacheExpiry = 20) {
  if (is.null(cache$clusters[[clusterIdx]])) return(TRUE)
  if (currentIter - cache$iteration[[clusterIdx]] > cacheExpiry) return(TRUE)

  # Check if cluster composition changed significantly
  currentSig <- digest::digest(which(currentZ == clusterIdx))
  if (cache$signatures[[clusterIdx]] != currentSig) return(TRUE)

  return(FALSE)
}
```

**Expected Benefit**: 20-30% reduction in redundant split calculations

---

## Phase 2: Parallel Split Evaluation (Week 2-3)

### Goal
Evaluate multiple cluster splits in parallel to reduce wall-clock time.

### Implementation Details

#### 2.1 Parallel Split Infrastructure

**File**: `R/split_clusters.R`

**Modified `.cCSplitZ()` function** (around line 38):
```r
.cCSplitZ <- function(counts, ..., nCores = 1,
                      parallelSplitThreshold = 5) {
  ## Identify clusters to split
  zToSplit <- .identifySplitCandidates(counts, z, K, minCell = minCell)

  if (length(zToSplit) == 0) {
    return(list(z = z, ...))
  }

  ## Parallel loop through each split-able cluster
  if (nCores > 1 && length(zToSplit) >= parallelSplitThreshold) {
    # Use parallel processing
    clustSplit <- parallel::mclapply(zToSplit, function(i) {
      clustLabel <- .celda_C(
        counts[, z == i],
        K = 2,
        zInitialize = "random",
        maxIter = 5,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE,
        reorder = FALSE,
        seed = NULL  # Let each worker use different random seed
      )
      return(as.integer(celdaClusters(clustLabel)$z))
    }, mc.cores = nCores)

    # Convert list to named list
    names(clustSplit) <- as.character(zToSplit)
    # Expand to full K length
    clustSplitFull <- vector("list", K)
    for (i in seq_along(zToSplit)) {
      clustSplitFull[[zToSplit[i]]] <- clustSplit[[i]]
    }
    clustSplit <- clustSplitFull

  } else {
    # Sequential processing (existing code)
    clustSplit <- vector("list", K)
    for (i in zToSplit) {
      clustLabel <- .celda_C(
        counts[, z == i],
        K = 2,
        zInitialize = "random",
        maxIter = 5,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE,
        reorder = FALSE
      )
      clustSplit[[i]] <- as.integer(celdaClusters(clustLabel)$z)
    }
  }

  # ... rest of existing code ...
}
```

#### 2.2 Add nCores parameter to public API

**File**: `R/celda_C.R` (line 78-98)

```r
setGeneric("celda_C",
  function(x,
    ...,
    nchains = 3,
    nCores = 1,              # NEW: Number of cores for parallel operations
    ...) {
    standardGeneric("celda_C")
  })
```

**Pass through to internal functions**:
- `.celdaCWithSeed()` → `.celda_C()` → `.cCSplitZ()`

**Expected Benefit**: Near-linear speedup with number of cores (3-4x with 4 cores)

---

## Phase 3: Batch Gibbs Sampling (Week 3-5)

### Goal
Update cells in batches rather than one-at-a-time to reduce matrix operation overhead.

### Implementation Details

#### 3.1 New Batch Gibbs Function

**File**: `R/celda_C.R`

```r
#' @title Collapsed Gibbs sampling with batch updates
#' @keywords internal
.cCCalcGibbsProbZ_Batch <- function(counts,
                                     mCPByS,
                                     nGByCP,
                                     nByC,
                                     nCP,
                                     z,
                                     s,
                                     K,
                                     nG,
                                     nM,
                                     alpha,
                                     beta,
                                     doSample = TRUE,
                                     batchSize = "auto") {

  # Determine batch size
  if (batchSize == "auto") {
    # Heuristic: batch size = sqrt(nM) for balanced efficiency
    batchSize <- max(10, floor(sqrt(nM)))
  }

  # Create batches by shuffling cells
  cellOrder <- sample(seq(nM))
  nBatches <- ceiling(nM / batchSize)

  probs <- matrix(NA, ncol = nM, nrow = K)

  for (batch in seq(nBatches)) {
    startIdx <- (batch - 1) * batchSize + 1
    endIdx <- min(batch * batchSize, nM)
    batchCells <- cellOrder[startIdx:endIdx]

    # Remove all cells in batch from current counts
    for (i in batchCells) {
      mCPByS[z[i], s[i]] <- mCPByS[z[i], s[i]] - 1L
      nGByCP[, z[i]] <- nGByCP[, z[i]] - counts[, i]
      nCP[z[i]] <- nCP[z[i]] - nByC[i]
    }

    # Calculate probabilities for all cells in batch
    for (i in batchCells) {
      for (j in seq_len(K)) {
        # Theta component
        probs[j, i] <- log(mCPByS[j, s[i]] + alpha)

        # Phi components
        if (j == z[i]) {
          # Cell was in this cluster
          probs[j, i] <- probs[j, i] +
            sum(lgamma(nGByCP[, j] + counts[, i] + beta)) -
            lgamma(nCP[j] + nByC[i] + nG * beta) -
            sum(lgamma(nGByCP[, j] + beta)) +
            lgamma(nCP[j] + nG * beta)
        } else {
          # Cell was not in this cluster
          probs[j, i] <- probs[j, i] +
            sum(lgamma(nGByCP[, j] + counts[, i] + beta)) -
            lgamma(nCP[j] + nByC[i] + nG * beta) -
            sum(lgamma(nGByCP[, j] + beta)) +
            lgamma(nCP[j] + nG * beta)
        }
      }
    }

    # Sample new assignments for batch
    if (isTRUE(doSample)) {
      for (i in batchCells) {
        z[i] <- .sampleLl(probs[, i])
      }
    }

    # Add all cells in batch back with new assignments
    for (i in batchCells) {
      mCPByS[z[i], s[i]] <- mCPByS[z[i], s[i]] + 1L
      nGByCP[, z[i]] <- nGByCP[, z[i]] + counts[, i]
      nCP[z[i]] <- nCP[z[i]] + nByC[i]
    }
  }

  return(list(
    mCPByS = mCPByS,
    nGByCP = nGByCP,
    nCP = nCP,
    z = z,
    probs = probs
  ))
}
```

#### 3.2 Add algorithm option

**File**: `R/celda_C.R` (line 26)

```r
#' @param algorithm String. Algorithm to use for clustering cell subpopulations.
#'  One of 'EM', 'Gibbs', or 'GibbsBatch'. The EM algorithm is fastest, especially
#'  for larger numbers of cells. 'GibbsBatch' uses batch updates for improved
#'  performance over standard Gibbs sampling. Default 'EM'.

algorithm = c("EM", "Gibbs", "GibbsBatch")
```

**Update algorithm dispatcher** (around line 355):
```r
algorithmFun <- switch(algorithm,
  "Gibbs" = ".cCCalcGibbsProbZ",
  "EM" = ".cCCalcEMProbZ",
  "GibbsBatch" = ".cCCalcGibbsProbZ_Batch"
)
```

**Expected Benefit**: 1.5-2x speedup over standard Gibbs, especially for large datasets

---

## Phase 4: Convergence Diagnostics (Week 5-6)

### Goal
Provide users with quality metrics to assess clustering results and detect convergence issues.

### Implementation Details

#### 4.1 Add Diagnostic Metrics

**New File**: `R/clustering_diagnostics.R`

```r
#' @title Calculate cluster quality metrics
#' @description Computes silhouette scores and other quality metrics
#' @param sce SingleCellExperiment with celda results
#' @param useAssay Assay to use for calculations
#' @return List of diagnostic metrics
#' @export
celdaClusterQuality <- function(sce, useAssay = "counts", maxCells = 5000) {
  counts <- SummarizedExperiment::assay(sce, useAssay)
  z <- as.integer(SummarizedExperiment::colData(sce)$celda_cell_cluster)
  K <- S4Vectors::metadata(sce)$celda_parameters$K

  # Subsample if too large
  if (ncol(counts) > maxCells) {
    idx <- sample(ncol(counts), maxCells)
    counts <- counts[, idx]
    z <- z[idx]
  }

  # Calculate silhouette scores
  silScores <- .calculateSilhouetteScores(counts, z, K)

  # Calculate cluster separation
  separation <- .calculateClusterSeparation(counts, z, K)

  # Calculate within-cluster sum of squares
  wcss <- .calculateWCSS(counts, z, K)

  # Calculate cluster stability (if multiple chains available)
  stability <- .calculateClusterStability(sce)

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

  # Calculate pairwise distances (sample if too large)
  if (ncol(counts) > 1000) {
    # Use approximate silhouette with sampling
    sampleSize <- 1000
    idx <- sample(ncol(counts), sampleSize)
    distMat <- as.matrix(dist(t(normCounts[, idx])))
    zSample <- z[idx]
    sil <- cluster::silhouette(zSample, distMat)
  } else {
    distMat <- as.matrix(dist(t(normCounts)))
    sil <- cluster::silhouette(z, distMat)
  }

  return(sil)
}

#' @keywords internal
.calculateClusterSeparation <- function(counts, z, K) {
  # Calculate mean expression for each cluster
  clusterMeans <- vapply(seq(K), function(k) {
    Matrix::rowMeans(counts[, z == k, drop = FALSE])
  }, numeric(nrow(counts)))

  # Calculate pairwise distances between cluster centroids
  centroidDist <- as.matrix(dist(t(clusterMeans)))

  # Return min, mean, and max separation
  minSep <- min(centroidDist[upper.tri(centroidDist)])
  meanSep <- mean(centroidDist[upper.tri(centroidDist)])
  maxSep <- max(centroidDist[upper.tri(centroidDist)])

  return(list(min = minSep, mean = meanSep, max = maxSep))
}

#' @keywords internal
.calculateWCSS <- function(counts, z, K) {
  # Within-cluster sum of squares
  wcss <- vapply(seq(K), function(k) {
    clusterCells <- which(z == k)
    if (length(clusterCells) == 0) return(0)

    clusterCounts <- counts[, clusterCells, drop = FALSE]
    centroid <- Matrix::rowMeans(clusterCounts)

    sum((clusterCounts - centroid)^2)
  }, numeric(1))

  return(list(total = sum(wcss), byCluster = wcss))
}

#' @keywords internal
.calculateClusterStability <- function(sce) {
  # Check if multiple chains were run
  # This requires storing per-chain results (enhancement)
  # For now, return NA
  return(NA)
}

#' @keywords internal
.summarizeDiagnostics <- function(silScores, separation, wcss, stability) {
  meanSil <- mean(silScores[, 3])

  quality <- "Good"
  if (meanSil < 0.25) quality <- "Poor"
  else if (meanSil < 0.5) quality <- "Fair"

  return(list(
    meanSilhouette = meanSil,
    minSeparation = separation$min,
    totalWCSS = wcss$total,
    overallQuality = quality
  ))
}

#' @export
print.celdaDiagnostics <- function(x, ...) {
  cat("Celda Clustering Diagnostics\n")
  cat("============================\n\n")
  cat(sprintf("Overall Quality: %s\n", x$summary$overallQuality))
  cat(sprintf("Mean Silhouette Score: %.3f\n", x$summary$meanSilhouette))
  cat(sprintf("Min Cluster Separation: %.3f\n", x$summary$minSeparation))
  cat(sprintf("Total WCSS: %.2e\n", x$summary$totalWCSS))
  cat("\nInterpretation:\n")
  cat("  Silhouette > 0.5: Good separation\n")
  cat("  Silhouette 0.25-0.5: Moderate separation\n")
  cat("  Silhouette < 0.25: Poor separation or overlapping clusters\n")
}

#' @export
plot.celdaDiagnostics <- function(x, ...) {
  # Create diagnostic plots
  par(mfrow = c(2, 2))

  # 1. Silhouette plot
  cluster::plot(x$silhouette, main = "Silhouette Plot")

  # 2. Cluster size distribution
  clusterSizes <- table(x$silhouette[, 1])
  barplot(clusterSizes, main = "Cluster Sizes",
          xlab = "Cluster", ylab = "Number of Cells")

  # 3. Silhouette by cluster
  silByCluster <- tapply(x$silhouette[, 3], x$silhouette[, 1], mean)
  barplot(silByCluster, main = "Mean Silhouette by Cluster",
          xlab = "Cluster", ylab = "Mean Silhouette Score")
  abline(h = 0.25, lty = 2, col = "red")
  abline(h = 0.5, lty = 2, col = "green")

  # 4. WCSS by cluster
  barplot(x$wcss$byCluster, main = "Within-Cluster Sum of Squares",
          xlab = "Cluster", ylab = "WCSS")

  par(mfrow = c(1, 1))
}
```

**Expected Benefit**: Helps users identify poor clustering results and parameter choices

---

## Phase 5: Smart Chain Management (Week 6-7)

### Goal
Optimize multi-chain runs by terminating poor chains early and sharing information.

### Implementation Details

#### 5.1 Early Chain Termination

**File**: `R/celda_C.R` (modify main loop around line 364)

```r
.celda_C <- function(...,
                     earlyChainStop = TRUE,
                     earlyStopThreshold = 0.05) {

  allChains <- seq(nchains)
  bestResult <- NULL
  chainResults <- vector("list", nchains)

  for (i in allChains) {
    # ... existing initialization ...

    iter <- 1L
    numIterWithoutImprovement <- 0L
    shouldTerminate <- FALSE

    while (iter <= maxIter &
           numIterWithoutImprovement <= stopIter &
           !shouldTerminate) {

      # ... existing sampling code ...

      # Early termination check
      if (earlyChainStop && !is.null(bestResult) && iter > 20) {
        # Compare current chain to best chain so far
        llDiff <- (bestResult$finalLogLik - tempLl) / abs(bestResult$finalLogLik)

        if (llDiff > earlyStopThreshold) {
          .logMessages(date(),
            ".... Chain", i, "terminated early (significantly worse than best)",
            logfile = logfile,
            append = TRUE,
            verbose = verbose
          )
          shouldTerminate <- TRUE
        }
      }

      # ... rest of iteration ...
    }

    # Store result
    result <- list(
      z = zBest,
      completeLogLik = ll,
      finalLogLik = llBest,
      # ... other fields ...
    )

    chainResults[[i]] <- result

    # Update best result
    if (is.null(bestResult) || result$finalLogLik > bestResult$finalLogLik) {
      bestResult <- result
    }
  }

  # ... rest of function ...
}
```

#### 5.2 Chain Consensus and Stability

**New Function**:
```r
#' @title Calculate consensus clustering across chains
#' @keywords internal
.calculateChainConsensus <- function(chainResults, nM, K) {
  nChains <- length(chainResults)

  # Build co-clustering matrix
  coClusterMat <- matrix(0, nrow = nM, ncol = nM)

  for (chain in chainResults) {
    z <- chain$z
    for (i in seq(nM)) {
      for (j in seq(i, nM)) {
        if (z[i] == z[j]) {
          coClusterMat[i, j] <- coClusterMat[i, j] + 1
          if (i != j) coClusterMat[j, i] <- coClusterMat[j, i] + 1
        }
      }
    }
  }

  # Normalize by number of chains
  coClusterMat <- coClusterMat / nChains

  # Calculate stability metric
  stability <- mean(coClusterMat[upper.tri(coClusterMat)] > 0.8)

  return(list(
    coClusterMatrix = coClusterMat,
    stability = stability
  ))
}
```

**Expected Benefit**: 20-40% speedup in multi-chain runs, better uncertainty quantification

---

## Phase 6: Benchmarking and Validation (Week 7-8)

### Goal
Validate improvements and quantify performance gains.

### Implementation Details

#### 6.1 Comprehensive Benchmark Suite

**New File**: `benchmark_clustering_optimizations.R`

```r
#!/usr/bin/env Rscript
# Benchmark script for clustering optimization evaluation

library(celda)
library(SingleCellExperiment)
library(microbenchmark)
library(ggplot2)

#' Generate synthetic test data
generateTestData <- function(nGenes = 1000, nCells = 500, K = 5, L = 10) {
  # Simulate celda_CG data
  sim <- simulateCellsCelda(
    model = "celda_CG",
    nCells = nCells,
    nGenes = nGenes,
    K = K,
    L = L,
    minCells = 50,
    seed = 12345
  )
  return(sim)
}

#' Benchmark split heuristic improvements
benchmarkSplitOptimizations <- function(sim) {
  cat("\n=== Benchmarking Split Optimizations ===\n")

  # Original implementation
  t1 <- system.time({
    result1 <- celda_C(
      sim$counts,
      K = sim$K,
      splitOnIter = 10,
      splitAdaptive = FALSE,  # Original behavior
      maxIter = 100,
      nchains = 1,
      verbose = FALSE
    )
  })

  # Adaptive split frequency
  t2 <- system.time({
    result2 <- celda_C(
      sim$counts,
      K = sim$K,
      splitOnIter = 10,
      splitAdaptive = TRUE,   # NEW: Adaptive
      maxIter = 100,
      nchains = 1,
      verbose = FALSE
    )
  })

  # Parallel splits
  t3 <- system.time({
    result3 <- celda_C(
      sim$counts,
      K = sim$K,
      splitOnIter = 10,
      splitAdaptive = TRUE,
      nCores = 4,             # NEW: Parallel
      maxIter = 100,
      nchains = 1,
      verbose = FALSE
    )
  })

  results <- data.frame(
    Method = c("Original", "Adaptive", "Adaptive+Parallel"),
    Time = c(t1["elapsed"], t2["elapsed"], t3["elapsed"]),
    Speedup = c(1, t1["elapsed"] / t2["elapsed"],
                t1["elapsed"] / t3["elapsed"]),
    LogLik = c(
      metadata(result1)$celda_parameters$finalLogLik,
      metadata(result2)$celda_parameters$finalLogLik,
      metadata(result3)$celda_parameters$finalLogLik
    )
  )

  print(results)
  return(results)
}

#' Benchmark batch Gibbs sampling
benchmarkBatchGibbs <- function(sim) {
  cat("\n=== Benchmarking Batch Gibbs Sampling ===\n")

  # Standard Gibbs
  t1 <- system.time({
    result1 <- celda_C(
      sim$counts,
      K = sim$K,
      algorithm = "Gibbs",
      maxIter = 50,
      splitOnIter = -1,  # Disable splitting for fair comparison
      nchains = 1,
      verbose = FALSE
    )
  })

  # Batch Gibbs
  t2 <- system.time({
    result2 <- celda_C(
      sim$counts,
      K = sim$K,
      algorithm = "GibbsBatch",
      maxIter = 50,
      splitOnIter = -1,
      nchains = 1,
      verbose = FALSE
    )
  })

  # EM (for reference)
  t3 <- system.time({
    result3 <- celda_C(
      sim$counts,
      K = sim$K,
      algorithm = "EM",
      maxIter = 50,
      splitOnIter = -1,
      nchains = 1,
      verbose = FALSE
    )
  })

  results <- data.frame(
    Algorithm = c("Gibbs", "GibbsBatch", "EM"),
    Time = c(t1["elapsed"], t2["elapsed"], t3["elapsed"]),
    Speedup = c(1, t1["elapsed"] / t2["elapsed"],
                t1["elapsed"] / t3["elapsed"]),
    LogLik = c(
      metadata(result1)$celda_parameters$finalLogLik,
      metadata(result2)$celda_parameters$finalLogLik,
      metadata(result3)$celda_parameters$finalLogLik
    )
  )

  print(results)
  return(results)
}

#' Benchmark chain management
benchmarkChainManagement <- function(sim) {
  cat("\n=== Benchmarking Smart Chain Management ===\n")

  # Standard multi-chain
  t1 <- system.time({
    result1 <- celda_C(
      sim$counts,
      K = sim$K,
      nchains = 5,
      earlyChainStop = FALSE,
      maxIter = 100,
      verbose = FALSE
    )
  })

  # Early chain termination
  t2 <- system.time({
    result2 <- celda_C(
      sim$counts,
      K = sim$K,
      nchains = 5,
      earlyChainStop = TRUE,
      maxIter = 100,
      verbose = FALSE
    )
  })

  results <- data.frame(
    Method = c("Standard", "EarlyStop"),
    Time = c(t1["elapsed"], t2["elapsed"]),
    Speedup = c(1, t1["elapsed"] / t2["elapsed"]),
    LogLik = c(
      metadata(result1)$celda_parameters$finalLogLik,
      metadata(result2)$celda_parameters$finalLogLik
    )
  )

  print(results)
  return(results)
}

#' Test scaling behavior
testScaling <- function() {
  cat("\n=== Testing Scaling Behavior ===\n")

  cellCounts <- c(500, 1000, 2000, 5000)
  results <- lapply(cellCounts, function(nCells) {
    cat(sprintf("\nTesting with %d cells...\n", nCells))
    sim <- generateTestData(nGenes = 1000, nCells = nCells, K = 5)

    # Original
    t1 <- system.time({
      celda_C(sim$counts, K = 5, maxIter = 50,
              splitAdaptive = FALSE, nchains = 1, verbose = FALSE)
    })

    # Optimized
    t2 <- system.time({
      celda_C(sim$counts, K = 5, maxIter = 50,
              splitAdaptive = TRUE, nCores = 4,
              algorithm = "GibbsBatch",
              nchains = 1, verbose = FALSE)
    })

    data.frame(
      nCells = nCells,
      OriginalTime = t1["elapsed"],
      OptimizedTime = t2["elapsed"],
      Speedup = t1["elapsed"] / t2["elapsed"]
    )
  })

  results <- do.call(rbind, results)
  print(results)

  # Plot scaling
  p <- ggplot(results, aes(x = nCells)) +
    geom_line(aes(y = OriginalTime, color = "Original")) +
    geom_line(aes(y = OptimizedTime, color = "Optimized")) +
    geom_point(aes(y = OriginalTime, color = "Original")) +
    geom_point(aes(y = OptimizedTime, color = "Optimized")) +
    scale_y_log10() +
    labs(title = "Scaling Behavior",
         x = "Number of Cells",
         y = "Time (seconds, log scale)",
         color = "Method") +
    theme_minimal()

  print(p)
  ggsave("scaling_comparison.pdf", p, width = 8, height = 6)

  return(results)
}

#' Main benchmark execution
main <- function() {
  cat("Celda Clustering Optimization Benchmarks\n")
  cat("========================================\n")

  # Generate test data
  cat("\nGenerating test data...\n")
  sim <- generateTestData(nGenes = 1000, nCells = 1000, K = 8, L = 15)

  # Run benchmarks
  splitResults <- benchmarkSplitOptimizations(sim)
  gibbsResults <- benchmarkBatchGibbs(sim)
  chainResults <- benchmarkChainManagement(sim)
  scalingResults <- testScaling()

  # Save results
  allResults <- list(
    split = splitResults,
    gibbs = gibbsResults,
    chains = chainResults,
    scaling = scalingResults
  )

  saveRDS(allResults, "clustering_optimization_benchmarks.rds")

  cat("\n=== Summary ===\n")
  cat(sprintf("Split optimizations: %.1fx speedup\n",
              mean(splitResults$Speedup[-1])))
  cat(sprintf("Batch Gibbs: %.1fx speedup over standard Gibbs\n",
              gibbsResults$Speedup[2]))
  cat(sprintf("Chain management: %.1fx speedup\n",
              chainResults$Speedup[2]))
  cat(sprintf("Overall combined speedup: %.1f-%.1fx\n",
              min(scalingResults$Speedup),
              max(scalingResults$Speedup)))
}

if (!interactive()) {
  main()
}
```

#### 6.2 Quality Validation Tests

**New File**: `tests/testthat/test-clustering-optimizations.R`

```r
context("Clustering Optimizations")

test_that("Adaptive split produces equivalent results", {
  data(celdaCMod)

  # Run with and without adaptive splits
  set.seed(12345)
  result1 <- celda_C(celdaCMod$counts, K = 5,
                     splitAdaptive = FALSE, nchains = 1)

  set.seed(12345)
  result2 <- celda_C(celdaCMod$counts, K = 5,
                     splitAdaptive = TRUE, nchains = 1)

  # Log-likelihoods should be very similar
  ll1 <- metadata(result1)$celda_parameters$finalLogLik
  ll2 <- metadata(result2)$celda_parameters$finalLogLik
  expect_lt(abs((ll1 - ll2) / ll1), 0.01)  # Within 1%

  # Clustering should be similar (Adjusted Rand Index)
  z1 <- colData(result1)$celda_cell_cluster
  z2 <- colData(result2)$celda_cell_cluster
  ari <- mclust::adjustedRandIndex(z1, z2)
  expect_gt(ari, 0.9)  # High agreement
})

test_that("Batch Gibbs produces equivalent results to standard Gibbs", {
  data(celdaCMod)

  set.seed(12345)
  result1 <- celda_C(celdaCMod$counts, K = 5,
                     algorithm = "Gibbs",
                     splitOnIter = -1,
                     maxIter = 50,
                     nchains = 1)

  set.seed(12345)
  result2 <- celda_C(celdaCMod$counts, K = 5,
                     algorithm = "GibbsBatch",
                     splitOnIter = -1,
                     maxIter = 50,
                     nchains = 1)

  # Check log-likelihood convergence
  ll1 <- metadata(result1)$celda_parameters$finalLogLik
  ll2 <- metadata(result2)$celda_parameters$finalLogLik
  expect_lt(abs((ll1 - ll2) / ll1), 0.02)  # Within 2%
})

test_that("Parallel split evaluation produces deterministic results", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 2)

  data(celdaCMod)

  # Run twice with same seed
  set.seed(12345)
  result1 <- celda_C(celdaCMod$counts, K = 5,
                     nCores = 2, nchains = 1)

  set.seed(12345)
  result2 <- celda_C(celdaCMod$counts, K = 5,
                     nCores = 2, nchains = 1)

  z1 <- colData(result1)$celda_cell_cluster
  z2 <- colData(result2)$celda_cell_cluster

  # Results should be identical
  expect_identical(z1, z2)
})

test_that("Cluster quality diagnostics run without error", {
  data(celdaCMod)
  result <- celda_C(celdaCMod$counts, K = 5, nchains = 1)

  # Calculate diagnostics
  diag <- expect_no_error(celdaClusterQuality(result))

  # Check structure
  expect_s3_class(diag, "celdaDiagnostics")
  expect_true("silhouette" %in% names(diag))
  expect_true("separation" %in% names(diag))
  expect_true("wcss" %in% names(diag))

  # Check validity of metrics
  meanSil <- mean(diag$silhouette[, 3])
  expect_true(meanSil >= -1 && meanSil <= 1)
  expect_true(diag$separation$min > 0)
})
```

---

## Implementation Timeline

### Week 1-2: Phase 1 - Adaptive Split Heuristic
- **Tasks**:
  - Implement adaptive split frequency
  - Add split candidate pre-filtering
  - Add split result caching
  - Unit tests
- **Deliverable**: Faster split evaluation with maintained quality

### Week 2-3: Phase 2 - Parallel Split Evaluation
- **Tasks**:
  - Add parallel infrastructure to `.cCSplitZ()`
  - Propagate `nCores` parameter through API
  - Handle edge cases (single core, small K)
  - Unit tests
- **Deliverable**: Near-linear speedup with multiple cores

### Week 3-5: Phase 3 - Batch Gibbs Sampling
- **Tasks**:
  - Implement `.cCCalcGibbsProbZ_Batch()`
  - Add to algorithm dispatcher
  - Tune default batch size
  - Extensive validation tests
- **Deliverable**: 1.5-2x speedup over standard Gibbs

### Week 5-6: Phase 4 - Convergence Diagnostics
- **Tasks**:
  - Implement silhouette calculations
  - Add cluster separation metrics
  - Create visualization functions
  - Documentation
- **Deliverable**: Quality assessment tools for users

### Week 6-7: Phase 5 - Smart Chain Management
- **Tasks**:
  - Implement early chain termination
  - Add consensus clustering
  - Store per-chain results
  - Update metadata structure
- **Deliverable**: 20-40% speedup in multi-chain runs

### Week 7-8: Phase 6 - Benchmarking and Validation
- **Tasks**:
  - Create comprehensive benchmark suite
  - Run scaling tests
  - Validate quality preservation
  - Write performance report
- **Deliverable**: Validated performance improvements

### Week 8: Documentation and Release Preparation
- **Tasks**:
  - Update all roxygen documentation
  - Update vignettes with new features
  - Update NEWS.md
  - Update BENCHMARKING.md
  - Create migration guide for users
- **Deliverable**: Release-ready documentation

---

## Testing Strategy

### Unit Tests
- Each new function has dedicated tests in `tests/testthat/`
- Test both correctness and edge cases
- Use small synthetic datasets for speed

### Integration Tests
- Test full workflow with optimizations enabled
- Compare results to baseline (without optimizations)
- Validate log-likelihood convergence

### Performance Tests
- Benchmark on datasets of varying sizes
- Test scaling behavior (cells: 500 → 10,000)
- Profile memory usage

### Quality Tests
- Compare clustering quality metrics (ARI, Silhouette)
- Ensure optimizations don't degrade results
- Test on real single-cell datasets

---

## Expected Performance Improvements

| Optimization | Expected Speedup | Conditions |
|-------------|------------------|------------|
| Adaptive splits | 1.3-1.5x | Typical workflows |
| Parallel splits | 2-3x | With 4+ cores, K ≥ 8 |
| Batch Gibbs | 1.5-2x | Large datasets (>5000 cells) |
| Chain management | 1.2-1.4x | Multi-chain runs (nchains ≥ 3) |
| **Combined** | **3-5x** | All optimizations enabled |

### Scaling Improvements
- **Current**: O(K × M × iterations)
- **Optimized**: Effective O(K/p × M/b × iterations × 0.7) where:
  - p = parallelization factor
  - b = batch size factor
  - 0.7 = adaptive split reduction

---

## Backward Compatibility

### API Changes
All new parameters have sensible defaults to maintain backward compatibility:
- `splitAdaptive = TRUE` (can set to FALSE for old behavior)
- `nCores = 1` (sequential by default)
- `algorithm = c("EM", "Gibbs", "GibbsBatch")` (GibbsBatch is new)
- `earlyChainStop = TRUE` (can disable)

### Migration Path
Users can:
1. Continue using existing code without changes
2. Gradually adopt new features by setting parameters
3. Use `celdaClusterQuality()` to assess improvements

---

## Success Metrics

1. **Performance**: 3-5x overall speedup on representative datasets
2. **Quality**: ARI > 0.95 compared to baseline results
3. **Scalability**: Linear scaling with number of cores
4. **Usability**: Positive user feedback on diagnostics
5. **Stability**: All existing tests pass
6. **Coverage**: >90% code coverage for new functions

---

## Risk Mitigation

### Risk: Performance improvements don't materialize
- **Mitigation**: Benchmark early and often; pivot if needed

### Risk: Quality degradation
- **Mitigation**: Extensive validation tests; abort if ARI < 0.9

### Risk: Complexity increases maintenance burden
- **Mitigation**: Comprehensive documentation; modular design

### Risk: Breaking changes upset users
- **Mitigation**: Maintain backward compatibility; provide migration guide

---

## Future Enhancements (Post-Implementation)

1. **Variational Bayes inference** (celda_G)
2. **Multi-resolution clustering** (coarse-to-fine)
3. **GPU acceleration** for matrix operations
4. **Online learning** for streaming data
5. **Hierarchical clustering** within celda framework
6. **Automatic parameter selection** (K, L tuning)

---

## References and Resources

### Internal Documentation
- `BENCHMARKING.md`: Current performance characteristics
- `MODULE_DECISION_TREE_README.md`: Recursive splitting approach
- `CLAUDE.md`: Development guidelines

### External Resources
- Collapsed Gibbs sampling: Griffiths & Steyvers (2004)
- Batch updates: Newman et al. (2009)
- Parallel MCMC: Wilkinson (2006)
- Cluster validation: Rousseeuw (1987) - Silhouette method

---

## Contact and Questions

For questions about this plan:
- Open an issue at: https://github.com/campbio/celda/issues
- Tag with: `enhancement`, `performance`, `optimization`

---

**Document Version**: 1.0
**Last Updated**: 2025-11-13
**Next Review**: After Phase 2 completion
