.initializeCluster <- function(N,
                               len,
                               z = NULL,
                               initial = NULL,
                               fixed = NULL) {

  # If initial values are given, then they will not be randomly initialized
  if (!is.null(initial)) {
    initValues <- sort(unique(initial))
    if (length(unique(initial)) != N || length(initial) != len ||
      !all(initValues %in% seq(N))) {
      stop(
        "'initial' needs to be a vector of length 'len'",
        " containing N unique values."
      )
    }
    z <- as.integer(as.factor(initial))
  } else {
    z <- rep(NA, len)
  }

  # Set any values that need to be fixed during sampling
  if (!is.null(fixed)) {
    fixedValues <- sort(unique(fixed))
    if (length(fixed) != len || !all(fixedValues %in% seq(N))) {
      stop(
        "'fixed' to be a vector of length 'len' where each entry is",
        " one of N unique values or NA."
      )
    }
    fixedIx <- !is.na(fixed)
    z[fixedIx] <- fixed[fixedIx]
    zNotUsed <- setdiff(seq(N), unique(fixed[fixedIx]))
  } else {
    zNotUsed <- seq(N)
    fixedIx <- rep(FALSE, len)
  }

  # Randomly sample remaining values
  zNa <- which(is.na(z))
  if (length(zNa) > 0) {
    z[zNa] <- sample(zNotUsed, length(zNa), replace = TRUE)
  }

  # Check to ensure each value is in the vector at least once
  missing <- setdiff(seq(N), z)
  for (i in missing) {
    ta <- sort(table(z[!fixedIx]), decreasing = TRUE)
    if (ta[1] == 1) {
      stop("'len' is not long enough to accomodate 'N' unique values")
    }
    ix <- which(z == as.integer(names(ta))[1] & !fixedIx)
    z[sample(ix, 1)] <- i
  }

  return(z)
}


# Adaptive subcluster selection for cell clustering initialization
# Uses silhouette scores from quick hierarchical clustering to determine
# optimal number of subclusters for split initialization
.adaptiveKSubcluster <- function(counts, K, samplingSize = 1000) {
  nCells <- ncol(counts)

  # For very small K, just use the fixed heuristic
  if (K < 4) {
    return(ceiling(sqrt(K)))
  }

  # Handle small datasets - can't subsample
  if (nCells <= samplingSize || nCells < 100) {
    return(ceiling(sqrt(K)))
  }

  # Sample cells for quick assessment
  sampleIdx <- sample(seq_len(nCells), min(samplingSize, nCells))
  countsSample <- counts[, sampleIdx, drop = FALSE]

  # Filter to most variable genes for faster distance calculation
  nGenesUse <- min(500, nrow(countsSample))
  if (nrow(countsSample) > nGenesUse) {
    geneVars <- apply(countsSample, 1, var)
    topGenes <- order(geneVars, decreasing = TRUE)[seq_len(nGenesUse)]
    countsSample <- countsSample[topGenes, , drop = FALSE]
  }

  # Quick hierarchical clustering with sqrt(K) clusters
  kTest <- ceiling(sqrt(K))

  # Transpose for cell-cell distance
  distMatrix <- tryCatch({
    dist(t(as.matrix(countsSample)))
  }, error = function(e) {
    # If distance calculation fails, fall back to default
    return(NULL)
  })

  if (is.null(distMatrix)) {
    return(kTest)
  }

  # Use fastcluster if available, else stats::hclust
  hcResult <- tryCatch({
    if (requireNamespace("fastcluster", quietly = TRUE)) {
      fastcluster::hclust(distMatrix, method = "average")
    } else {
      stats::hclust(distMatrix, method = "average")
    }
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(hcResult)) {
    return(kTest)
  }

  # Cut tree to get sqrt(K) clusters
  clusters <- cutree(hcResult, k = kTest)

  # Calculate average silhouette score
  avgSilhouette <- tryCatch({
    if (requireNamespace("cluster", quietly = TRUE)) {
      silResult <- cluster::silhouette(clusters, distMatrix)
      mean(silResult[, "sil_width"])
    } else {
      # If cluster package not available, use default
      NA
    }
  }, error = function(e) {
    NA
  })

  # Adaptive logic based on silhouette score
  # High silhouette (>0.5) = clear structure, use fewer subclusters
  # Low silhouette (<0.2) = diffuse structure, use more subclusters
  # Medium silhouette = use default sqrt(K)
  if (!is.na(avgSilhouette)) {
    if (avgSilhouette > 0.5) {
      # Clear structure - fewer subclusters needed
      result <- ceiling(sqrt(K) * 0.7)
    } else if (avgSilhouette < 0.2) {
      # Diffuse structure - more subclusters helpful
      result <- ceiling(sqrt(K) * 1.3)
    } else {
      # Medium structure - use default
      result <- kTest
    }
  } else {
    result <- kTest
  }

  # Bound result to reasonable range
  result <- max(2, min(K, result))

  return(as.integer(result))
}


# Adaptive subcluster selection for gene module initialization
# Similar to .adaptiveKSubcluster but operates on genes instead of cells
.adaptiveLSubcluster <- function(counts, L, samplingSize = 100) {
  nGenes <- nrow(counts)

  # For very small L, just use the fixed heuristic
  if (L < 4) {
    return(ceiling(sqrt(L)))
  }

  # Handle small datasets
  if (nGenes <= samplingSize || nGenes < 50) {
    return(ceiling(sqrt(L)))
  }

  # Sample genes for quick assessment
  sampleIdx <- sample(seq_len(nGenes), min(samplingSize, nGenes))
  countsSample <- counts[sampleIdx, , drop = FALSE]

  # Filter to cells with non-zero expression for sampled genes
  cellSums <- colSums(countsSample)
  if (sum(cellSums > 0) < 10) {
    # Too few informative cells, fall back to default
    return(ceiling(sqrt(L)))
  }
  countsSample <- countsSample[, cellSums > 0, drop = FALSE]

  # Further subsample cells if too many
  if (ncol(countsSample) > 500) {
    cellIdx <- sample(seq_len(ncol(countsSample)), 500)
    countsSample <- countsSample[, cellIdx, drop = FALSE]
  }

  # Quick hierarchical clustering with sqrt(L) clusters
  lTest <- ceiling(sqrt(L))

  # Gene-gene distance
  distMatrix <- tryCatch({
    dist(as.matrix(countsSample))
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(distMatrix)) {
    return(lTest)
  }

  # Use fastcluster if available
  hcResult <- tryCatch({
    if (requireNamespace("fastcluster", quietly = TRUE)) {
      fastcluster::hclust(distMatrix, method = "average")
    } else {
      stats::hclust(distMatrix, method = "average")
    }
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(hcResult)) {
    return(lTest)
  }

  # Cut tree to get sqrt(L) clusters
  clusters <- cutree(hcResult, k = lTest)

  # Calculate average silhouette score
  avgSilhouette <- tryCatch({
    if (requireNamespace("cluster", quietly = TRUE)) {
      silResult <- cluster::silhouette(clusters, distMatrix)
      mean(silResult[, "sil_width"])
    } else {
      NA
    }
  }, error = function(e) {
    NA
  })

  # Adaptive logic
  if (!is.na(avgSilhouette)) {
    if (avgSilhouette > 0.5) {
      result <- ceiling(sqrt(L) * 0.7)
    } else if (avgSilhouette < 0.2) {
      result <- ceiling(sqrt(L) * 1.3)
    } else {
      result <- lTest
    }
  } else {
    result <- lTest
  }

  # Bound result
  result <- max(2, min(L, result))

  return(as.integer(result))
}


.initializeSplitZ <- function(counts,
                              K,
                              KSubcluster = NULL,
                              alpha = 1,
                              beta = 1,
                              minCell = 3,
                              adaptiveSubclusters = FALSE,
                              markerGenes = NULL,
                              priorClustering = NULL) {

  # If marker genes provided, use marker-guided initialization
  if (!is.null(markerGenes)) {
    return(.initializeSplitZ_MarkerGuided(
      counts = counts,
      K = K,
      markerGenes = markerGenes,
      minCell = minCell
    ))
  }

  # If prior clustering provided, refine it
  if (!is.null(priorClustering)) {
    return(.initializeSplitZ_PriorClustering(
      counts = counts,
      K = K,
      priorClustering = priorClustering,
      minCell = minCell
    ))
  }

  # Otherwise, use existing split initialization method
  s <- rep(1, ncol(counts))

  # Use adaptive subcluster selection if enabled and KSubcluster not specified
  if (is.null(KSubcluster)) {
    if (adaptiveSubclusters) {
      KSubcluster <- .adaptiveKSubcluster(counts, K)
    } else {
      KSubcluster <- ceiling(sqrt(K))
    }
  }

  # Initialize the model with KSubcluster clusters
  res <- .celda_C(
    counts,
    K = min(KSubcluster, ncol(counts)),
    maxIter = 20,
    zInitialize = "random",
    alpha = alpha,
    beta = beta,
    splitOnIter = -1,
    splitOnLast = FALSE,
    verbose = FALSE,
    reorder = FALSE
  )
  overallZ <- as.integer(as.factor(celdaClusters(res)$z))
  currentK <- max(overallZ)

  counter <- 0
  while (currentK < K & counter < 25) {
    # Determine which clusters are split-able
    # KRemaining <- K - currentK
    KPerCluster <- min(ceiling(K / currentK), KSubcluster)
    KToUse <- ifelse(KPerCluster < 2, 2, KPerCluster)

    zTa <- tabulate(overallZ, max(overallZ))

    zToSplit <- which(zTa > minCell & zTa > KToUse)
    if (length(zToSplit) > 1) {
      zToSplit <- sample(zToSplit)
    } else if (length(zToSplit) == 0) {
      break
    }

    # Cycle through each splitable cluster and split it up into
    # K.sublcusters
    for (i in zToSplit) {
      clustLabel <- .celda_C(counts[, overallZ == i, drop = FALSE],
        K = KToUse,
        zInitialize = "random",
        alpha = alpha,
        beta = beta,
        maxIter = 20,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE,
        reorder = FALSE)
      tempZ <- as.integer(as.factor(celdaClusters(clustLabel)$z))

      # Reassign clusters with label > 1
      splitIx <- tempZ > 1
      ix <- overallZ == i
      newZ <- overallZ[ix]
      newZ[splitIx] <- currentK + tempZ[splitIx] - 1

      overallZ[ix] <- newZ
      currentK <- max(overallZ)

      # Ensure that the maximum number of clusters does not get too large'
      if (currentK > K + 10) {
        break
      }
    }
    counter <- counter + 1
  }

  # Decompose counts for likelihood calculation
  p <- .cCDecomposeCounts(counts, s, overallZ, currentK)
  nS <- p$nS
  nG <- p$nG
  nM <- p$nM
  mCPByS <- p$mCPByS
  nGByCP <- p$nGByCP
  nCP <- p$nCP
  nByC <- p$nByC

  # Remove clusters 1-by-1 until K is reached
  while (currentK > K) {
    # Find second best assignment give current assignments for each cell
    probs <- .cCCalcEMProbZ(counts,
      s = s,
      z = overallZ,
      K = currentK,
      mCPByS = mCPByS,
      nGByCP = nGByCP,
      nByC = nByC,
      nCP = nCP,
      nG = nG,
      nM = nM,
      alpha = alpha,
      beta = beta,
      doSample = FALSE)
    zProb <- t(as.matrix(probs$probs))
    zProb[cbind(seq(nrow(zProb)), overallZ)] <- NA
    zSecond <- apply(zProb, 1, which.max)

    zTa <- tabulate(overallZ, currentK)
    zNonEmpty <- which(zTa > 0)

    # Find worst cluster by logLik to remove
    previousZ <- overallZ
    llShuffle <- rep(NA, currentK)
    for (i in zNonEmpty) {
      ix <- overallZ == i
      newZ <- overallZ
      newZ[ix] <- zSecond[ix]

      p <- .cCReDecomposeCounts(
        counts,
        s,
        newZ,
        previousZ,
        nGByCP,
        currentK)
      nGByCP <- p$nGByCP
      mCPByS <- p$mCPByS
      llShuffle[i] <- .cCCalcLL(
        mCPByS,
        nGByCP,
        s,
        newZ,
        currentK,
        nS,
        nG,
        alpha,
        beta)
      previousZ <- newZ
    }

    # Remove the cluster which had the the largest likelihood after removal
    zToRemove <- which.max(llShuffle)

    ix <- overallZ == zToRemove
    overallZ[ix] <- zSecond[ix]

    p <- .cCReDecomposeCounts(counts,
      s,
      overallZ,
      previousZ,
      nGByCP,
      currentK)
    nGByCP <- p$nGByCP[, -zToRemove, drop = FALSE]
    mCPByS <- p$mCPByS[-zToRemove, , drop = FALSE]
    overallZ <- as.integer(as.factor(overallZ))
    currentK <- currentK - 1
  }
  return(overallZ)
}


.initializeSplitY <- function(counts,
                              L,
                              LSubcluster = NULL,
                              tempK = 100,
                              beta = 1,
                              delta = 1,
                              gamma = 1,
                              minFeature = 3,
                              adaptiveSubclusters = FALSE) {

  # Use adaptive subcluster selection if enabled and LSubcluster not specified
  if (is.null(LSubcluster)) {
    if (adaptiveSubclusters) {
      LSubcluster <- .adaptiveLSubcluster(counts, L)
    } else {
      LSubcluster <- ceiling(sqrt(L))
    }
  }

  # Collapse cells to managable number of clusters
  if (!is.null(tempK) && ncol(counts) > tempK) {
    z <- .initializeSplitZ(counts, K = tempK)
    counts <- .colSumByGroup(counts, z, length(unique(z)))
  }

  # Initialize the model with KSubcluster clusters
  res <- .celda_G(counts,
    L = LSubcluster,
    maxIter = 10,
    yInitialize = "random",
    beta = beta,
    delta = delta,
    gamma = gamma,
    splitOnIter = -1,
    splitOnLast = FALSE,
    verbose = FALSE,
    reorder = FALSE)
  overallY <- as.integer(as.factor(celdaClusters(res)$y))
  currentL <- max(overallY)

  counter <- 0
  while (currentL < L & counter < 25) {
    # Determine which clusters are split-able
    yTa <- tabulate(overallY, max(overallY))
    yToSplit <- sample(which(yTa > minFeature & yTa > LSubcluster))

    if (length(yToSplit) == 0) {
      break
    }

    # Cycle through each splitable cluster and split it up into
    # LSublcusters
    for (i in yToSplit) {
      # make sure the colSums of subset counts is not 0
      countsY <- counts[overallY == i, , drop = FALSE]
      countsY <- countsY[, !(colSums(countsY) == 0)]

      if (ncol(countsY) == 0) {
        next
      }

      clustLabel <- .celda_G(
        countsY,
        L = min(LSubcluster, nrow(countsY)),
        yInitialize = "random",
        beta = beta,
        delta = delta,
        gamma = gamma,
        maxIter = 20,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE,
        reorder = FALSE
      )
      tempY <- as.integer(as.factor(celdaClusters(clustLabel)$y))

      # Reassign clusters with label > 1
      splitIx <- tempY > 1
      ix <- overallY == i
      newY <- overallY[ix]
      newY[splitIx] <- currentL + tempY[splitIx] - 1

      overallY[ix] <- newY
      currentL <- max(overallY)

      # Ensure that the maximum number of clusters does not get too large
      if (currentL > L + 10) {
        break
      }
    }
    counter <- counter + 1
  }

  ## Decompose counts for likelihood calculation
  p <- .cGDecomposeCounts(counts = counts, y = overallY, L = currentL)
  nTSByC <- p$nTSByC
  nByG <- p$nByG
  nByTS <- p$nByTS
  nGByTS <- p$nGByTS
  nM <- p$nM
  nG <- p$nG
  rm(p)

  # Pre-compute lgamma values
  lgbeta <- lgamma((seq(0, max(colSums(counts)))) + beta)
  lggamma <- lgamma(seq(0, nrow(counts) + L) + gamma)
  lgdelta <- c(NA, lgamma(seq(nrow(counts) + L) * delta))

  # Remove clusters 1-by-1 until L is reached
  while (currentL > L) {
    # Find second best assignment give current assignments for each cell
    probs <- .cGCalcGibbsProbY(
      counts = counts,
      y = overallY,
      L = currentL,
      nTSByC = nTSByC,
      nByTS = nByTS,
      nGByTS = nGByTS,
      nByG = nByG,
      nG = nG,
      beta = beta,
      delta = delta,
      gamma = gamma,
      lgbeta = lgbeta,
      lggamma = lggamma,
      lgdelta = lgdelta,
      doSample = FALSE)
    yProb <- t(probs$probs)
    yProb[cbind(seq(nrow(yProb)), overallY)] <- NA
    ySecond <- apply(yProb, 1, which.max)

    yTa <- tabulate(overallY, currentL)
    yNonEmpty <- which(yTa > 0)

    # Find worst cluster by logLik to remove
    previousY <- overallY
    llShuffle <- rep(NA, currentL)
    for (i in yNonEmpty) {
      ix <- overallY == i
      newY <- overallY
      newY[ix] <- ySecond[ix]

      # Move arounds counts for likelihood calculation
      p <- .cGReDecomposeCounts(
        counts,
        newY,
        previousY,
        nTSByC,
        nByG,
        currentL)
      nTSByC <- p$nTSByC
      nGByTS <- p$nGByTS
      nByTS <- p$nByTS
      llShuffle[i] <- .cGCalcLL(
        nTSByC,
        nByTS,
        nByG,
        nGByTS,
        nM,
        nG,
        currentL,
        beta,
        delta,
        gamma)
      previousY <- newY
    }

    # Remove the cluster which had the the largest likelihood after removal
    yToRemove <- which.max(llShuffle)

    ix <- overallY == yToRemove
    overallY[ix] <- ySecond[ix]

    # Move around counts and remove module
    p <- .cGReDecomposeCounts(
      counts,
      overallY,
      previousY,
      nTSByC,
      nByG,
      currentL)
    nTSByC <- p$nTSByC[-yToRemove, , drop = FALSE]
    nGByTS <- p$nGByTS[-yToRemove]
    nByTS <- p$nByTS[-yToRemove]
    overallY <- as.integer(as.factor(overallY))
    currentL <- currentL - 1
  }
  return(overallY)
}
