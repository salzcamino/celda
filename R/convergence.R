#' Advanced Convergence Detection for Celda Models
#'
#' This file contains functions for detecting convergence in celda clustering
#' algorithms using both log-likelihood and cluster stability criteria.


#' @title Calculate Adjusted Rand Index
#' @description Calculate the Adjusted Rand Index (ARI) between two clusterings.
#'   ARI measures similarity between cluster assignments, adjusted for chance.
#'   Returns 1 for identical clusterings, 0 for random clusterings, and
#'   negative values for worse than random agreement.
#' @param x Integer vector of cluster assignments
#' @param y Integer vector of cluster assignments
#' @return Numeric ARI value
#' @keywords internal
.calculateARI <- function(x, y) {
    # Handle edge cases
    if (length(x) != length(y)) {
        stop("x and y must have the same length")
    }
    if (length(x) == 0) {
        return(NA_real_)
    }

    # If all same cluster, ARI is 1
    if (length(unique(x)) == 1 && length(unique(y)) == 1) {
        return(1)
    }

    # Try to use mclust's implementation if available
    if (requireNamespace("mclust", quietly = TRUE)) {
        return(mclust::adjustedRandIndex(x, y))
    }

    # Otherwise, implement ARI calculation
    # Create contingency table
    tab <- table(x, y)

    # Calculate choose(n, 2) for each cell
    n <- sum(tab)

    # Sum of choose(n_ij, 2) over all cells
    sumChoose2 <- sum(choose(tab, 2))

    # Sum of choose(a_i, 2) over rows
    rowSums <- rowSums(tab)
    sumRowChoose2 <- sum(choose(rowSums, 2))

    # Sum of choose(b_j, 2) over columns
    colSums <- colSums(tab)
    sumColChoose2 <- sum(choose(colSums, 2))

    # Expected index
    expectedIndex <- sumRowChoose2 * sumColChoose2 / choose(n, 2)

    # Max index
    maxIndex <- (sumRowChoose2 + sumColChoose2) / 2

    # ARI formula
    if (maxIndex - expectedIndex == 0) {
        return(1)
    }

    ari <- (sumChoose2 - expectedIndex) / (maxIndex - expectedIndex)
    return(ari)
}


#' @title Check Advanced Convergence
#' @description Advanced convergence detection using both log-likelihood plateau
#'   and cluster stability criteria. Provides more accurate convergence detection
#'   than simple log-likelihood checking alone.
#'
#' @param llHistory Numeric vector of log-likelihood values over iterations
#' @param zHistory List of cell cluster assignments (z vectors) from recent
#'   iterations. Only the last 5-10 iterations need to be stored.
#' @param yHistory List of gene module assignments (y vectors) from recent
#'   iterations. Can be NULL for celda_C model. Only last 5-10 iterations needed.
#' @param iter Integer. Current iteration number
#' @param stopIter Integer. Number of iterations without improvement to consider
#'   convergence. Default 10.
#' @param relTol Numeric. Relative tolerance for log-likelihood convergence.
#'   Convergence detected when (max(recent LL) - min(recent LL)) / abs(max) < relTol.
#'   Default 1e-5.
#' @param checkStability Logical. Whether to check cluster stability using ARI.
#'   Default TRUE.
#' @param ariThreshold Numeric. Minimum ARI value to consider clusters stable.
#'   Default 0.99.
#' @param minIter Integer. Minimum number of iterations before convergence can
#'   be declared. Default 10.
#'
#' @return List with components:
#'   \describe{
#'     \item{converged}{Logical. TRUE if model has converged}
#'     \item{reason}{Character. Explanation of convergence decision}
#'     \item{llConverged}{Logical. TRUE if log-likelihood has converged}
#'     \item{zStable}{Logical. TRUE if cell clusters are stable (NA if not checked)}
#'     \item{yStable}{Logical. TRUE if gene modules are stable (NA if not checked)}
#'     \item{iterations}{Integer. Current iteration number}
#'   }
#' @keywords internal
.checkConvergence_Advanced <- function(llHistory,
                                        zHistory = NULL,
                                        yHistory = NULL,
                                        iter,
                                        stopIter = 10,
                                        relTol = 1e-5,
                                        checkStability = TRUE,
                                        ariThreshold = 0.99,
                                        minIter = 10) {

    # Initialize return values
    converged <- FALSE
    reason <- ""
    llConverged <- FALSE
    zStable <- NA
    yStable <- NA

    # Minimum iterations check
    if (iter < minIter) {
        return(list(
            converged = FALSE,
            reason = paste("Iteration", iter, "< minimum required iterations",
                         minIter),
            llConverged = FALSE,
            zStable = NA,
            yStable = NA,
            iterations = iter
        ))
    }

    # Check if we have enough history
    if (length(llHistory) < stopIter + 1) {
        return(list(
            converged = FALSE,
            reason = paste("Insufficient log-likelihood history:",
                         length(llHistory), "values (need", stopIter + 1, ")"),
            llConverged = FALSE,
            zStable = NA,
            yStable = NA,
            iterations = iter
        ))
    }

    # ===== Log-likelihood convergence check =====
    # Get recent log-likelihood values
    recentLL <- tail(llHistory, stopIter + 1)
    maxLL <- max(recentLL)
    minLL <- min(recentLL)

    # Calculate relative change
    relChange <- if (abs(maxLL) > 0) {
        (maxLL - minLL) / abs(maxLL)
    } else {
        Inf
    }

    llConverged <- relChange < relTol

    # ===== Cluster stability check =====
    if (checkStability && !is.null(zHistory)) {
        # Filter out NULL entries
        zHistoryClean <- Filter(Negate(is.null), zHistory)

        if (length(zHistoryClean) >= 5) {
            # Check ARI for last 4 consecutive transitions
            # We need at least 5 z vectors to check 4 transitions
            recentZ <- tail(zHistoryClean, 5)

            ariValues <- numeric(4)
            for (i in 1:4) {
                ariValues[i] <- .calculateARI(recentZ[[i]], recentZ[[i + 1]])
            }

            # Check if all ARI values exceed threshold
            zStable <- all(ariValues > ariThreshold, na.rm = TRUE)
        } else {
            zStable <- FALSE
        }
    } else {
        zStable <- NA
    }

    # Check gene module stability if yHistory provided
    if (checkStability && !is.null(yHistory)) {
        # Filter out NULL entries
        yHistoryClean <- Filter(Negate(is.null), yHistory)

        if (length(yHistoryClean) >= 5) {
            recentY <- tail(yHistoryClean, 5)

            ariValuesY <- numeric(4)
            for (i in 1:4) {
                ariValuesY[i] <- .calculateARI(recentY[[i]], recentY[[i + 1]])
            }

            yStable <- all(ariValuesY > ariThreshold, na.rm = TRUE)
        } else {
            yStable <- FALSE
        }
    } else {
        yStable <- NA
    }

    # ===== Combined convergence decision =====

    # Case 1: Both LL and clusters stable (best case)
    if (llConverged &&
        (is.na(zStable) || isTRUE(zStable)) &&
        (is.na(yStable) || isTRUE(yStable))) {
        converged <- TRUE
        stabilityMsg <- ""
        if (!is.na(zStable) && !is.na(yStable)) {
            stabilityMsg <- " Cell clusters and gene modules are stable (ARI > 0.99)."
        } else if (!is.na(zStable)) {
            stabilityMsg <- " Cell clusters are stable (ARI > 0.99)."
        } else if (!is.na(yStable)) {
            stabilityMsg <- " Gene modules are stable (ARI > 0.99)."
        }
        reason <- paste0("Converged: Log-likelihood stable (relative change < ",
                        relTol, ").", stabilityMsg)
    }

    # Case 2: LL converged for extended period (2x stopIter), even if clusters
    # still moving slightly
    else if (length(llHistory) >= 2 * stopIter + 1) {
        extendedLL <- tail(llHistory, 2 * stopIter + 1)
        maxLLExt <- max(extendedLL)
        minLLExt <- min(extendedLL)
        relChangeExt <- if (abs(maxLLExt) > 0) {
            (maxLLExt - minLLExt) / abs(maxLLExt)
        } else {
            Inf
        }

        if (relChangeExt < relTol) {
            converged <- TRUE
            reason <- paste0("Converged: Log-likelihood stable for ",
                           2 * stopIter,
                           " iterations (relative change < ", relTol,
                           "), despite minor cluster reassignments.")
        }
    }

    # Case 3: Not converged - provide diagnostic info
    if (!converged) {
        reasons <- c()
        if (!llConverged) {
            reasons <- c(reasons,
                        paste0("LL not stable (relative change = ",
                              sprintf("%.2e", relChange), ")"))
        }
        if (isTRUE(!zStable)) {
            reasons <- c(reasons, "Cell clusters still changing")
        }
        if (isTRUE(!yStable)) {
            reasons <- c(reasons, "Gene modules still changing")
        }
        reason <- paste("Not converged:", paste(reasons, collapse = "; "))
    }

    return(list(
        converged = converged,
        reason = reason,
        llConverged = llConverged,
        zStable = zStable,
        yStable = yStable,
        iterations = iter
    ))
}


#' @title Simple Convergence Check
#' @description Simple convergence detection based only on log-likelihood
#'   improvement. This is the original celda convergence criterion.
#'
#' @param llHistory Numeric vector of log-likelihood values
#' @param currentLL Numeric. Current log-likelihood value
#' @param stopIter Integer. Number of iterations without improvement
#' @param numIterWithoutImprovement Integer. Current count of iterations
#'   without improvement
#'
#' @return List with:
#'   \describe{
#'     \item{improved}{Logical. TRUE if LL improved}
#'     \item{numIterWithoutImprovement}{Integer. Updated counter}
#'   }
#' @keywords internal
.checkConvergence_Simple <- function(llHistory,
                                      currentLL,
                                      stopIter,
                                      numIterWithoutImprovement) {

    # Check if current LL is better than all previous
    improved <- all(currentLL > llHistory)

    if (improved) {
        numIterWithoutImprovement <- 1L
    } else {
        numIterWithoutImprovement <- numIterWithoutImprovement + 1L
    }

    return(list(
        improved = improved,
        numIterWithoutImprovement = numIterWithoutImprovement
    ))
}
