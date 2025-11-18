#' @title Cross-platform parallel lapply
#' @description Internal helper function for cross-platform parallel processing
#'  using the future framework. Automatically configures the appropriate backend
#'  based on the operating system and number of cores requested.
#' @param X A vector (atomic or list) to apply the function over
#' @param FUN The function to apply to each element of X
#' @param nCores Integer. Number of cores to use. If 1, uses sequential processing.
#'  If > 1, uses parallel processing with future framework.
#' @param ... Additional arguments passed to FUN
#' @return A list of the same length as X
#' @keywords internal
#' @importFrom future plan multicore multisession sequential
#' @importFrom future.apply future_lapply
.parallelLapply <- function(X, FUN, nCores = 1, ...) {
  # If sequential processing requested or only one item, use regular lapply
  if (nCores <= 1 || length(X) <= 1) {
    return(lapply(X, FUN, ...))
  }

  # Determine appropriate future plan based on platform
  # multicore: Uses forking (Unix/macOS only, most efficient)
  # multisession: Uses separate R sessions (works on all platforms including Windows)
  oldPlan <- future::plan()
  on.exit(future::plan(oldPlan), add = TRUE)

  # Use multicore on Unix-like systems (Linux/macOS), multisession on Windows
  if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = nCores)
  } else {
    future::plan(future::multisession, workers = nCores)
  }

  # Run parallel lapply using future framework
  result <- future.apply::future_lapply(X, FUN, ...,
                                        future.seed = TRUE,
                                        future.scheduling = TRUE)

  return(result)
}


#' @title Check if future package is available
#' @description Internal helper to check if future framework is available
#' @return Logical indicating if future and future.apply are installed
#' @keywords internal
.hasFuture <- function() {
  requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)
}


#' @title Safe parallel lapply with fallback
#' @description Internal wrapper that falls back to parallel::mclapply if
#'  future is not available, and to regular lapply if neither works.
#' @param X A vector (atomic or list) to apply the function over
#' @param FUN The function to apply to each element of X
#' @param nCores Integer. Number of cores to use.
#' @param ... Additional arguments passed to FUN
#' @return A list of the same length as X
#' @keywords internal
.safeParallelLapply <- function(X, FUN, nCores = 1, ...) {
  # Try future-based approach first (cross-platform)
  if (.hasFuture()) {
    return(.parallelLapply(X, FUN, nCores = nCores, ...))
  }

  # Fallback to mclapply (Unix only, but no additional dependencies)
  if (nCores > 1 && length(X) > 1) {
    if (.Platform$OS.type == "unix") {
      return(parallel::mclapply(X, FUN, mc.cores = nCores, ...))
    } else {
      warning("Parallel processing requested but future package not available ",
              "and mclapply not supported on Windows. Using sequential processing.")
    }
  }

  # Final fallback to sequential lapply
  return(lapply(X, FUN, ...))
}
