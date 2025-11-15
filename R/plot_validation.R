# Visualization Functions for Cluster Validation Metrics
#
# This file contains functions for plotting cluster validation metrics
# to help with model selection and quality assessment.


#' @title Plot Cluster Validation Metrics
#' @description Visualizes cluster validation metrics across different models
#'   from celdaGridSearch results or a list of celda models.
#' @param celdaList A celdaList object from celdaGridSearch, or a list of
#'   celda model objects
#' @param metrics Character vector. Which metrics to plot. Options:
#'   "silhouette", "calinskiHarabasz", "daviesBouldin", "moduleCoherence",
#'   "all". Default c("calinskiHarabasz", "daviesBouldin").
#' @param useAssay String. Name of assay to use for calculations if metrics
#'   need to be computed. Default "counts".
#' @return A ggplot object (if ggplot2 available) or base plot
#' @keywords internal
.plotClusterMetrics <- function(celdaList,
                                metrics = c("calinskiHarabasz", "daviesBouldin"),
                                useAssay = "counts") {

    # Determine if metrics need to be calculated
    # or if they're already stored in the models
    needsCalculation <- TRUE

    # Extract or calculate metrics for each model
    modelMetrics <- list()
    modelParams <- list()

    if (is(celdaList, "celdaList")) {
        # celdaGridSearch results
        nModels <- length(runParams(celdaList)$index)

        for (i in seq_len(nModels)) {
            model <- celdaList[[i]]

            # Check if metrics already stored
            if (!is.null(model$validationMetrics)) {
                modelMetrics[[i]] <- model$validationMetrics
            } else {
                # Calculate metrics
                if (is(model, "SingleCellExperiment")) {
                    counts <- SummarizedExperiment::assay(model, useAssay)
                    z <- SummarizedExperiment::colData(model)$celda_cell_cluster
                    y <- SummarizedExperiment::rowData(model)$celda_feature_module

                    if (!is.null(z)) {
                        modelMetrics[[i]] <- .calculateClusterMetrics(
                            counts = counts,
                            z = as.integer(z),
                            y = if (!is.null(y)) as.integer(y) else NULL,
                            metrics = metrics
                        )
                    }
                }
            }

            # Extract K/L parameters
            params <- runParams(celdaList)
            modelParams[[i]] <- list(
                K = params$K[i],
                L = params$L[i]
            )
        }
    } else if (is.list(celdaList)) {
        # List of models
        for (i in seq_along(celdaList)) {
            model <- celdaList[[i]]

            # Similar logic as above
            if (!is.null(model$validationMetrics)) {
                modelMetrics[[i]] <- model$validationMetrics
            } else if (is(model, "celda_CG")) {
                # Extract from S4 object
                # This would require the counts matrix to be passed separately
                warning("Metrics not pre-calculated. Please provide counts.")
                modelMetrics[[i]] <- NULL
            }

            # Extract parameters
            if (methods::isS4(model)) {
                modelParams[[i]] <- list(
                    K = model@params$K,
                    L = model@params$L
                )
            }
        }
    }

    # Create plot data
    plotData <- data.frame()
    for (i in seq_along(modelMetrics)) {
        if (!is.null(modelMetrics[[i]])) {
            for (metricName in names(modelMetrics[[i]])) {
                plotData <- rbind(plotData, data.frame(
                    model = i,
                    K = modelParams[[i]]$K,
                    L = if (!is.null(modelParams[[i]]$L)) modelParams[[i]]$L else NA,
                    metric = metricName,
                    value = modelMetrics[[i]][[metricName]]
                ))
            }
        }
    }

    # Create plot
    if (nrow(plotData) == 0) {
        warning("No metrics to plot")
        return(NULL)
    }

    # Use ggplot2 if available
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        p <- ggplot2::ggplot(plotData,
                             ggplot2::aes(x = K, y = value, color = metric)) +
            ggplot2::geom_line() +
            ggplot2::geom_point() +
            ggplot2::facet_wrap(~ metric, scales = "free_y") +
            ggplot2::theme_bw() +
            ggplot2::labs(
                title = "Cluster Validation Metrics",
                x = "Number of Clusters (K)",
                y = "Metric Value"
            )

        # Add L to x-axis if present
        if (any(!is.na(plotData$L))) {
            p <- p + ggplot2::aes(x = interaction(K, L))
        }

        return(p)
    } else {
        # Base R plotting
        graphics::par(mfrow = c(ceiling(length(unique(plotData$metric)) / 2), 2))

        for (metricName in unique(plotData$metric)) {
            subData <- plotData[plotData$metric == metricName, ]

            graphics::plot(subData$K, subData$value,
                          type = "b",
                          xlab = "K",
                          ylab = metricName,
                          main = paste("Cluster Validation:", metricName),
                          pch = 19,
                          col = "steelblue")
        }

        graphics::par(mfrow = c(1, 1))
        return(invisible(NULL))
    }
}


#' @title Plot Model Comparison
#' @description Creates a comparison plot showing log-likelihood and validation
#'   metrics for model selection.
#' @param models List of celda model results, each containing finalLogLik
#'   and optionally validationMetrics
#' @param selectionCriterion Character. Which criterion was used for selection.
#' @param validationWeight Numeric. Weight used for validation metrics.
#' @return ggplot object or base plot
#' @keywords internal
.plotModelComparison <- function(models,
                                 selectionCriterion = "combined",
                                 validationWeight = 0.3) {

    nModels <- length(models)

    # Extract log-likelihoods
    logLiks <- sapply(models, function(x) x$finalLogLik)

    # Extract validation metrics if available
    hasMetrics <- sapply(models, function(x) !is.null(x$validationMetrics))

    plotData <- data.frame(
        model = seq_len(nModels),
        logLik = logLiks
    )

    if (any(hasMetrics)) {
        # Extract CH and DB indices
        plotData$CH <- sapply(models, function(x) {
            if (!is.null(x$validationMetrics$calinskiHarabasz)) {
                x$validationMetrics$calinskiHarabasz
            } else {
                NA
            }
        })

        plotData$DB <- sapply(models, function(x) {
            if (!is.null(x$validationMetrics$daviesBouldin)) {
                x$validationMetrics$daviesBouldin
            } else {
                NA
            }
        })
    }

    # Create plot
    if (requireNamespace("ggplot2", quietly = TRUE) &&
        requireNamespace("tidyr", quietly = TRUE)) {

        # Reshape for plotting
        plotDataLong <- tidyr::pivot_longer(plotData,
                                             cols = -model,
                                             names_to = "metric",
                                             values_to = "value")

        p <- ggplot2::ggplot(plotDataLong,
                             ggplot2::aes(x = model, y = value, color = metric)) +
            ggplot2::geom_line() +
            ggplot2::geom_point() +
            ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
            ggplot2::theme_bw() +
            ggplot2::labs(
                title = paste("Model Comparison -", selectionCriterion),
                x = "Chain",
                y = "Value"
            )

        return(p)
    } else {
        # Base R plotting
        nMetrics <- sum(!is.na(plotData[1, -1]))
        graphics::par(mfrow = c(nMetrics, 1))

        # Plot log-likelihood
        graphics::plot(plotData$model, plotData$logLik,
                      type = "b",
                      xlab = "Chain",
                      ylab = "Log-Likelihood",
                      main = "Log-Likelihood by Chain",
                      pch = 19,
                      col = "darkgreen")

        # Plot validation metrics if available
        if ("CH" %in% names(plotData) && any(!is.na(plotData$CH))) {
            graphics::plot(plotData$model, plotData$CH,
                          type = "b",
                          xlab = "Chain",
                          ylab = "Calinski-Harabasz",
                          main = "CH Index by Chain",
                          pch = 19,
                          col = "steelblue")
        }

        if ("DB" %in% names(plotData) && any(!is.na(plotData$DB))) {
            graphics::plot(plotData$model, plotData$DB,
                          type = "b",
                          xlab = "Chain",
                          ylab = "Davies-Bouldin",
                          main = "DB Index by Chain (lower is better)",
                          pch = 19,
                          col = "coral")
        }

        graphics::par(mfrow = c(1, 1))
        return(invisible(NULL))
    }
}


#' @title Summary Plot for Validation Metrics
#' @description Creates a summary visualization of all validation metrics
#'   for a single model.
#' @param validationMetrics Named list of validation metric values
#' @return Base R plot
#' @keywords internal
.plotMetricsSummary <- function(validationMetrics) {

    # Remove NA values
    validMetrics <- validationMetrics[!is.na(unlist(validationMetrics))]

    if (length(validMetrics) == 0) {
        warning("No valid metrics to plot")
        return(invisible(NULL))
    }

    # Normalize metrics to 0-1 scale for visualization
    metricValues <- unlist(validMetrics)
    metricNames <- names(validMetrics)

    # Different metrics have different scales and interpretations
    # CH: higher is better
    # DB: lower is better
    # Silhouette: -1 to 1, higher is better
    # Module coherence: -1 to 1, higher is better

    # Create bar plot
    colors <- c("calinskiHarabasz" = "steelblue",
                "daviesBouldin" = "coral",
                "silhouette" = "darkgreen",
                "moduleCoherence" = "purple",
                "clusterSizeCV" = "orange")

    barColors <- colors[metricNames]
    barColors[is.na(barColors)] <- "gray"

    graphics::barplot(metricValues,
                     names.arg = metricNames,
                     col = barColors,
                     main = "Cluster Validation Metrics Summary",
                     ylab = "Value",
                     las = 2,  # Vertical labels
                     cex.names = 0.8)

    # Add reference lines for known thresholds
    if ("silhouette" %in% metricNames) {
        graphics::abline(h = 0.5, lty = 2, col = "darkgreen")
        graphics::abline(h = 0.25, lty = 2, col = "orange")
    }

    return(invisible(NULL))
}
