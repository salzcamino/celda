#' @title Plot consensus clustering confidence scores
#' @description Visualizes per-cell or per-gene confidence scores from consensus
#'   clustering, highlighting low-confidence assignments.
#' @param sce A \linkS4class{SingleCellExperiment} object containing consensus
#'   clustering results
#' @param type Character. Either "cells" or "genes" to plot cell or gene
#'   confidence scores
#' @param reducedDimName Character. Name of reduced dimension to use for plotting
#'   cells (e.g., "celda_UMAP", "UMAP", "TSNE"). Only used when type = "cells".
#'   Default "celda_UMAP".
#' @param minConfidence Numeric. Threshold for flagging low-confidence points.
#'   Default 0.7.
#' @param size Numeric. Point size for plotting. Default 0.5.
#' @param alpha Numeric. Point transparency. Default 0.8.
#' @return A ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' # Assuming sce has consensus clustering results
#' plotConsensusConfidence(sce, type = "cells")
#' plotConsensusConfidence(sce, type = "genes")
#' }
plotConsensusConfidence <- function(sce,
                                     type = c("cells", "genes"),
                                     reducedDimName = "celda_UMAP",
                                     minConfidence = 0.7,
                                     size = 0.5,
                                     alpha = 0.8) {
    type <- match.arg(type)

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. ",
             "Please install it with install.packages('ggplot2')")
    }

    if (type == "cells") {
        # Check for confidence scores in colData
        if (!"celda_z_confidence" %in% colnames(SummarizedExperiment::colData(sce))) {
            stop("No cell confidence scores found. ",
                 "Run celda_CG with useConsensus=TRUE to generate them.")
        }

        # Get confidence scores
        confidence <- SummarizedExperiment::colData(sce)$celda_z_confidence

        # Get reduced dimensions
        if (!reducedDimName %in% SingleCellExperiment::reducedDimNames(sce)) {
            stop("Reduced dimension '", reducedDimName, "' not found in sce. ",
                 "Available: ", paste(SingleCellExperiment::reducedDimNames(sce), collapse = ", "))
        }

        coords <- SingleCellExperiment::reducedDim(sce, reducedDimName)
        if (ncol(coords) < 2) {
            stop("Reduced dimension must have at least 2 dimensions")
        }

        # Create plot data
        plotData <- data.frame(
            x = coords[, 1],
            y = coords[, 2],
            confidence = confidence,
            lowConfidence = confidence < minConfidence
        )

        # Create plot
        p <- ggplot2::ggplot(plotData, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_point(ggplot2::aes(color = confidence),
                                size = size, alpha = alpha) +
            ggplot2::scale_color_gradient2(
                low = "red", mid = "yellow", high = "blue",
                midpoint = minConfidence,
                name = "Consensus\nConfidence"
            ) +
            ggplot2::labs(
                title = "Cell Consensus Clustering Confidence",
                subtitle = paste0(
                    sum(plotData$lowConfidence),
                    " cells with confidence < ",
                    minConfidence
                ),
                x = paste0(reducedDimName, "_1"),
                y = paste0(reducedDimName, "_2")
            ) +
            ggplot2::theme_minimal()

    } else if (type == "genes") {
        # Check for confidence scores in rowData
        if (!"celda_y_confidence" %in% colnames(SummarizedExperiment::rowData(sce))) {
            stop("No gene confidence scores found. ",
                 "Run celda_CG with useConsensus=TRUE to generate them.")
        }

        # Get confidence scores
        confidence <- SummarizedExperiment::rowData(sce)$celda_y_confidence

        # Create histogram
        plotData <- data.frame(
            confidence = confidence,
            lowConfidence = confidence < minConfidence
        )

        p <- ggplot2::ggplot(plotData, ggplot2::aes(x = confidence)) +
            ggplot2::geom_histogram(
                ggplot2::aes(fill = lowConfidence),
                bins = 50,
                color = "black",
                alpha = 0.7
            ) +
            ggplot2::scale_fill_manual(
                values = c("TRUE" = "red", "FALSE" = "blue"),
                name = "Low Confidence",
                labels = c("FALSE" = paste0(">= ", minConfidence),
                          "TRUE" = paste0("< ", minConfidence))
            ) +
            ggplot2::labs(
                title = "Gene Module Consensus Confidence",
                subtitle = paste0(
                    sum(plotData$lowConfidence),
                    " genes with confidence < ",
                    minConfidence
                ),
                x = "Consensus Confidence",
                y = "Number of Genes"
            ) +
            ggplot2::theme_minimal()
    }

    return(p)
}


#' @title Plot chain agreement heatmap
#' @description Creates a heatmap showing how often cells or genes cluster
#'   together across multiple chains.
#' @param chainResults List of chain results from celda clustering
#' @param type Character. Either "z" for cells or "y" for genes
#' @param maxElements Integer. Maximum number of elements to include in heatmap.
#'   For large datasets, a random sample will be taken. Default 500.
#' @param clusterColors Logical. Whether to color rows/columns by cluster
#'   assignment. Default TRUE.
#' @return A plot object (using pheatmap or similar)
#' @keywords internal
.plotChainAgreement <- function(chainResults,
                                 type = c("z", "y"),
                                 maxElements = 500,
                                 clusterColors = TRUE) {
    type <- match.arg(type)

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required for this plot. ",
             "Please install it with install.packages('pheatmap')")
    }

    # Build co-occurrence matrix
    cooccur <- .buildCooccurrenceMatrix(chainResults, type = type)

    # Sample if too large
    n <- nrow(cooccur)
    if (n > maxElements) {
        idx <- sample(seq_len(n), maxElements)
        cooccur <- cooccur[idx, idx]

        # Get cluster assignments for annotation
        assignments <- chainResults[[1]][[type]][idx]
    } else {
        assignments <- chainResults[[1]][[type]]
    }

    # Create annotation
    if (clusterColors) {
        annotation <- data.frame(
            Cluster = as.factor(assignments)
        )
        rownames(annotation) <- rownames(cooccur)
    } else {
        annotation <- NULL
    }

    # Create heatmap
    p <- pheatmap::pheatmap(
        cooccur,
        color = grDevices::colorRampPalette(c("white", "blue"))(100),
        annotation_row = annotation,
        annotation_col = annotation,
        show_rownames = FALSE,
        show_colnames = FALSE,
        main = paste0(
            "Co-clustering Frequency Across Chains (",
            ifelse(type == "z", "Cells", "Genes"),
            ")"
        ),
        legend = TRUE
    )

    return(p)
}
