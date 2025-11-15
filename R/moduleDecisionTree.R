#' @title Build decision tree from recursive module splitting
#' @description Constructs a hierarchical decision tree structure that traces
#'  how feature modules are split during recursive module splitting. Each node
#'  in the tree represents a module at a specific L value, with edges showing
#'  parent-child relationships when modules are split.
#' @param x Either a \linkS4class{celdaList} object or a
#'  \linkS4class{SingleCellExperiment} object returned by
#'  \link{recursiveSplitModule}. If a SingleCellExperiment is provided, the
#'  celdaList will be extracted from the metadata.
#' @param altExpName The name for the \link{altExp} slot
#'  to use if \code{x} is a SingleCellExperiment. Default "featureSubset".
#' @param minL Integer. Minimum number of modules to include in the tree.
#'  If NULL, uses the minimum L from the celdaList. Default NULL.
#' @param maxL Integer. Maximum number of modules to include in the tree.
#'  If NULL, uses the maximum L from the celdaList. Default NULL.
#' @return A list with class "moduleDecisionTree" containing:
#'  \itemize{
#'  \item \code{nodes}: Data frame with columns:
#'    \itemize{
#'      \item \code{nodeId}: Unique identifier for each node
#'      \item \code{L}: Number of modules at this level
#'      \item \code{moduleLabel}: Module label at this L
#'      \item \code{parentNode}: Parent node ID (NA for root nodes)
#'      \item \code{nGenes}: Number of genes in this module
#'      \item \code{logLik}: Log-likelihood of the model at this L
#'    }
#'  \item \code{edges}: Data frame with parent-child relationships
#'  \item \code{geneAssignments}: List mapping genes to modules at each L
#'  \item \code{celdaList}: Original celdaList object
#'  }
#' @seealso \link{recursiveSplitModule} for generating the celdaList,
#'  \link{plotModuleDecisionTree} for visualizing the tree
#' @examples
#' data(celdaGSim)
#' ## Create models with recursive splitting
#' moduleSplit <- recursiveSplitModule(celdaGSim$counts,
#'   initialL = 3, maxL = 10)
#'
#' ## Build decision tree
#' modTree <- buildModuleDecisionTree(moduleSplit)
#'
#' ## View tree structure
#' print(modTree)
#' @export
buildModuleDecisionTree <- function(x,
    altExpName = "featureSubset",
    minL = NULL,
    maxL = NULL) {

    # Handle SingleCellExperiment input
    if (methods::is(x, "SingleCellExperiment")) {
        if (!altExpName %in% SingleCellExperiment::altExpNames(x)) {
            stop(altExpName, " not in 'altExpNames(x)'")
        }
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        if (!"celda_grid_search" %in% names(S4Vectors::metadata(altExp))) {
            stop("No celda_grid_search found in metadata. ",
                "Make sure x is the result of recursiveSplitModule()")
        }
        celdaList <- S4Vectors::metadata(altExp)$celda_grid_search
    } else if (methods::is(x, "celdaList")) {
        celdaList <- x
    } else {
        stop("x must be a 'celdaList' or 'SingleCellExperiment' object ",
            "from recursiveSplitModule")
    }

    runParams <- runParams(celdaList)
    if (!"L" %in% colnames(runParams)) {
        stop("celdaList must contain models with L parameter (from ",
            "recursiveSplitModule)")
    }

    # Set L range
    availL <- sort(unique(runParams$L))
    if (is.null(minL)) {
        minL <- min(availL)
    }
    if (is.null(maxL)) {
        maxL <- max(availL)
    }

    # Filter to requested L range
    lValues <- availL[availL >= minL & availL <= maxL]
    if (length(lValues) < 2) {
        stop("Need at least 2 L values to build a tree. ",
            "Available L values: ", paste(availL, collapse = ", "))
    }

    # Extract models and module assignments
    models <- list()
    moduleAssignments <- list()
    logLiks <- numeric(length(lValues))

    for (i in seq_along(lValues)) {
        L <- lValues[i]
        idx <- which(runParams$L == L)[1]
        models[[i]] <- celdaList@resList[[idx]]
        moduleAssignments[[i]] <- celdaClusters(models[[i]])$y
        logLiks[i] <- runParams$log_likelihood[idx]
    }
    names(moduleAssignments) <- paste0("L", lValues)

    # Build tree structure
    treeData <- .buildModuleTreeStructure(
        moduleAssignments = moduleAssignments,
        lValues = lValues,
        logLiks = logLiks
    )

    result <- list(
        nodes = treeData$nodes,
        edges = treeData$edges,
        geneAssignments = moduleAssignments,
        lValues = lValues,
        celdaList = celdaList
    )
    class(result) <- "moduleDecisionTree"

    return(result)
}


#' @keywords internal
.buildModuleTreeStructure <- function(moduleAssignments,
    lValues,
    logLiks) {

    nodes <- data.frame(
        nodeId = character(),
        L = integer(),
        moduleLabel = integer(),
        parentNode = character(),
        nGenes = integer(),
        logLik = numeric(),
        stringsAsFactors = FALSE
    )
    edges <- data.frame(
        from = character(),
        to = character(),
        stringsAsFactors = FALSE
    )

    nodeCounter <- 1

    # Create nodes for first L
    L1 <- lValues[1]
    y1 <- moduleAssignments[[1]]
    modules1 <- sort(unique(y1))

    for (mod in modules1) {
        nodeId <- paste0("L", L1, "_M", mod)
        nodes <- rbind(nodes, data.frame(
            nodeId = nodeId,
            L = L1,
            moduleLabel = mod,
            parentNode = NA_character_,
            nGenes = sum(y1 == mod),
            logLik = logLiks[1],
            stringsAsFactors = FALSE
        ))
    }

    # Process subsequent L values
    for (i in seq(2, length(lValues))) {
        Lprev <- lValues[i - 1]
        Lcurr <- lValues[i]
        yPrev <- moduleAssignments[[i - 1]]
        yCurr <- moduleAssignments[[i]]
        modulesCurr <- sort(unique(yCurr))

        # For each module in current L, find which module(s) from previous L
        # it came from
        for (mod in modulesCurr) {
            genesInMod <- which(yCurr == mod)
            prevModules <- unique(yPrev[genesInMod])

            # Find the parent module (the one with most genes in common)
            geneCounts <- table(yPrev[genesInMod])
            parentMod <- as.integer(names(geneCounts)[which.max(geneCounts)])

            nodeId <- paste0("L", Lcurr, "_M", mod)
            parentNodeId <- paste0("L", Lprev, "_M", parentMod)

            nodes <- rbind(nodes, data.frame(
                nodeId = nodeId,
                L = Lcurr,
                moduleLabel = mod,
                parentNode = parentNodeId,
                nGenes = length(genesInMod),
                logLik = logLiks[i],
                stringsAsFactors = FALSE
            ))

            edges <- rbind(edges, data.frame(
                from = parentNodeId,
                to = nodeId,
                stringsAsFactors = FALSE
            ))
        }
    }

    return(list(nodes = nodes, edges = edges))
}


#' @title Plot module decision tree
#' @description Visualizes the hierarchical decision tree of module splitting
#'  using an improved dendrogram-style plot. The visualization shows how modules
#'  split across different L values with a top-down tree layout, curved edges,
#'  color-coded nodes by gene density, and clear labeling.
#' @param x A "moduleDecisionTree" object from
#'  \link{buildModuleDecisionTree}.
#' @param labelModules Logical. Whether to label modules with their module IDs
#'  (e.g., "M1", "M2"). Default TRUE.
#' @param labelGenes Logical. Whether to show gene counts below each node.
#'  Default TRUE.
#' @param plotType Character. Type of plot: "dendrogram" for a hierarchical
#'  tree plot, or "network" for a network-style visualization.
#'  Default "dendrogram".
#' @param ... Additional parameters passed to plotting functions.
#' @return A plot of the module decision tree. The plot features:
#'  \itemize{
#'    \item Top-down layout (higher L values at top)
#'    \item Curved edges showing parent-child relationships
#'    \item Color gradient indicating gene density (light to dark blue)
#'    \item Node sizes proportional to gene counts
#'    \item Grid lines for easy reading of L values
#'    \item Module IDs and gene counts labeled
#'  }
#' @seealso \link{buildModuleDecisionTree}
#' @examples
#' data(celdaGSim)
#' ## Create models with recursive splitting
#' moduleSplit <- recursiveSplitModule(celdaGSim$counts,
#'   initialL = 3, maxL = 10)
#'
#' ## Build and plot decision tree
#' modTree <- buildModuleDecisionTree(moduleSplit)
#' plotModuleDecisionTree(modTree)
#' @export
plotModuleDecisionTree <- function(x,
    labelModules = TRUE,
    labelGenes = TRUE,
    plotType = c("dendrogram", "network"),
    ...) {

    if (!methods::is(x, "moduleDecisionTree")) {
        stop("x must be a 'moduleDecisionTree' object from ",
            "buildModuleDecisionTree")
    }

    plotType <- match.arg(plotType)

    if (plotType == "dendrogram") {
        .plotModuleTreeDendrogram(x,
            labelModules = labelModules,
            labelGenes = labelGenes, ...)
    } else {
        .plotModuleTreeNetwork(x,
            labelModules = labelModules,
            labelGenes = labelGenes, ...)
    }
}


#' @keywords internal
.plotModuleTreeDendrogram <- function(tree,
    labelModules = TRUE,
    labelGenes = TRUE,
    ...) {

    nodes <- tree$nodes
    edges <- tree$edges
    lValues <- tree$lValues

    # Reverse L values so tree flows top to bottom (high L at top)
    lValues <- rev(lValues)

    # Calculate hierarchical positions using a better layout algorithm
    nodePositions <- .calculateTreeLayout(nodes, edges, lValues)

    # Set up plot with better margins
    oldPar <- graphics::par(mar = c(5, 5, 4, 2))
    on.exit(graphics::par(oldPar))

    # Calculate plot dimensions
    xRange <- range(nodePositions$x)
    xPad <- diff(xRange) * 0.1
    yRange <- range(nodePositions$y)
    yPad <- diff(yRange) * 0.15

    graphics::plot.new()
    graphics::plot.window(
        xlim = c(xRange[1] - xPad, xRange[2] + xPad),
        ylim = c(yRange[1] - yPad, yRange[2] + yPad)
    )

    # Add title and labels
    graphics::title(
        main = "Gene Module Decision Tree",
        xlab = "",
        ylab = "Number of Modules (L)",
        cex.main = 1.3,
        font.main = 2
    )

    # Add horizontal grid lines for each L value
    for (L in lValues) {
        yPos <- nodePositions[nodePositions$L == L, "y"][1]
        graphics::abline(h = yPos, col = "gray90", lty = 2, lwd = 0.5)
    }

    # Add L value labels on left
    for (L in lValues) {
        yPos <- nodePositions[nodePositions$L == L, "y"][1]
        graphics::mtext(paste0("L=", L),
            side = 2,
            at = yPos,
            line = 3,
            cex = 0.8,
            las = 1)
    }

    # Color palette - use a better gradient
    colorPal <- grDevices::colorRampPalette(
        c("#E3F2FD", "#2196F3", "#0D47A1")
    )(100)

    # Calculate node sizes based on gene count
    minSize <- 1.5
    maxSize <- 4
    nodeSizes <- minSize + (maxSize - minSize) *
        (nodes$nGenes - min(nodes$nGenes)) /
        (max(nodes$nGenes) - min(nodes$nGenes))
    names(nodeSizes) <- nodes$nodeId

    # Draw edges with curved lines for better visibility
    if (nrow(edges) > 0) {
        for (i in seq_len(nrow(edges))) {
            fromPos <- nodePositions[nodePositions$nodeId == edges$from[i], ]
            toPos <- nodePositions[nodePositions$nodeId == edges$to[i], ]

            # Draw curved edge
            .drawCurvedEdge(
                x0 = fromPos$x, y0 = fromPos$y,
                x1 = toPos$x, y1 = toPos$y,
                col = "gray60",
                lwd = 1.5
            )
        }
    }

    # Draw nodes
    for (i in seq_len(nrow(nodePositions))) {
        nodeId <- nodePositions$nodeId[i]
        nodeData <- nodes[nodes$nodeId == nodeId, ]
        x <- nodePositions$x[i]
        y <- nodePositions$y[i]

        # Color based on gene count
        nGenes <- nodeData$nGenes
        colorIdx <- min(100, max(1, floor(
            (nGenes - min(nodes$nGenes)) /
            (max(nodes$nGenes) - min(nodes$nGenes)) * 99) + 1
        ))
        nodeCol <- colorPal[colorIdx]

        # Draw node
        nodeSize <- nodeSizes[nodeId]
        graphics::points(x, y,
            pch = 21,
            bg = nodeCol,
            cex = nodeSize,
            col = "gray30",
            lwd = 1.5)

        # Add module label inside node if requested
        if (labelModules) {
            modLabel <- paste0("M", nodeData$moduleLabel)
            graphics::text(x, y,
                labels = modLabel,
                cex = 0.5 + nodeSize * 0.15,
                col = if (colorIdx > 50) "white" else "black",
                font = 2)
        }

        # Add gene count as annotation below node
        if (labelGenes) {
            graphics::text(x, y - 0.15,
                labels = paste0(nGenes, " genes"),
                cex = 0.6,
                col = "gray30",
                pos = 1)
        }
    }

    # Add informative legend
    legendText <- c(
        paste0("Modules: ", min(lValues), " → ", max(lValues)),
        paste0("Total genes: ", sum(nodes[nodes$L == min(lValues), "nGenes"])),
        "Node size ~ # genes",
        "Color: gene density"
    )

    graphics::legend("topright",
        legend = legendText,
        bty = "n",
        cex = 0.8,
        text.col = "gray30")

    invisible(nodePositions)
}


#' @keywords internal
.calculateTreeLayout <- function(nodes, edges, lValues) {
    # Better tree layout algorithm
    nodePositions <- data.frame(
        nodeId = character(),
        L = integer(),
        x = numeric(),
        y = numeric(),
        stringsAsFactors = FALSE
    )

    # Assign Y positions (higher L at top)
    yPositions <- seq(length(lValues), 1, by = -1)
    names(yPositions) <- lValues

    # Calculate X positions using hierarchical layout
    for (i in seq_along(lValues)) {
        L <- lValues[i]
        nodesAtL <- nodes[nodes$L == L, ]
        nNodes <- nrow(nodesAtL)

        if (i == 1) {
            # First level - evenly spaced
            xPos <- seq(0, nNodes - 1, length.out = nNodes)
        } else {
            # Position based on parent positions
            xPos <- numeric(nNodes)
            for (j in seq_len(nNodes)) {
                nodeId <- nodesAtL$nodeId[j]
                parentId <- nodesAtL$parentNode[j]

                if (!is.na(parentId)) {
                    # Find parent position
                    parentX <- nodePositions[nodePositions$nodeId == parentId, "x"]
                    # Find siblings (nodes with same parent)
                    siblings <- nodesAtL[nodesAtL$parentNode == parentId, ]
                    nSiblings <- nrow(siblings)
                    siblingIdx <- which(siblings$nodeId == nodeId)

                    # Position relative to parent
                    if (nSiblings == 1) {
                        xPos[j] <- parentX
                    } else {
                        offset <- (siblingIdx - (nSiblings + 1) / 2) * 0.5
                        xPos[j] <- parentX + offset
                    }
                } else {
                    xPos[j] <- j - 1
                }
            }
        }

        # Add to positions
        for (j in seq_len(nNodes)) {
            nodePositions <- rbind(nodePositions, data.frame(
                nodeId = nodesAtL$nodeId[j],
                L = L,
                x = xPos[j],
                y = yPositions[as.character(L)],
                stringsAsFactors = FALSE
            ))
        }
    }

    rownames(nodePositions) <- nodePositions$nodeId
    return(nodePositions)
}


#' @keywords internal
.drawCurvedEdge <- function(x0, y0, x1, y1, col, lwd) {
    # Draw a curved bezier-like edge
    nSteps <- 20

    # Control point for curve
    midY <- (y0 + y1) / 2

    # Generate curve points
    t <- seq(0, 1, length.out = nSteps)

    # Quadratic bezier curve
    xCurve <- (1 - t)^2 * x0 + 2 * (1 - t) * t * ((x0 + x1) / 2) + t^2 * x1
    yCurve <- (1 - t)^2 * y0 + 2 * (1 - t) * t * midY + t^2 * y1

    # Draw curve as line segments
    graphics::lines(xCurve, yCurve, col = col, lwd = lwd)
}


#' @keywords internal
.plotModuleTreeNetwork <- function(tree,
    labelModules = TRUE,
    labelGenes = TRUE,
    ...) {

    message("Network visualization requires 'igraph' package. ",
        "Using dendrogram plot instead.")
    .plotModuleTreeDendrogram(tree,
        labelModules = labelModules,
        labelGenes = labelGenes, ...)
}


#' @title Get module lineage
#' @description Trace the lineage of a specific module back through the
#'  decision tree to its ancestral modules at lower L values.
#' @param tree A "moduleDecisionTree" object from
#'  \link{buildModuleDecisionTree}.
#' @param L Integer. The L value where the module exists.
#' @param module Integer. The module label to trace.
#' @return A data frame showing the lineage of the module, with columns:
#'  \itemize{
#'    \item \code{L}: Number of modules at each level
#'    \item \code{moduleLabel}: Module label at this level
#'    \item \code{nGenes}: Number of genes in the module
#'    \item \code{logLik}: Log-likelihood at this level
#'  }
#' @examples
#' data(celdaGSim)
#' moduleSplit <- recursiveSplitModule(celdaGSim$counts,
#'   initialL = 3, maxL = 10)
#' modTree <- buildModuleDecisionTree(moduleSplit)
#'
#' ## Trace lineage of module 5 at L=10
#' lineage <- getModuleLineage(modTree, L = 10, module = 5)
#' print(lineage)
#' @export
getModuleLineage <- function(tree, L, module) {
    if (!methods::is(tree, "moduleDecisionTree")) {
        stop("tree must be a 'moduleDecisionTree' object")
    }

    nodes <- tree$nodes
    nodeId <- paste0("L", L, "_M", module)

    if (!nodeId %in% nodes$nodeId) {
        stop("Module ", module, " at L=", L, " not found in tree")
    }

    lineage <- data.frame(
        L = integer(),
        moduleLabel = integer(),
        nGenes = integer(),
        logLik = numeric(),
        stringsAsFactors = FALSE
    )

    currentNode <- nodeId
    while (!is.na(currentNode)) {
        nodeData <- nodes[nodes$nodeId == currentNode, ]
        lineage <- rbind(data.frame(
            L = nodeData$L,
            moduleLabel = nodeData$moduleLabel,
            nGenes = nodeData$nGenes,
            logLik = nodeData$logLik,
            stringsAsFactors = FALSE
        ), lineage)

        currentNode <- nodeData$parentNode
    }

    return(lineage)
}


#' @title Get module descendants
#' @description Find all descendant modules that arose from splitting a
#'  specific module at higher L values.
#' @param tree A "moduleDecisionTree" object from
#'  \link{buildModuleDecisionTree}.
#' @param L Integer. The L value where the ancestral module exists.
#' @param module Integer. The module label to find descendants for.
#' @return A data frame showing all descendant modules, with columns:
#'  \itemize{
#'    \item \code{L}: Number of modules at each level
#'    \item \code{moduleLabel}: Module label at this level
#'    \item \code{nGenes}: Number of genes in the module
#'    \item \code{generationsFromAncestor}: How many splits from the ancestor
#'  }
#' @examples
#' data(celdaGSim)
#' moduleSplit <- recursiveSplitModule(celdaGSim$counts,
#'   initialL = 3, maxL = 10)
#' modTree <- buildModuleDecisionTree(moduleSplit)
#'
#' ## Find descendants of module 1 at L=3
#' descendants <- getModuleDescendants(modTree, L = 3, module = 1)
#' print(descendants)
#' @export
getModuleDescendants <- function(tree, L, module) {
    if (!methods::is(tree, "moduleDecisionTree")) {
        stop("tree must be a 'moduleDecisionTree' object")
    }

    nodes <- tree$nodes
    edges <- tree$edges
    ancestorId <- paste0("L", L, "_M", module)

    if (!ancestorId %in% nodes$nodeId) {
        stop("Module ", module, " at L=", L, " not found in tree")
    }

    descendants <- data.frame(
        L = integer(),
        moduleLabel = integer(),
        nGenes = integer(),
        generationsFromAncestor = integer(),
        stringsAsFactors = FALSE
    )

    # Breadth-first search for descendants
    toExplore <- data.frame(
        nodeId = ancestorId,
        generation = 0,
        stringsAsFactors = FALSE
    )

    while (nrow(toExplore) > 0) {
        currentNode <- toExplore[1, ]
        toExplore <- toExplore[-1, , drop = FALSE]

        nodeData <- nodes[nodes$nodeId == currentNode$nodeId, ]
        if (currentNode$generation > 0) {  # Don't include the ancestor itself
            descendants <- rbind(descendants, data.frame(
                L = nodeData$L,
                moduleLabel = nodeData$moduleLabel,
                nGenes = nodeData$nGenes,
                generationsFromAncestor = currentNode$generation,
                stringsAsFactors = FALSE
            ))
        }

        # Add children to explore
        children <- edges[edges$from == currentNode$nodeId, "to"]
        if (length(children) > 0) {
            toExplore <- rbind(toExplore, data.frame(
                nodeId = children,
                generation = currentNode$generation + 1,
                stringsAsFactors = FALSE
            ))
        }
    }

    return(descendants)
}


#' @title Print method for moduleDecisionTree
#' @param x A moduleDecisionTree object
#' @param ... Additional arguments (not used)
#' @export
print.moduleDecisionTree <- function(x, ...) {
    cat("Module Decision Tree\n")
    cat("====================\n\n")
    cat("L values:", paste(x$lValues, collapse = ", "), "\n")
    cat("Total nodes:", nrow(x$nodes), "\n")
    cat("Total edges:", nrow(x$edges), "\n\n")

    cat("Modules at each L:\n")
    for (L in x$lValues) {
        nMods <- sum(x$nodes$L == L)
        cat("  L =", L, ":", nMods, "modules\n")
    }

    cat("\nUse plotModuleDecisionTree() to visualize\n")
    cat("Use getModuleLineage() to trace module ancestry\n")
    cat("Use getModuleDescendants() to find module descendants\n")
}
