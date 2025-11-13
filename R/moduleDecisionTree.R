#' @title Build decision tree from recursive module splitting
#' @description Constructs a hierarchical decision tree structure that traces
#'  how feature modules are split during recursive module splitting. Each node
#'  in the tree represents a module at a specific L value, with edges showing
#'  parent-child relationships when modules are split.
#' @param celdaList A \linkS4class{celdaList} object returned by
#'  \link{recursiveSplitModule}.
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
buildModuleDecisionTree <- function(celdaList,
    minL = NULL,
    maxL = NULL) {

    if (!methods::is(celdaList, "celdaList")) {
        stop("celdaList must be a 'celdaList' object from recursiveSplitModule")
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
#'  using a dendrogram-style plot.
#' @param x A "moduleDecisionTree" object from
#'  \link{buildModuleDecisionTree}.
#' @param labelModules Logical. Whether to label modules at terminal nodes.
#'  Default TRUE.
#' @param labelGenes Logical. Whether to show gene counts at nodes.
#'  Default TRUE.
#' @param plotType Character. Type of plot: "dendrogram" for a hierarchical
#'  tree plot, or "network" for a network-style visualization.
#'  Default "dendrogram".
#' @param ... Additional parameters passed to plotting functions.
#' @return A plot of the module decision tree.
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

    # Set up plot area
    nLevels <- length(lValues)
    maxModules <- max(table(nodes$L))

    graphics::plot.new()
    graphics::plot.window(xlim = c(0, maxModules + 1),
        ylim = c(min(lValues) - 1, max(lValues) + 1))
    graphics::title(main = "Module Decision Tree",
        xlab = "Module Position",
        ylab = "Number of Modules (L)")
    graphics::axis(2, at = lValues)

    # Calculate positions for each node
    nodePositions <- data.frame(
        nodeId = character(),
        x = numeric(),
        y = numeric(),
        stringsAsFactors = FALSE
    )

    for (L in lValues) {
        nodesAtL <- nodes[nodes$L == L, ]
        nNodes <- nrow(nodesAtL)
        xPos <- seq(1, maxModules, length.out = nNodes)

        for (i in seq_len(nNodes)) {
            nodePositions <- rbind(nodePositions, data.frame(
                nodeId = nodesAtL$nodeId[i],
                x = xPos[i],
                y = L,
                stringsAsFactors = FALSE
            ))
        }
    }
    rownames(nodePositions) <- nodePositions$nodeId

    # Draw edges
    for (i in seq_len(nrow(edges))) {
        fromPos <- nodePositions[edges$from[i], ]
        toPos <- nodePositions[edges$to[i], ]
        graphics::segments(fromPos$x, fromPos$y, toPos$x, toPos$y,
            col = "gray40", lwd = 1.5)
    }

    # Draw nodes
    for (i in seq_len(nrow(nodePositions))) {
        nodeId <- nodePositions$nodeId[i]
        nodeData <- nodes[nodes$nodeId == nodeId, ]
        x <- nodePositions$x[i]
        y <- nodePositions$y[i]

        # Node color based on size
        nGenes <- nodeData$nGenes
        nodeCol <- grDevices::colorRampPalette(c("lightblue", "darkblue"))(100)[
            min(100, floor(nGenes / max(nodes$nGenes) * 100) + 1)
        ]

        graphics::points(x, y, pch = 21, bg = nodeCol, cex = 2, col = "black")

        # Labels
        if (labelModules && y == max(lValues)) {
            graphics::text(x, y, labels = nodeData$moduleLabel,
                pos = 3, cex = 0.7)
        }
        if (labelGenes) {
            graphics::text(x, y, labels = nGenes, cex = 0.6, col = "white")
        }
    }

    # Add legend
    graphics::legend("topright",
        legend = c("Node size = # genes"),
        pch = 21,
        pt.bg = "lightblue",
        pt.cex = 1.5,
        bty = "n")
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
