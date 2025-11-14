library(celda)
context("Testing module decision tree functions")

# Note: These tests require running recursiveSplitModule which can be slow
# Tests are kept minimal to ensure the basic functionality works

test_that(desc = "buildModuleDecisionTree requires valid input", {
    expect_error(
        buildModuleDecisionTree(matrix(0)),
        "x must be a 'celdaList' or 'SingleCellExperiment' object"
    )
})

test_that(desc = "buildModuleDecisionTree requires L parameter", {
    # Create a mock celdaList without L parameter
    mockList <- methods::new("celdaList",
        runParams = data.frame(K = c(3, 4)),
        resList = list(),
        countChecksum = "test")

    expect_error(
        buildModuleDecisionTree(mockList),
        "celdaList must contain models with L parameter"
    )
})

test_that(desc = "buildModuleDecisionTree handles SCE without grid search", {
    # Create a mock SCE without celda_grid_search in metadata
    mockSCE <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 10, ncol = 10))
    )
    SingleCellExperiment::altExp(mockSCE, "featureSubset") <- mockSCE

    expect_error(
        buildModuleDecisionTree(mockSCE),
        "No celda_grid_search found in metadata"
    )
})

test_that(desc = "getModuleLineage requires valid inputs", {
    # Create a minimal mock tree
    mockTree <- list(
        nodes = data.frame(
            nodeId = c("L3_M1", "L4_M1"),
            L = c(3, 4),
            moduleLabel = c(1, 1),
            parentNode = c(NA, "L3_M1"),
            nGenes = c(10, 5),
            logLik = c(-100, -95),
            stringsAsFactors = FALSE
        ),
        edges = data.frame(
            from = "L3_M1",
            to = "L4_M1",
            stringsAsFactors = FALSE
        ),
        geneAssignments = list(),
        lValues = c(3, 4),
        celdaList = NULL
    )
    class(mockTree) <- "moduleDecisionTree"

    expect_error(
        getModuleLineage(matrix(0), L = 3, module = 1),
        "tree must be a 'moduleDecisionTree' object"
    )

    expect_error(
        getModuleLineage(mockTree, L = 5, module = 1),
        "Module 1 at L=5 not found in tree"
    )

    # Valid call should work
    lineage <- getModuleLineage(mockTree, L = 4, module = 1)
    expect_true(is.data.frame(lineage))
    expect_equal(nrow(lineage), 2)
    expect_equal(lineage$L, c(3, 4))
})

test_that(desc = "getModuleDescendants requires valid inputs", {
    # Create a minimal mock tree
    mockTree <- list(
        nodes = data.frame(
            nodeId = c("L3_M1", "L4_M1", "L4_M2"),
            L = c(3, 4, 4),
            moduleLabel = c(1, 1, 2),
            parentNode = c(NA, "L3_M1", "L3_M1"),
            nGenes = c(10, 5, 5),
            logLik = c(-100, -95, -95),
            stringsAsFactors = FALSE
        ),
        edges = data.frame(
            from = c("L3_M1", "L3_M1"),
            to = c("L4_M1", "L4_M2"),
            stringsAsFactors = FALSE
        ),
        geneAssignments = list(),
        lValues = c(3, 4),
        celdaList = NULL
    )
    class(mockTree) <- "moduleDecisionTree"

    expect_error(
        getModuleDescendants(matrix(0), L = 3, module = 1),
        "tree must be a 'moduleDecisionTree' object"
    )

    expect_error(
        getModuleDescendants(mockTree, L = 5, module = 1),
        "Module 1 at L=5 not found in tree"
    )

    # Valid call should work
    descendants <- getModuleDescendants(mockTree, L = 3, module = 1)
    expect_true(is.data.frame(descendants))
    expect_equal(nrow(descendants), 2)  # Should have 2 descendants
})

test_that(desc = "print.moduleDecisionTree works", {
    mockTree <- list(
        nodes = data.frame(
            nodeId = c("L3_M1", "L4_M1"),
            L = c(3, 4),
            moduleLabel = c(1, 1),
            parentNode = c(NA, "L3_M1"),
            nGenes = c(10, 5),
            logLik = c(-100, -95),
            stringsAsFactors = FALSE
        ),
        edges = data.frame(
            from = "L3_M1",
            to = "L4_M1",
            stringsAsFactors = FALSE
        ),
        geneAssignments = list(),
        lValues = c(3, 4),
        celdaList = NULL
    )
    class(mockTree) <- "moduleDecisionTree"

    # Should not error
    expect_output(print(mockTree), "Module Decision Tree")
})

test_that(desc = "plotModuleDecisionTree requires moduleDecisionTree input", {
    expect_error(
        plotModuleDecisionTree(matrix(0)),
        "x must be a 'moduleDecisionTree' object"
    )
})
