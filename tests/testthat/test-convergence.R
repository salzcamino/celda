context("Advanced Convergence Detection")

# Test helper functions
test_that(".calculateARI works correctly", {
    # Identical clusterings should have ARI = 1
    x <- c(1, 1, 2, 2, 3, 3)
    y <- c(1, 1, 2, 2, 3, 3)
    expect_equal(celda:::.calculateARI(x, y), 1)

    # Completely different clusterings should have ARI < 1
    x <- c(1, 1, 1, 2, 2, 2)
    y <- c(1, 2, 1, 2, 1, 2)
    ari <- celda:::.calculateARI(x, y)
    expect_true(ari < 1)
    expect_true(ari > -1)

    # All same cluster should give ARI = 1
    x <- c(1, 1, 1, 1)
    y <- c(2, 2, 2, 2)
    expect_equal(celda:::.calculateARI(x, y), 1)

    # Random clusterings should have ARI near 0
    set.seed(12345)
    x <- sample(1:5, 100, replace = TRUE)
    y <- sample(1:5, 100, replace = TRUE)
    ari <- celda:::.calculateARI(x, y)
    expect_true(abs(ari) < 0.3)  # Should be close to 0

    # Test error on mismatched lengths
    expect_error(celda:::.calculateARI(c(1, 2, 3), c(1, 2)),
                "must have the same length")

    # Test empty vectors
    expect_true(is.na(celda:::.calculateARI(integer(0), integer(0))))
})


test_that(".checkConvergence_Advanced detects LL convergence", {
    # Create converged log-likelihood series
    llHistory <- c(100, 150, 180, 190, 195, 197, 198, 198.5,
                  198.7, 198.8, 198.85)
    zHistory <- NULL  # Not checking stability

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = zHistory,
        yHistory = NULL,
        iter = 11,
        stopIter = 10,
        relTol = 1e-3,
        checkStability = FALSE,
        minIter = 10
    )

    expect_true(result$converged)
    expect_true(result$llConverged)
    expect_true(is.na(result$zStable))
    expect_true(is.na(result$yStable))
    expect_match(result$reason, "Log-likelihood stable")
})


test_that(".checkConvergence_Advanced does not converge prematurely", {
    # Create non-converged log-likelihood series (still improving)
    llHistory <- c(100, 150, 180, 190, 200, 210, 220, 230, 240, 250, 260)
    zHistory <- NULL

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = zHistory,
        yHistory = NULL,
        iter = 11,
        stopIter = 10,
        relTol = 1e-5,
        checkStability = FALSE,
        minIter = 10
    )

    expect_false(result$converged)
    expect_false(result$llConverged)
    expect_match(result$reason, "Not converged")
})


test_that(".checkConvergence_Advanced detects cluster stability", {
    # Create stable cluster history (clusters not changing)
    z1 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
    z2 <- z1  # Identical
    z3 <- z1
    z4 <- z1
    z5 <- z1

    zHistory <- list(z1, z2, z3, z4, z5)

    # Converged LL
    llHistory <- c(100, 150, 180, 190, 195, 197, 198, 198.5,
                  198.7, 198.8, 198.85)

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = zHistory,
        yHistory = NULL,
        iter = 11,
        stopIter = 10,
        relTol = 1e-3,
        checkStability = TRUE,
        ariThreshold = 0.99,
        minIter = 10
    )

    expect_true(result$converged)
    expect_true(result$llConverged)
    expect_true(result$zStable)
    expect_match(result$reason, "Cell clusters are stable")
})


test_that(".checkConvergence_Advanced detects unstable clusters", {
    # Create unstable cluster history (clusters still changing)
    z1 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
    z2 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
    z3 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
    z4 <- c(1, 2, 2, 2, 3, 3, 4, 4, 5, 5)  # One cell moved
    z5 <- c(1, 2, 3, 2, 3, 3, 4, 4, 5, 5)  # Another cell moved

    zHistory <- list(z1, z2, z3, z4, z5)

    # Converged LL but unstable clusters
    llHistory <- c(100, 150, 180, 190, 195, 197, 198, 198.5,
                  198.7, 198.8, 198.85)

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = zHistory,
        yHistory = NULL,
        iter = 11,
        stopIter = 10,
        relTol = 1e-3,
        checkStability = TRUE,
        ariThreshold = 0.99,
        minIter = 10
    )

    # Should not converge because clusters unstable
    expect_false(result$converged)
    expect_true(result$llConverged)
    expect_false(result$zStable)
    expect_match(result$reason, "Cell clusters still changing")
})


test_that(".checkConvergence_Advanced handles both z and y stability", {
    # Create stable z and y histories
    z1 <- c(1, 1, 2, 2, 3, 3)
    y1 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)

    zHistory <- list(z1, z1, z1, z1, z1)
    yHistory <- list(y1, y1, y1, y1, y1)

    llHistory <- c(100, 150, 180, 190, 195, 197, 198, 198.5,
                  198.7, 198.8, 198.85)

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = zHistory,
        yHistory = yHistory,
        iter = 11,
        stopIter = 10,
        relTol = 1e-3,
        checkStability = TRUE,
        ariThreshold = 0.99,
        minIter = 10
    )

    expect_true(result$converged)
    expect_true(result$llConverged)
    expect_true(result$zStable)
    expect_true(result$yStable)
    expect_match(result$reason, "Cell clusters and gene modules are stable")
})


test_that(".checkConvergence_Advanced handles extended LL convergence", {
    # LL converged for 2*stopIter iterations, even with minor cluster changes
    # Create a very stable LL for 20+ iterations
    llHistory <- rep(200, 25)
    llHistory[1:5] <- c(100, 150, 180, 190, 195)

    # Clusters slightly unstable
    z1 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
    z2 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
    z3 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
    z4 <- c(1, 2, 2, 2, 3, 3, 4, 4, 5, 5)  # One cell moved
    z5 <- c(1, 2, 3, 2, 3, 3, 4, 4, 5, 5)  # Another cell moved
    zHistory <- list(z1, z2, z3, z4, z5)

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = zHistory,
        yHistory = NULL,
        iter = 25,
        stopIter = 10,
        relTol = 1e-5,
        checkStability = TRUE,
        ariThreshold = 0.99,
        minIter = 10
    )

    # Should converge due to extended LL stability
    expect_true(result$converged)
    expect_match(result$reason, "stable for 20 iterations")
})


test_that(".checkConvergence_Advanced respects minimum iterations", {
    # Even if LL looks converged, don't stop before minIter
    llHistory <- rep(200, 8)

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = NULL,
        yHistory = NULL,
        iter = 8,
        stopIter = 5,
        relTol = 1e-5,
        checkStability = FALSE,
        minIter = 10
    )

    expect_false(result$converged)
    expect_match(result$reason, "minimum required iterations")
})


test_that(".checkConvergence_Advanced handles insufficient history", {
    # Too few iterations
    llHistory <- c(100, 150, 180)

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = NULL,
        yHistory = NULL,
        iter = 3,
        stopIter = 10,
        relTol = 1e-5,
        checkStability = FALSE,
        minIter = 1
    )

    expect_false(result$converged)
    expect_match(result$reason, "Insufficient log-likelihood history")
})


test_that(".checkConvergence_Simple works as expected", {
    # Improved LL
    llHistory <- c(100, 150, 180, 190)
    currentLL <- 195

    result <- celda:::.checkConvergence_Simple(
        llHistory = llHistory,
        currentLL = currentLL,
        stopIter = 10,
        numIterWithoutImprovement = 3
    )

    expect_true(result$improved)
    expect_equal(result$numIterWithoutImprovement, 1L)

    # Not improved LL
    llHistory <- c(100, 150, 180, 190, 195)
    currentLL <- 194

    result <- celda:::.checkConvergence_Simple(
        llHistory = llHistory,
        currentLL = currentLL,
        stopIter = 10,
        numIterWithoutImprovement = 3
    )

    expect_false(result$improved)
    expect_equal(result$numIterWithoutImprovement, 4L)
})


test_that("Advanced convergence catches oscillating LL", {
    # LL oscillates but doesn't really improve
    llHistory <- c(100, 150, 180, 190, 189, 191, 190, 192, 191, 193, 192)

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = NULL,
        yHistory = NULL,
        iter = 11,
        stopIter = 10,
        relTol = 1e-4,
        checkStability = FALSE,
        minIter = 10
    )

    # Should converge because relative change is small
    expect_true(result$converged)
    expect_true(result$llConverged)
})


test_that("Cluster history handles NULL entries", {
    # Sometimes history list may have NULL entries
    z1 <- c(1, 1, 2, 2, 3, 3)
    zHistory <- list(NULL, z1, NULL, z1, z1, z1, z1)

    llHistory <- rep(200, 11)

    result <- celda:::.checkConvergence_Advanced(
        llHistory = llHistory,
        zHistory = zHistory,
        yHistory = NULL,
        iter = 11,
        stopIter = 10,
        relTol = 1e-5,
        checkStability = TRUE,
        minIter = 10
    )

    # Should handle NULL entries gracefully
    expect_true(result$converged)
})


test_that("ARI calculation matches mclust when available", {
    skip_if_not_installed("mclust")

    x <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
    y <- c(1, 1, 2, 2, 2, 3, 3, 3, 3)

    celdaARI <- celda:::.calculateARI(x, y)
    mclustARI <- mclust::adjustedRandIndex(x, y)

    expect_equal(celdaARI, mclustARI, tolerance = 1e-10)
})


test_that("Edge cases for ARI", {
    # Single element
    expect_equal(celda:::.calculateARI(1, 1), 1)

    # Two elements, same cluster
    expect_equal(celda:::.calculateARI(c(1, 1), c(2, 2)), 1)

    # Two elements, different clusters
    expect_equal(celda:::.calculateARI(c(1, 2), c(1, 2)), 1)
    expect_equal(celda:::.calculateARI(c(1, 2), c(2, 1)), -1)
})
