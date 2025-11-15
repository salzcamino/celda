# Celda Clustering Improvements - Integration Testing Checklist

**Branch**: `claude/optimize-recursive-split-module-011CUmo1eBNfyKg2tRcVFYew`
**Status**: All implementations complete, ready for R testing
**Date**: 2025-11-15

---

## Overview

This checklist covers testing for **7 major clustering algorithm improvements** plus **performance optimizations** and **report enhancements**. All code is committed and pushed to the remote branch.

---

## Pre-Testing Setup

### 1. Environment Setup

```r
# Install/update package from branch
devtools::install_github("salzcamino/celda",
  ref = "claude/optimize-recursive-split-module-011CUmo1eBNfyKg2tRcVFYew")

# Load package
library(celda)
library(SingleCellExperiment)
library(testthat)

# Check version
packageVersion("celda")  # Should be 1.19.0 or development version
```

### 2. Generate Documentation

```r
# Generate roxygen2 documentation
setwd("/path/to/celda")
roxygen2::roxygenise()

# Check for warnings
devtools::check_man()
```

### 3. Verify New Files Exist

```r
# R functions (should all exist)
file.exists(c(
  "R/feature_weights.R",
  "R/marker_initialization.R",
  "R/graph_split.R",
  "R/convergence.R",
  "R/cluster_validation.R",
  "R/consensus_clustering.R"
))

# Test files (should all exist)
file.exists(c(
  "tests/testthat/test-feature_weights.R",
  "tests/testthat/test-marker_initialization.R",
  "tests/testthat/test-adaptive_subclusters.R",
  "tests/testthat/test-graph_split.R",
  "tests/testthat/test-convergence.R",
  "tests/testthat/test-cluster_validation.R",
  "tests/testthat/test-consensus_clustering.R"
))
```

---

## Unit Testing

### 4. Run All Tests

```r
# Run full test suite
devtools::test()

# Expected: All tests pass (>66 new tests added)
# No errors or warnings
```

### 5. Run Individual Test Files

Test each improvement separately:

```r
# Phase 1.1: Adaptive Feature Weighting
testthat::test_file("tests/testthat/test-feature_weights.R")
# Expected: 11 tests pass

# Phase 1.2: Marker-Guided Initialization
testthat::test_file("tests/testthat/test-marker_initialization.R")
# Expected: 16 tests pass

# Phase 1.3: Adaptive K/L Subcluster Selection
testthat::test_file("tests/testthat/test-adaptive_subclusters.R")
# Expected: 16+ tests pass

# Phase 2.1: Graph-Based Split Heuristic
testthat::test_file("tests/testthat/test-graph_split.R")
# Expected: 10 tests pass

# Phase 2.2: Advanced Convergence Detection
testthat::test_file("tests/testthat/test-convergence.R")
# Expected: 13 tests pass

# Phase 2.3: Internal Cluster Validation Metrics
testthat::test_file("tests/testthat/test-cluster_validation.R")
# Expected: 15+ tests pass

# Phase 3.1: Consensus Clustering
testthat::test_file("tests/testthat/test-consensus_clustering.R")
# Expected: 15 tests pass
```

### 6. Test Coverage

```r
# Generate coverage report
covr::package_coverage()

# Expected: >90% coverage for new code
# Check that all new functions are tested
```

---

## Functional Testing

### 7. Test Backward Compatibility

**Critical**: Existing code should work unchanged

```r
# Load example data
data(celdaCGSim)

# OLD API - Should work identically
sce_old <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  nchains = 1,
  maxIter = 50
)

# Verify clusters assigned
table(celdaClusters(sce_old)$z)  # Should have 5 cell clusters
table(celdaClusters(sce_old)$y)  # Should have 10 gene modules

# Check log-likelihood
metadata(sce_old)$celda_parameters$finalLogLik  # Should be numeric
```

### 8. Test Each New Feature Individually

#### 8.1 Adaptive Feature Weighting

```r
data(celdaCGSim)

sce_weighted <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  featureReweighting = TRUE,
  reweightInterval = 5,
  nchains = 1,
  maxIter = 50
)

# Verify:
# - Completes without errors
# - Gene weights are stored in metadata
# - Clustering quality is similar or better
```

#### 8.2 Marker-Guided Initialization

```r
data(celdaCGSim)

# Define marker genes (example for simulated data)
markers <- list(
  "CellType1" = c("Gene_1", "Gene_2", "Gene_3"),
  "CellType2" = c("Gene_10", "Gene_11", "Gene_12"),
  "CellType3" = c("Gene_20", "Gene_21", "Gene_22")
)

sce_marker <- celda_C(
  celdaCGSim$counts,
  K = 5,
  markerGenes = markers,
  zInitialize = "split",
  nchains = 1,
  maxIter = 50
)

# Verify:
# - Completes without errors
# - Converges faster than random initialization
# - Cluster assignments make biological sense
```

#### 8.3 Adaptive Subclusters

```r
data(celdaCGSim)

sce_adaptive <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  adaptiveSubclusters = TRUE,
  nchains = 1,
  maxIter = 50
)

# Verify:
# - Completes without errors
# - Uses data-driven subcluster selection
# - Initialization quality is good
```

#### 8.4 Graph-Based Splitting

```r
data(celdaCGSim)

# First, create reduced dimensions (UMAP)
library(scater)
sce <- SingleCellExperiment(assays = list(counts = celdaCGSim$counts))
sce <- logNormCounts(sce)
sce <- runUMAP(sce)
umap_coords <- reducedDim(sce, "UMAP")

# Run with graph-based splitting
sce_graph <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  useGraphBasedSplit = TRUE,
  reducedDimForSplit = umap_coords,
  nchains = 1,
  maxIter = 50
)

# Verify:
# - Completes without errors
# - Uses graph-based split heuristic
# - Subcluster detection is improved
```

#### 8.5 Advanced Convergence Detection

```r
data(celdaCGSim)

sce_advanced_conv <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  convergenceMethod = "advanced",
  convergenceRelTol = 1e-5,
  checkClusterStability = TRUE,
  nchains = 1,
  maxIter = 200,
  verbose = TRUE
)

# Verify:
# - Completes without errors
# - Converges earlier than "simple" method
# - Prints convergence reason
# - Results are high quality
```

#### 8.6 Cluster Validation Metrics

```r
data(celdaCGSim)

sce <- celda_CG(celdaCGSim$counts, K = 5, L = 10, nchains = 1)

# Calculate validation metrics
metrics <- .calculateClusterMetrics(
  counts = celdaCGSim$counts,
  z = celdaClusters(sce)$z,
  y = celdaClusters(sce)$y,
  metrics = "all"
)

# Verify:
# - Returns silhouette scores
# - Returns Calinski-Harabasz index
# - Returns Davies-Bouldin index
# - Returns module coherence
# - All metrics are numeric and reasonable
```

#### 8.7 Consensus Clustering

```r
data(celdaCGSim)

sce_consensus <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  nchains = 5,  # Multiple chains required
  useConsensus = TRUE,
  consensusMethod = "cooccurrence",
  minConsensusAgreement = 0.7,
  maxIter = 50
)

# Verify:
# - Completes without errors
# - Produces consensus clusters
# - Confidence scores are provided
# - Results are more robust than single chain
```

### 9. Test Combined Features

**Real-world usage with all improvements enabled:**

```r
data(celdaCGSim)

# Prepare reduced dimensions
library(scater)
sce <- SingleCellExperiment(assays = list(counts = celdaCGSim$counts))
sce <- logNormCounts(sce)
sce <- runUMAP(sce)
umap_coords <- reducedDim(sce, "UMAP")

# Define markers (if available)
markers <- list(
  "CellType1" = c("Gene_1", "Gene_2"),
  "CellType2" = c("Gene_10", "Gene_11")
)

# Run with ALL improvements
sce_full <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  # Initialization
  markerGenes = markers,
  adaptiveSubclusters = TRUE,
  # Clustering
  featureReweighting = TRUE,
  useGraphBasedSplit = TRUE,
  reducedDimForSplit = umap_coords,
  # Convergence
  convergenceMethod = "advanced",
  # Robustness
  nchains = 3,
  useConsensus = TRUE,
  maxIter = 100,
  verbose = TRUE
)

# Verify:
# - Completes without errors
# - Converges efficiently
# - High-quality results
# - All features work together harmoniously
```

---

## Performance Testing

### 10. Run Benchmarking Scripts

```r
# Comprehensive benchmark
source("benchmark_optimizations.R")

# Expected:
# - recursiveSplitModule: 2-4x speedup with nCores=4
# - recursiveSplitCell: 2-4x speedup with nCores=4
# - decontX: 3-5x speedup with nCores=4, nThreads=4
```

```r
# Simple workflow benchmark
source("benchmark_simple.R")

# Expected:
# - Overall workflow speedup with parallelization
# - Results identical between serial and parallel
```

### 11. Performance Profiling

```r
library(profvis)

data(celdaCGSim)

# Profile with improvements disabled
profvis({
  sce_old <- celda_CG(celdaCGSim$counts, K = 5, L = 10, nchains = 1, maxIter = 50)
})

# Profile with improvements enabled
profvis({
  sce_new <- celda_CG(
    celdaCGSim$counts,
    K = 5,
    L = 10,
    featureReweighting = TRUE,
    convergenceMethod = "advanced",
    nchains = 1,
    maxIter = 50
  )
})

# Compare:
# - Total execution time (should be 10-30% faster)
# - Memory usage (should be comparable or better)
# - Identify any bottlenecks
```

---

## Integration Testing

### 12. Test with Real Datasets

#### PBMC 3K (Small, Well-Characterized)

```r
# Load PBMC 3K dataset
library(TENxPBMCData)
sce <- TENxPBMCData("pbmc3k")

# Select features
sce <- selectFeatures(sce, minCount = 3)

# Run with improvements
sce <- celda_CG(
  sce,
  K = 8,  # Expected: ~8 cell types in PBMC
  L = 15,
  featureReweighting = TRUE,
  convergenceMethod = "advanced",
  nchains = 3,
  useConsensus = TRUE,
  maxIter = 200,
  verbose = TRUE
)

# Verify:
# - Identifies known cell types (T cells, B cells, monocytes, etc.)
# - Cluster quality metrics are good
# - Biological interpretation makes sense
```

#### Simulated Data with Ground Truth

```r
# Generate simulated data
sim <- simulateCellsV2(
  C = 10,
  G = 200,
  L = 20,
  N = 1000,
  K = 10
)

# Run celda with improvements
sce_sim <- celda_CG(
  sim$counts,
  K = 10,
  L = 20,
  featureReweighting = TRUE,
  convergenceMethod = "advanced",
  nchains = 3,
  maxIter = 100
)

# Calculate ARI vs. ground truth
library(mclust)
ari_cells <- adjustedRandIndex(sim$z, celdaClusters(sce_sim)$z)
ari_genes <- adjustedRandIndex(sim$y, celdaClusters(sce_sim)$y)

# Expected:
# - ARI > 0.8 for both cells and genes
# - Better than baseline celda without improvements
```

### 13. Test celda_C and celda_G Separately

```r
data(celdaCGSim)

# Test celda_C with convergence detection
sce_c <- celda_C(
  celdaCGSim$counts,
  K = 5,
  convergenceMethod = "advanced",
  nchains = 1,
  maxIter = 200
)

# Test celda_G with convergence detection
sce_g <- celda_G(
  celdaCGSim$counts,
  L = 10,
  convergenceMethod = "advanced",
  nchains = 1,
  maxIter = 200
)

# Verify both complete successfully
```

---

## Report Testing

### 14. Test Enhanced Reports

```r
library(celda)
data(celdaCGSim)

# Generate dataset
sce <- SingleCellExperiment(assays = list(counts = celdaCGSim$counts))
sce <- selectFeatures(sce)

# Run celda
sce <- celda_CG(sce, K = 5, L = 10, nchains = 1)

# Generate report with new features
reportCeldaCGPlotResults(
  sce = sce,
  runDiffExp = TRUE,
  diffExpMethod = "wilcox",
  topMarkers = 10,
  runCellTypeAnnotation = TRUE,
  outputFile = "test_report.html",
  outputDir = tempdir()
)

# Verify:
# - Report generates without errors
# - Module plots have correct aspect ratios
# - Differential expression section is present
# - Cell type annotation section is present (if ref data available)
# - All plots render correctly
```

---

## Package Check

### 15. R CMD check

```bash
# From command line
cd /path/to/celda
R CMD build .
R CMD check --as-cran celda_*.tar.gz
```

Expected results:
- 0 Errors
- 0 Warnings
- Notes are acceptable (package size, etc.)

### 16. BiocCheck

```r
# Bioconductor compliance check
BiocCheck::BiocCheck(".")
```

Expected results:
- Passes all Bioconductor requirements
- No critical issues

---

## Memory and Performance Validation

### 17. Memory Profiling

```r
library(profmem)

data(celdaCGSim)

# Profile memory usage
mem_profile <- profmem({
  sce <- celda_CG(
    celdaCGSim$counts,
    K = 5,
    L = 10,
    featureReweighting = TRUE,
    convergenceMethod = "advanced",
    nchains = 1,
    maxIter = 50
  )
})

# Analyze peak memory
total_mem <- sum(mem_profile$bytes, na.rm = TRUE)
print(paste("Total memory allocated:", format(total_mem, big.mark = ",")))

# Verify: No excessive memory allocation
```

### 18. Scalability Testing

```r
# Test with increasing dataset sizes
sizes <- c(500, 1000, 2000, 5000)
times <- numeric(length(sizes))

for (i in seq_along(sizes)) {
  # Generate data
  sim <- simulateCellsV2(C = 5, G = 100, L = 10, N = sizes[i], K = 5)

  # Time execution
  times[i] <- system.time({
    sce <- celda_CG(
      sim$counts,
      K = 5,
      L = 10,
      convergenceMethod = "advanced",
      nchains = 1,
      maxIter = 50,
      verbose = FALSE
    )
  })["elapsed"]
}

# Plot scaling
plot(sizes, times, type = "b",
     xlab = "Number of cells",
     ylab = "Time (seconds)",
     main = "celda_CG Scalability")

# Verify: Time scales reasonably (should be sub-quadratic)
```

---

## Expected Quality Improvements

Based on implementations, we expect to see:

### Clustering Quality Metrics

| Metric | Expected Improvement |
|--------|---------------------|
| ARI vs. ground truth | +15-25% |
| Silhouette score | +10-20% |
| Calinski-Harabasz index | +15-25% |
| Module coherence | +10-15% |

### Performance Metrics

| Operation | Expected Speedup |
|-----------|-----------------|
| Iterations to convergence | -20-40% |
| Total runtime | +10-30% faster |
| Parallel recursive splits | 2-4x |
| decontX with batches | 3-5x |

### Robustness Metrics

| Metric | Expected Improvement |
|--------|---------------------|
| Chain-to-chain consistency (ARI) | +20-30% |
| Stability across runs | +15-25% |

---

## Troubleshooting

### Common Issues

1. **Missing dependencies**
   ```r
   # Install optional dependencies for graph-based features
   install.packages(c("mclust", "diptest"))
   ```

2. **Slow convergence**
   - Try `convergenceMethod = "advanced"` for faster detection
   - Increase `stopIter` if clusters are unstable
   - Use `markerGenes` for better initialization

3. **High memory usage**
   - Use sparse matrix input (`dgCMatrix`)
   - Reduce `nchains` if memory is limited
   - Consider mini-batch approach (Phase 3.2, optional)

4. **Poor clustering quality**
   - Enable `featureReweighting = TRUE`
   - Provide `markerGenes` if known
   - Use `useGraphBasedSplit = TRUE` with UMAP/t-SNE
   - Increase `nchains` and use `useConsensus = TRUE`

---

## Sign-Off Checklist

Before merging to main/master:

- [ ] All unit tests pass (`devtools::test()`)
- [ ] All backward compatibility tests pass
- [ ] All new features work individually
- [ ] All new features work together
- [ ] Performance benchmarks show expected improvements
- [ ] Real dataset tests show biological validity
- [ ] R CMD check passes with no errors/warnings
- [ ] BiocCheck passes
- [ ] Memory usage is acceptable
- [ ] Scalability is reasonable
- [ ] Documentation is complete (`roxygen2::roxygenise()`)
- [ ] NEWS.md is updated
- [ ] Vignettes demonstrate new features (optional but recommended)

---

## Next Steps After Testing

1. **If all tests pass:**
   - Create pull request to main/master branch
   - Request code review from maintainers
   - Address any review comments
   - Merge when approved

2. **If issues found:**
   - Document issues in GitHub Issues
   - Fix critical bugs immediately
   - Reassess non-critical issues for future releases
   - Re-run testing after fixes

3. **Future enhancements (optional):**
   - Phase 3.2: Stochastic mini-batch for very large datasets
   - Phase 3.3: scRNA-seq biological priors (zero-inflation, batch-aware)
   - Additional vignettes and tutorials
   - Performance optimization for specific use cases

---

## Contact

For questions or issues during testing:
- GitHub Issues: https://github.com/salzcamino/celda/issues
- Branch: `claude/optimize-recursive-split-module-011CUmo1eBNfyKg2tRcVFYew`

---

**End of Checklist**
