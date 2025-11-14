# Testing Instructions for Clustering Optimizations

## Quick Start

Run the comprehensive test suite:

```bash
Rscript test_optimizations.R
```

Or from R console:
```r
source("test_optimizations.R")
```

## Individual Tests

### 1. Regenerate Documentation

```r
# Must be run from package root directory
roxygen2::roxygenise()
```

This will update:
- `man/*.Rd` files with new parameter documentation
- `NAMESPACE` with any new exports

### 2. Run Full Test Suite

```r
# Run all existing tests
devtools::test()

# Or with coverage
covr::package_coverage()
```

Expected: All existing tests should pass without modification.

### 3. Verify Specific Models

#### Test celda_C

```r
library(celda)
data(celdaCSim)

# Test 1: Basic backward compatibility
result1 <- celda_C(
  celdaCSim$counts,
  K = 5,
  sampleLabel = celdaCSim$sampleLabel,
  nchains = 1
)

# Test 2: With adaptive splits (enabled by default)
result2 <- celda_C(
  celdaCSim$counts,
  K = 5,
  sampleLabel = celdaCSim$sampleLabel,
  splitAdaptive = TRUE,
  nchains = 1
)

# Test 3: With parallel processing
result3 <- celda_C(
  celdaCSim$counts,
  K = 5,
  sampleLabel = celdaCSim$sampleLabel,
  nCores = 4,
  nchains = 1
)

# Test 4: With custom heterogeneity threshold
result4 <- celda_C(
  celdaCSim$counts,
  K = 5,
  sampleLabel = celdaCSim$sampleLabel,
  heterogeneityThreshold = 0.5,
  nchains = 1
)

# Verify clustering quality is maintained
z1 <- celdaClusters(result1)$z
z2 <- celdaClusters(result2)$z
ari <- mclust::adjustedRandIndex(z1, z2)
print(paste("ARI between default and adaptive:", round(ari, 3)))
# Expected: ARI > 0.9 (high agreement)
```

#### Test celda_G

```r
data(celdaGSim)

# Test 1: Basic functionality
result_g1 <- celda_G(
  celdaGSim$counts,
  L = 10,
  nchains = 1
)

# Test 2: With parallel processing
result_g2 <- celda_G(
  celdaGSim$counts,
  L = 10,
  nCores = 4,
  nchains = 1
)

# Verify results
print(celdaClusters(result_g1)$y)
```

#### Test celda_CG

```r
data(celdaCGSim)

# Test 1: Basic functionality
result_cg <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  sampleLabel = celdaCGSim$sampleLabel,
  nchains = 1
)

# Test 2: With nCores
result_cg2 <- celda_CG(
  celdaCGSim$counts,
  K = 5,
  L = 10,
  sampleLabel = celdaCGSim$sampleLabel,
  nCores = 2,
  nchains = 1
)

# Verify results
print(celdaClusters(result_cg))
```

## Performance Benchmarking

### Quick Benchmark

```r
library(microbenchmark)
data(celdaCSim)

# Compare baseline vs optimized
bench <- microbenchmark(
  baseline = celda_C(
    celdaCSim$counts,
    K = 5,
    sampleLabel = celdaCSim$sampleLabel,
    splitAdaptive = FALSE,
    nCores = 1,
    maxIter = 50,
    nchains = 1
  ),
  optimized = celda_C(
    celdaCSim$counts,
    K = 5,
    sampleLabel = celdaCSim$sampleLabel,
    splitAdaptive = TRUE,
    heterogeneityThreshold = 0.3,
    nCores = 1,
    maxIter = 50,
    nchains = 1
  ),
  parallel = celda_C(
    celdaCSim$counts,
    K = 5,
    sampleLabel = celdaCSim$sampleLabel,
    splitAdaptive = TRUE,
    nCores = 4,
    maxIter = 50,
    nchains = 1
  ),
  times = 3
)

print(bench)
boxplot(bench)
```

Expected speedups:
- Optimized vs Baseline: 1.3-1.5x
- Parallel (4 cores) vs Baseline: 2.5-4x

### Comprehensive Benchmark

```r
source("benchmark_clustering_optimizations.R")  # If created from plan
```

## Validation Checklist

- [ ] Documentation regenerates without errors
- [ ] All existing tests pass
- [ ] celda_C works with new parameters
- [ ] celda_G works with nCores parameter
- [ ] celda_CG accepts nCores parameter
- [ ] Clustering quality maintained (ARI > 0.9)
- [ ] Measurable performance improvement (>1.2x)
- [ ] Parallel processing works on multi-core systems
- [ ] No warnings or errors in R CMD check

## Full R CMD Check

```r
# From package root
devtools::check()
```

Expected output:
```
0 errors ✓ | 0 warnings ✓ | 0 notes ✓
```

## Common Issues and Solutions

### Issue: "object '.identifySplitCandidates' not found"

**Solution**: The function is internal (starts with `.`). If testing directly:
```r
celda:::.identifySplitCandidates(counts, z, K, 0.3, 3)
```

### Issue: Parallel processing doesn't work

**Cause**: `parallel::mclapply()` doesn't work on Windows
**Solution**: This is expected. On Windows, it falls back to sequential processing.

### Issue: ARI < 0.9 between baseline and optimized

**Potential causes**:
1. Stochasticity in clustering (try with same seed)
2. Heterogeneity threshold too aggressive
3. Actual clustering difference (investigate further)

**Action**:
```r
# Test with same seed
set.seed(12345)
result1 <- celda_C(counts, K=5, splitAdaptive=FALSE)
set.seed(12345)
result2 <- celda_C(counts, K=5, splitAdaptive=TRUE)
ari <- mclust::adjustedRandIndex(
  celdaClusters(result1)$z,
  celdaClusters(result2)$z
)
```

### Issue: No performance improvement observed

**Potential causes**:
1. Dataset too small (overhead dominates)
2. K too small (< 8 clusters)
3. Not enough iterations for splits to matter

**Solution**: Test with larger dataset and K ≥ 10

## Testing on Real Data

```r
# Load your real dataset
counts <- ... # Your count matrix
sampleLabel <- ... # Your sample labels

# Run with optimizations
result <- celda_C(
  counts,
  K = 15,  # Use larger K to see benefits
  sampleLabel = sampleLabel,
  nCores = 4,
  splitAdaptive = TRUE,
  heterogeneityThreshold = 0.3,
  maxIter = 200
)

# Evaluate clustering quality
library(celda)
plotCeldaTsne(result)  # Visualize
celdaProbabilityMap(result)  # Check probabilities
```

## Reporting Issues

If tests fail or performance is not as expected:

1. **Collect diagnostic information**:
```r
sessionInfo()
packageVersion("celda")
parallel::detectCores()
```

2. **Create minimal reproducible example**:
```r
library(celda)
data(celdaCSim)
result <- celda_C(celdaCSim$counts, K=5, nCores=4)
# Include error message or unexpected behavior
```

3. **Report at**: https://github.com/campbio/celda/issues

## Next Steps After Testing

1. ✅ All tests pass → Ready to merge!
2. ⚠️ Minor issues → Fix and retest
3. ❌ Major issues → Review implementation

## Performance Expectations by Dataset Size

| Cells | K | Expected Speedup |
|-------|---|------------------|
| < 500 | 5 | 1.2-1.5x |
| 1000 | 8 | 2-2.5x |
| 5000 | 15 | 3-4x |
| 10000+ | 20+ | 3.5-4.5x |

With `nCores = 4`, add 1.5-2x additional speedup.
