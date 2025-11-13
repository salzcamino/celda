# CLAUDE.md - AI Assistant Guide for celda

**Last Updated**: 2025-11-13
**Version**: 1.18.2
**Package Type**: R/Bioconductor

This document provides comprehensive guidance for AI assistants working with the celda codebase.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Core Architecture](#core-architecture)
4. [Development Workflows](#development-workflows)
5. [Coding Conventions](#coding-conventions)
6. [Testing Strategy](#testing-strategy)
7. [Common Tasks](#common-tasks)
8. [Key Files Reference](#key-files-reference)
9. [Important Notes for AI Assistants](#important-notes-for-ai-assistants)

---

## Project Overview

### What is celda?

**celda** (CEllular Latent Dirichlet Allocation) is a Bioconductor package for:
- **Single-cell RNA-seq analysis** using Bayesian hierarchical models
- **Bi-clustering**: Simultaneous clustering of genes AND cells
- **DecontX**: Estimation and removal of ambient RNA contamination
- Extension of Latent Dirichlet Allocation (LDA) from text mining to genomics

### Key Features
- Three core clustering models: `celda_C`, `celda_G`, `celda_CG`
- High-performance C/C++ implementation for computationally intensive operations
- Comprehensive visualization suite
- Integration with SingleCellExperiment and Bioconductor ecosystem
- Parallel processing support for large datasets

### Repository Information
- **Organization**: campbio
- **Repository**: celda
- **License**: MIT + file LICENSE
- **Bug Reports**: https://github.com/campbio/celda/issues
- **Wiki**: https://github.com/campbio/celda/wiki

---

## Repository Structure

```
celda/
├── R/                          # R source code (40 files, ~20,498 lines)
│   ├── aaa.R                   # S4 class definitions (MUST load first)
│   ├── celda_C.R               # Cell clustering model
│   ├── celda_G.R               # Gene clustering model
│   ├── celda_CG.R              # Combined bi-clustering model
│   ├── decon.R                 # DecontX contamination removal
│   ├── celdaGridSearch.R       # Parameter optimization
│   ├── moduleDecisionTree.R    # Module hierarchy analysis
│   └── ...                     # Supporting functions
│
├── src/                        # C/C++ compiled code (8 files, ~1,716 lines)
│   ├── RcppExports.cpp         # Auto-generated Rcpp bindings
│   ├── DecontX.cpp             # DecontX C++ implementation
│   ├── cG_calcGibbsProbY.cpp   # Gibbs sampling for gene clustering
│   ├── matrixSums.c            # Matrix aggregation functions
│   └── ...                     # Performance-critical operations
│
├── tests/testthat/             # Unit tests
│   ├── test-celda_C.R
│   ├── test-celda_G.R
│   ├── test-celda_CG.R
│   ├── test-decon.R
│   ├── test-moduleDecisionTree.R
│   └── ...
│
├── man/                        # Auto-generated documentation (roxygen2)
├── vignettes/                  # Package vignettes
│   ├── celda.Rmd               # Main celda vignette
│   └── decontX.Rmd             # DecontX vignette
│
├── data/                       # Example datasets
├── inst/                       # Additional package resources
├── docs/                       # pkgdown website
│
├── .github/workflows/          # CI/CD workflows
│   ├── check-standard.yaml     # R CMD check (multi-OS)
│   └── BioC-check.yaml         # Bioconductor compliance check
│
├── DESCRIPTION                 # Package metadata
├── NAMESPACE                   # Exported functions and imports
├── README.md                   # User-facing documentation
├── NEWS.md                     # Change log
├── BENCHMARKING.md             # Performance optimization guide
├── MODULE_DECISION_TREE_README.md  # Module decision tree documentation
├── benchmark_optimizations.R   # Comprehensive benchmarking script
└── benchmark_simple.R          # Simple workflow benchmark
```

### File Organization Principle

**Model-Centric Organization**: The codebase is organized around the three core models:
- **celda_C**: Cell population clustering
- **celda_G**: Gene module clustering
- **celda_CG**: Combined bi-clustering

Supporting functions are grouped by category (initialization, computation, visualization, etc.) rather than by model.

---

## Core Architecture

### Three Core Models

#### 1. celda_C (Cell Clustering)
- **File**: `/home/user/celda/R/celda_C.R` (978 lines)
- **Purpose**: Clusters **cells** into K cell populations
- **Algorithm**: Gibbs sampling or EM
- **Key Parameters**:
  - `K`: Number of cell clusters
  - `alpha`, `beta`: Dirichlet concentration parameters
  - `sampleLabel`: Sample identifiers for batch effects
- **Output**: Vector `z` of cell cluster assignments
- **Key Functions**:
  - `.celda_C()`: Main clustering loop
  - `.cCCalcGibbsProbZ()`: Gibbs sampling probabilities
  - `.cCCalcEMProbZ()`: EM algorithm probabilities
  - `.cCCalcLL()`: Log-likelihood calculation
  - `.cCSplitZ()`: Cluster splitting heuristic

#### 2. celda_G (Gene Clustering)
- **File**: `/home/user/celda/R/celda_G.R` (853 lines)
- **Purpose**: Clusters **genes** into L feature modules
- **Algorithm**: Gibbs sampling
- **Key Parameters**:
  - `L`: Number of gene modules
  - `beta`, `delta`, `gamma`: Dirichlet concentration parameters
- **Output**: Vector `y` of gene module assignments
- **Key Functions**:
  - `.celda_G()`: Main clustering loop
  - `.cGCalcGibbsProbY()`: Gibbs sampling (calls C++)
  - `.cGCalcLL()`: Log-likelihood calculation
  - `.cGSplitY()`: Module splitting heuristic

#### 3. celda_CG (Bi-Clustering)
- **File**: `/home/user/celda/R/celda_CG.R` (1,102 lines)
- **Purpose**: Simultaneous clustering of **both cells AND genes**
- **Algorithm**: Iterative Gibbs sampling (alternates between cells and genes)
- **Key Parameters**: Combines all parameters from celda_C and celda_G
- **Output**: Both `z` (cell clusters) and `y` (gene modules)
- **Key Functions**:
  - `.celda_CG()`: Main integrated clustering
  - `.cCGCalcLL()`: Combined log-likelihood
  - `.cCGDecomposeCounts()`: Count matrix decomposition

### S4 Class Hierarchy

```r
# Defined in R/aaa.R (loaded first due to naming convention)

celdaModel (base S4 class)
    ├── params: List of model parameters
    ├── names: Feature and sample names
    ├── completeLogLik: Complete data log-likelihood history
    ├── finalLogLik: Final log-likelihood
    └── ...

celda_C extends celdaModel
    ├── clusters: Cell cluster assignments (z)
    └── sampleLabel: Sample labels

celda_G extends celdaModel
    └── clusters: Gene module assignments (y)

celda_CG extends both celda_C and celda_G
    └── clusters: Both z and y assignments
```

**Why S4?** S4 provides:
- Formal class definitions with type checking
- Multiple dispatch (methods can specialize on multiple arguments)
- Required for Bioconductor standards

### R and C/C++ Integration

**Performance-Critical Operations in C/C++**:

| C/C++ Function | File | Purpose | Called From (R) |
|----------------|------|---------|----------------|
| `_rowSumByGroup` | `matrixSums.c` | Aggregate rows by cluster | `.rowSumByGroup()` |
| `_colSumByGroup` | `matrixSums.c` | Aggregate columns by cluster | `.colSumByGroup()` |
| `cG_CalcGibbsProbY` | `cG_calcGibbsProbY.cpp` | Gibbs probabilities for genes | `.cGCalcGibbsProbY()` |
| `eigenMatMultInt` | `eigenMatMultInt.cpp` | Matrix multiplication | `.countsTimesProbs()` |
| `fastNormPropLog` | `matrixNorm.cpp` | Fast normalized log proportions | EM algorithms |
| `decontXEM` | `DecontX.cpp` | Contamination estimation | `decontX()` |

**Call Chain Example** (Cell Clustering):
```
celda_C()                    [R: Public API]
  ↓
.celdaCWithSeed()            [R: Seed wrapper]
  ↓
.celda_C()                   [R: Main algorithm loop]
  ↓
.cCCalcGibbsProbZ()          [R: Probability calculation]
  ↓
.countsTimesProbs()          [R: Dispatcher]
  ↓
eigenMatMultInt()            [C++: Via Rcpp]
```

**Dispatcher Pattern**: Many R functions check input type and dispatch to appropriate C/C++ implementation:
```r
.rowSumByGroup <- function(x, group) {
  if (is(x, "dgCMatrix")) {
    # Sparse matrix path
  } else if (is.integer(x)) {
    .Call("_rowSumByGroup", x, group)  # C implementation for integers
  } else {
    .Call("_rowSumByGroup_numeric", x, group)  # C implementation for numerics
  }
}
```

---

## Development Workflows

### Local Development Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/campbio/celda.git
   cd celda
   ```

2. **Install dependencies** (from R console):
   ```r
   # Install BiocManager if needed
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")

   # Install all dependencies
   BiocManager::install(c(
       "SingleCellExperiment", "Matrix", "Rcpp", "RcppEigen",
       "testthat", "roxygen2", "knitr", "rmarkdown"
   ))
   ```

3. **Build and install the package**:
   ```r
   # From project root
   devtools::document()     # Generate documentation
   devtools::build()        # Build package
   devtools::install()      # Install locally
   ```

### Running Tests

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-celda_C.R")

# Run tests with coverage
covr::package_coverage()
```

**Test Command Line**:
```bash
R CMD check --no-manual --no-build-vignettes .
```

### Continuous Integration

**Two CI Workflows**:

1. **R-CMD-check** (`.github/workflows/check-standard.yaml`)
   - Runs on: macOS, Windows, Ubuntu
   - Triggered on: Push/PR to `devel` or `master` branches
   - Tests: Standard R CMD check

2. **BioC-check** (`.github/workflows/BioC-check.yaml`)
   - Runs on: macOS only
   - Triggered on: Push/PR to `devel` or `master` branches
   - Tests: Bioconductor compliance (`BiocCheck`)

### Documentation Generation

**Roxygen2 Workflow**:
```r
# 1. Write roxygen2 comments in R files
#' @title Function Title
#' @description Description
#' @param x Input parameter
#' @return Return value description
#' @export
#' @examples
#' example_code()

# 2. Generate documentation
roxygen2::roxygenise()

# 3. Check documentation
?functionName
```

**Vignette Building**:
```r
# Build vignettes
devtools::build_vignettes()

# View vignettes
vignette("celda")
vignette("decontX")
```

**pkgdown Website**:
```r
# Build website
pkgdown::build_site()

# Preview locally
pkgdown::preview_site()
```

### Performance Benchmarking

The package includes comprehensive benchmarking scripts:

```bash
# Run comprehensive benchmarks
Rscript benchmark_optimizations.R

# Run simple workflow benchmark
Rscript benchmark_simple.R
```

See `BENCHMARKING.md` for detailed information on performance optimization and parallelization.

---

## Coding Conventions

### Naming Conventions

#### Function Names

| Pattern | Visibility | Example | Usage |
|---------|-----------|---------|-------|
| `functionName()` | Public/exported | `celda_C()`, `celdaGridSearch()` | User-facing API |
| `.functionName()` | Private/internal | `.celda_C()`, `.cCCalcLL()` | Internal implementation |
| `.cC*()` | Internal | `.cCCalcGibbsProbZ()` | celda_C model internals |
| `.cG*()` | Internal | `.cGCalcGibbsProbY()` | celda_G model internals |
| `.cCG*()` | Internal | `.cCGCalcLL()` | celda_CG model internals |
| `.cD*()` | Internal | `.cDCalcEMDecontamination()` | DecontX internals |
| `*WithSeed()` | Wrapper | `.celdaCWithSeed()` | Seed-wrapping functions |

**Convention Rules**:
1. **Public functions**: camelCase, no leading dot
2. **Internal functions**: camelCase with leading dot (`.`)
3. **Model-specific prefixes**: `.cC*` (celda_C), `.cG*` (celda_G), `.cCG*` (celda_CG)
4. **C/C++ functions**: Use underscores for C (`_rowSumByGroup`), camelCase for C++ (`eigenMatMultInt`)

#### Variable Names

| Variable | Meaning | Type |
|----------|---------|------|
| `z` | Cell cluster assignments | Integer vector |
| `y` | Gene module assignments | Integer vector |
| `K` | Number of cell clusters | Integer |
| `L` | Number of gene modules | Integer |
| `counts` | Count matrix | Matrix (features × cells) |
| `s` | Sample labels | Factor/character vector |
| `nG`, `nM` | Number of genes/cells | Integer |
| `alpha`, `beta`, `gamma`, `delta` | Dirichlet concentration parameters | Numeric |
| `ll` | Log-likelihood | Numeric |
| `m*` | Count aggregations by group | Matrix/vector (e.g., `mCPByS`) |
| `n*` | Count aggregations by group | Matrix/vector (e.g., `nGByCP`) |

### Code Style

**Follow the Celda Development Coding Style Guide** (see Wiki):
- Indentation: 4 spaces (no tabs)
- Line length: Aim for ≤80 characters
- Braces: Opening brace on same line, closing brace on new line
- Assignment: Use `<-` not `=`
- Spacing: Space after commas, around operators

**Example**:
```r
.cCCalcLL <- function(counts, s, z, K, alpha, beta) {
    nG <- nrow(counts)

    # Calculate log-likelihood
    ll <- sum(
        lgamma(nGByCP + beta) - lgamma(beta),
        lgamma(nCP + alpha) - lgamma(alpha)
    )

    return(ll)
}
```

### S4 Method Definition Pattern

```r
# Generic definition
setGeneric("functionName", function(x, ...) {
    standardGeneric("functionName")
})

# Method for SingleCellExperiment
setMethod("functionName", signature(x = "SingleCellExperiment"),
    function(x, ...) {
        # Implementation
    }
)

# Method for matrix
setMethod("functionName", signature(x = "matrix"),
    function(x, ...) {
        # Implementation
    }
)
```

### Documentation Standards

**Every exported function must have**:
1. `@title`: Brief title
2. `@description`: Detailed description
3. `@param`: All parameters documented
4. `@return`: Return value described
5. `@export`: Mark for export
6. `@examples`: Working examples

**Example**:
```r
#' @title Cluster cells using celda_C
#' @description Clusters cells into K populations using Bayesian hierarchical model
#' @param counts Integer matrix of counts (genes × cells)
#' @param K Integer. Number of cell clusters.
#' @param alpha Numeric. Dirichlet concentration parameter. Default 1.
#' @param beta Numeric. Dirichlet concentration parameter. Default 1.
#' @param maxIter Integer. Maximum number of Gibbs sampling iterations. Default 200.
#' @param ... Additional parameters
#' @return A celda_C S4 object containing cluster assignments and model parameters
#' @export
#' @examples
#' data(celdaCMod)
#' result <- celda_C(celdaCMod$counts, K = 5)
celda_C <- function(counts, K, alpha = 1, beta = 1, maxIter = 200, ...) {
    # Implementation
}
```

### Error Handling

**Use informative error messages**:
```r
if (!is.matrix(counts) && !is(counts, "dgCMatrix")) {
    stop("'counts' must be a matrix or dgCMatrix (sparse matrix)")
}

if (K < 2) {
    stop("'K' must be at least 2. Use K = 2 or higher.")
}

if (any(counts < 0)) {
    stop("'counts' cannot contain negative values")
}
```

---

## Testing Strategy

### Test Organization

Tests are organized by module in `tests/testthat/`:

| Test File | Coverage |
|-----------|----------|
| `test-celda_C.R` | Cell clustering model |
| `test-celda_G.R` | Gene clustering model |
| `test-celda_CG.R` | Bi-clustering model |
| `test-decon.R` | DecontX contamination removal |
| `test-celda-functions.R` | Utility functions |
| `test-matrixSums.R` | Matrix aggregation functions |
| `test-moduleDecisionTree.R` | Module decision tree functions |
| `test-intialize_cluster.R` | Cluster initialization |
| `test-with_seed.R` | Seed handling |

### Test Patterns

**1. Basic Functionality Test**:
```r
test_that("celda_C clusters cells correctly", {
    data(celdaCMod)
    result <- celda_C(celdaCMod$counts, K = 5, maxIter = 10)

    expect_s4_class(result, "celda_C")
    expect_equal(length(celdaClusters(result)), ncol(celdaCMod$counts))
    expect_true(all(celdaClusters(result) %in% 1:5))
})
```

**2. Input Validation Test**:
```r
test_that("celda_C validates input parameters", {
    data(celdaCMod)

    # K must be >= 2
    expect_error(
        celda_C(celdaCMod$counts, K = 1),
        "K must be at least 2"
    )

    # counts must be non-negative
    badCounts <- celdaCMod$counts
    badCounts[1, 1] <- -1
    expect_error(
        celda_C(badCounts, K = 5),
        "counts cannot contain negative values"
    )
})
```

**3. Reproducibility Test**:
```r
test_that("celda_C produces reproducible results with same seed", {
    data(celdaCMod)

    result1 <- celda_C(celdaCMod$counts, K = 5, seed = 123)
    result2 <- celda_C(celdaCMod$counts, K = 5, seed = 123)

    expect_identical(celdaClusters(result1), celdaClusters(result2))
})
```

**4. Performance Test**:
```r
test_that("Parallel execution produces same results as serial", {
    data(celdaCMod)

    resultSerial <- recursiveSplitModule(celdaCMod$counts, initialL = 5, maxL = 10, nCores = 1)
    resultParallel <- recursiveSplitModule(celdaCMod$counts, initialL = 5, maxL = 10, nCores = 4)

    # Results should be identical (same seed)
    expect_identical(resultSerial, resultParallel)
})
```

### Running Tests

```r
# Run all tests
devtools::test()

# Run with coverage
covr::package_coverage()

# Run specific file
testthat::test_file("tests/testthat/test-celda_C.R")

# Run single test
testthat::test_file("tests/testthat/test-celda_C.R",
                    filter = "celda_C clusters cells correctly")
```

---

## Common Tasks

### Task 1: Adding a New Function to an Existing Model

**Example**: Add a helper function to celda_C

1. **Add the function to the appropriate R file** (e.g., `R/celda_C.R`):
   ```r
   #' @title Calculate cell purity scores
   #' @description Calculates purity score for each cell cluster
   #' @param celdaMod A celda_C model object
   #' @return Numeric vector of purity scores
   #' @export
   calculateCellPurity <- function(celdaMod) {
       # Implementation
   }
   ```

2. **Add internal helper if needed** (prefix with `.`):
   ```r
   .calculatePurityHelper <- function(z, counts) {
       # Internal implementation
   }
   ```

3. **Document with roxygen2**:
   ```r
   roxygen2::roxygenise()
   ```

4. **Add unit tests** (`tests/testthat/test-celda_C.R`):
   ```r
   test_that("calculateCellPurity returns correct scores", {
       data(celdaCMod)
       result <- celda_C(celdaCMod$counts, K = 5)
       purity <- calculateCellPurity(result)

       expect_length(purity, 5)  # K clusters
       expect_true(all(purity >= 0 & purity <= 1))
   })
   ```

5. **Test and check**:
   ```r
   devtools::test()
   devtools::check()
   ```

### Task 2: Optimizing Performance with C/C++

**When to use C/C++**: Operations in tight loops, matrix operations, repetitive calculations

**Steps**:

1. **Profile to identify bottleneck**:
   ```r
   profvis::profvis({
       result <- celda_C(largeCounts, K = 10)
   })
   ```

2. **Write C++ function** (`src/myFunction.cpp`):
   ```cpp
   #include <Rcpp.h>
   using namespace Rcpp;

   // [[Rcpp::export]]
   NumericVector myFastFunction(NumericMatrix x, IntegerVector groups) {
       int n = x.nrow();
       NumericVector result(n);

       for (int i = 0; i < n; i++) {
           // Fast C++ implementation
       }

       return result;
   }
   ```

3. **Rebuild package** (auto-generates `RcppExports.R` and `RcppExports.cpp`):
   ```r
   devtools::clean_dll()
   devtools::load_all()
   ```

4. **Call from R** (wrapper for type checking):
   ```r
   .myFastFunctionWrapper <- function(x, groups) {
       if (!is.matrix(x)) {
           stop("x must be a matrix")
       }
       myFastFunction(x, groups)  # Calls C++ function
   }
   ```

5. **Benchmark improvement**:
   ```r
   microbenchmark::microbenchmark(
       old = oldSlowFunction(x, groups),
       new = .myFastFunctionWrapper(x, groups),
       times = 100
   )
   ```

### Task 3: Adding a New Model

**Example**: Create `celda_D` for a new clustering approach

1. **Define S4 class** (`R/aaa.R`):
   ```r
   setClass("celda_D",
       contains = "celdaModel",
       slots = c(
           clusters = "integer",  # New cluster assignments
           newParam = "numeric"   # Model-specific parameter
       )
   )
   ```

2. **Implement core model** (`R/celda_D.R`):
   ```r
   #' @export
   celda_D <- function(counts, D, ...) {
       # Public API
       .celdaDWithSeed(counts, D, ...)
   }

   .celdaDWithSeed <- function(counts, D, seed = 12345, ...) {
       # Seed wrapper
       withr::with_seed(seed, {
           .celda_D(counts, D, ...)
       })
   }

   .celda_D <- function(counts, D, maxIter = 200, ...) {
       # Main algorithm implementation
   }
   ```

3. **Add S4 methods**:
   ```r
   setMethod("celdaClusters", signature(celdaMod = "celda_D"),
       function(celdaMod) {
           return(celdaMod@clusters)
       }
   )
   ```

4. **Document thoroughly**:
   ```r
   roxygen2::roxygenise()
   ```

5. **Add comprehensive tests** (`tests/testthat/test-celda_D.R`)

6. **Add to vignettes** if user-facing

### Task 4: Working with Module Decision Trees

**The module decision tree functionality was recently added**. See `MODULE_DECISION_TREE_README.md` for details.

**Basic workflow**:
```r
library(celda)

# 1. Perform recursive splitting
moduleSplit <- recursiveSplitModule(
    counts,
    initialL = 5,
    maxL = 20,
    nCores = 4  # Parallel processing
)

# 2. Build decision tree
modTree <- buildModuleDecisionTree(moduleSplit)

# 3. Visualize
plotModuleDecisionTree(modTree)

# 4. Analyze lineage
lineage <- getModuleLineage(modTree, L = 15, module = 7)

# 5. Find descendants
descendants <- getModuleDescendants(modTree, L = 5, module = 2)
```

### Task 5: Updating Documentation After Changes

1. **Update roxygen2 comments** in R source files

2. **Regenerate documentation**:
   ```r
   roxygen2::roxygenise()
   ```

3. **Check for warnings**:
   ```r
   devtools::check()
   ```

4. **Rebuild pkgdown site**:
   ```r
   pkgdown::build_site()
   ```

5. **Update NEWS.md** with changes:
   ```markdown
   # celda v1.18.3 (2024-XX-XX)
   * Added new function `calculateCellPurity()` for cluster quality assessment
   * Fixed bug in `decontX()` when using sparse matrices
   * Performance improvement: 2x speedup in `recursiveSplitModule()`
   ```

---

## Key Files Reference

### Essential Files to Understand

| File | Purpose | When to Modify |
|------|---------|---------------|
| `DESCRIPTION` | Package metadata, dependencies | Adding dependencies, version bumps |
| `NAMESPACE` | Exported functions, imports | Auto-generated (don't edit manually) |
| `R/aaa.R` | S4 class definitions | Adding new model classes |
| `R/celda_C.R` | Cell clustering model | Modifying cell clustering |
| `R/celda_G.R` | Gene clustering model | Modifying gene clustering |
| `R/celda_CG.R` | Bi-clustering model | Modifying combined model |
| `R/decon.R` | DecontX implementation | Contamination removal changes |
| `R/RcppExports.R` | Rcpp-generated bindings | Auto-generated (don't edit) |
| `src/RcppExports.cpp` | Rcpp-generated C++ bindings | Auto-generated (don't edit) |
| `tests/testthat.R` | Test runner | Rarely (test configuration) |

### Configuration Files

| File | Purpose |
|------|---------|
| `.Rbuildignore` | Files to ignore during R CMD build |
| `.gitignore` | Files to ignore in git |
| `_pkgdown.yml` | pkgdown website configuration |
| `.github/workflows/*.yaml` | CI/CD configuration |

### Documentation Files

| File | Audience | Purpose |
|------|----------|---------|
| `README.md` | Users | Quick start, installation |
| `NEWS.md` | Users | Change log |
| `BENCHMARKING.md` | Developers | Performance optimization guide |
| `MODULE_DECISION_TREE_README.md` | Developers | Module decision tree implementation |
| `CLAUDE.md` (this file) | AI Assistants | Comprehensive development guide |

---

## Important Notes for AI Assistants

### Critical Rules

1. **NEVER edit auto-generated files**:
   - `NAMESPACE` (generated by roxygen2)
   - `R/RcppExports.R` (generated by Rcpp)
   - `src/RcppExports.cpp` (generated by Rcpp)
   - `man/*.Rd` files (generated by roxygen2)

2. **Always use roxygen2 for documentation**:
   - Edit roxygen2 comments in R source files
   - Run `roxygen2::roxygenise()` to regenerate documentation
   - Never edit `.Rd` files directly

3. **S4 classes must be defined in `aaa.R`**:
   - R loads files alphabetically
   - `aaa.R` ensures classes are defined before use
   - Other files depend on these class definitions

4. **Use consistent naming conventions**:
   - Public: `functionName()`
   - Internal: `.functionName()`
   - Model-specific: `.cC*()`, `.cG*()`, `.cCG*()`

5. **Test rigorously**:
   - Add tests for every new function
   - Test edge cases and error conditions
   - Ensure reproducibility with seeds

6. **Follow Bioconductor standards**:
   - Use S4 classes, not S3
   - Support SingleCellExperiment objects
   - Use BiocStyle for vignettes
   - Pass `BiocCheck()`

### Common Pitfalls

1. **Modifying C++ without recompiling**:
   ```r
   # Always run after C++ changes:
   devtools::clean_dll()
   devtools::load_all()
   ```

2. **Forgetting to export new functions**:
   ```r
   # Add to roxygen2 comment:
   #' @export
   ```

3. **Not handling sparse matrices**:
   ```r
   # Check for sparse matrix support:
   if (is(counts, "dgCMatrix")) {
       # Sparse matrix path
   } else {
       # Dense matrix path
   }
   ```

4. **Ignoring seed for reproducibility**:
   ```r
   # Always support seed parameter:
   .functionWithSeed <- function(..., seed = 12345) {
       withr::with_seed(seed, {
           .actualFunction(...)
       })
   }
   ```

5. **Not updating NEWS.md**:
   - Always document user-facing changes in `NEWS.md`
   - Follow existing format

### Understanding the Workflow

**For adding new features**:
1. Understand which model(s) are affected (C, G, CG, or shared)
2. Add function to appropriate R file
3. Add roxygen2 documentation
4. Add unit tests
5. Run `roxygen2::roxygenise()`
6. Run `devtools::test()`
7. Run `devtools::check()`
8. Update `NEWS.md`
9. Commit changes

**For fixing bugs**:
1. Write a failing test that reproduces the bug
2. Fix the bug in source code
3. Verify test now passes
4. Run full test suite
5. Update `NEWS.md`
6. Commit changes

**For performance optimization**:
1. Profile code to identify bottleneck
2. Implement optimization (consider C++ for inner loops)
3. Add benchmarking script or update existing ones
4. Verify results are identical (not just faster)
5. Document performance improvement in `BENCHMARKING.md`
6. Update `NEWS.md`
7. Commit changes

### Resources

**Internal Documentation**:
- Package Wiki: https://github.com/campbio/celda/wiki
- Developer's Coding Style: https://github.com/campbio/celda/wiki/Celda-Development-Coding-Style-Guide
- Robust and Efficient Code: https://github.com/campbio/celda/wiki/Celda-Development-Robust-and-Efficient-Code
- RStudio Configuration: https://github.com/campbio/celda/wiki/Celda-Development-Rstudio-configuration
- FAQ (Usage): https://github.com/campbio/celda/wiki/FAQ-on-how-to-use-celda
- FAQ (Development): https://github.com/campbio/celda/wiki/FAQ-on-package-development

**External Resources**:
- Bioconductor Packages: http://bioconductor.org/developers/package-guidelines/
- R Package Development: https://r-pkgs.org/
- Advanced R: https://adv-r.hadley.nz/
- Rcpp Documentation: https://www.rcpp.org/
- roxygen2 Documentation: https://roxygen2.r-lib.org/

### Quick Reference Commands

```r
# Development
devtools::load_all()        # Load package for testing
devtools::document()        # Generate documentation
devtools::test()            # Run tests
devtools::check()           # Full package check
devtools::build()           # Build package
devtools::install()         # Install locally

# C++ compilation
devtools::clean_dll()       # Clean compiled code
pkgbuild::compile_dll()     # Recompile C++ code

# Documentation
roxygen2::roxygenise()      # Generate .Rd files
pkgdown::build_site()       # Build website
devtools::build_vignettes() # Build vignettes

# Testing
testthat::test_file("tests/testthat/test-celda_C.R")
covr::package_coverage()    # Coverage report

# Benchmarking
Rscript benchmark_optimizations.R
Rscript benchmark_simple.R
```

---

## Version History

- **v1.0** (2025-11-13): Initial CLAUDE.md creation
  - Comprehensive codebase analysis
  - Documentation of architecture and conventions
  - Development workflow guidance

---

## Maintenance

This file should be updated when:
- New models are added
- Major architectural changes occur
- Development workflows change
- New conventions are established
- CI/CD pipeline is modified

**Last Updated**: 2025-11-13
**Maintained By**: AI Assistants working on celda
**Questions**: Open an issue at https://github.com/campbio/celda/issues
