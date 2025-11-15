# Celda Clustering Algorithm Improvements - Implementation Plan

**Date**: 2025-11-15
**Status**: In Progress
**Target Completion**: 2-3 weeks

---

## Overview

This document outlines the step-by-step implementation plan for the 8 major clustering algorithm improvements identified in `CLUSTERING_ALGORITHM_REVIEW.md`. Each improvement is broken down into concrete, testable tasks that can be implemented incrementally.

---

## Phase 1: Foundation & Quick Wins (Week 1)

### Task 1.1: Adaptive Feature Weighting During Clustering
**Priority**: HIGH
**Impact**: 10-15% better cluster purity
**Complexity**: MEDIUM
**Estimated Time**: 1-2 days

**Implementation Steps**:
1. Create `.calculateGeneWeights()` function to compute gene importance scores
   - Use between-cluster variance / total variance
   - Add optional marker gene boosting
   - Implement softmax normalization with clipping

2. Modify `.cGCalcGibbsProbY()` to accept gene weights
   - Add optional `geneWeights` parameter
   - Modify probability calculations to include weights
   - Ensure backward compatibility (weights = 1 when NULL)

3. Update `.celda_CG()` main loop
   - Add `featureReweighting = TRUE` parameter
   - Recalculate weights every 5 iterations
   - Apply weights to count matrix for next iteration

4. Add unit tests
   - Test weight calculation correctness
   - Test backward compatibility (weights = 1)
   - Test convergence with/without weighting

**Files to Modify**:
- `R/celda_CG.R` - Main loop
- `R/celda_G.R` - Gene clustering
- `R/feature_weights.R` - NEW FILE
- `tests/testthat/test-feature_weights.R` - NEW FILE

**Success Criteria**:
- ✓ All existing tests pass
- ✓ New tests for weighted clustering pass
- ✓ Benchmark shows improved cluster purity on test datasets
- ✓ Backward compatible (default behavior unchanged)

---

### Task 1.2: Marker-Guided Initialization
**Priority**: HIGH
**Impact**: 10-20% faster convergence
**Complexity**: MEDIUM
**Estimated Time**: 1-2 days

**Implementation Steps**:
1. Create `.initializeSplitZ_MarkerGuided()` function
   - Accept optional `markerGenes` list (named list of gene vectors)
   - Calculate marker scores for each cell
   - Use hierarchical clustering on marker scores
   - Return initial z assignments

2. Create `.calculateMarkerScores()` helper
   - For each cell, score against each marker set
   - Use mean/median expression of marker genes
   - Normalize scores across marker sets

3. Create `.initializeSplitZ_PriorClustering()` function
   - Accept optional `priorClustering` vector
   - Refine prior clustering to match target K
   - Use smart merging/splitting based on heterogeneity

4. Update `.initializeSplitZ()` dispatcher
   - Add `markerGenes` parameter
   - Add `priorClustering` parameter
   - Route to appropriate initialization method

5. Add unit tests
   - Test marker-guided initialization with known markers
   - Test prior clustering refinement
   - Test fallback to split initialization

**Files to Modify**:
- `R/initialize_clusters.R` - Add new initialization methods
- `R/marker_initialization.R` - NEW FILE
- `tests/testthat/test-initialize_cluster.R` - Extend existing tests

**Success Criteria**:
- ✓ Marker-guided initialization converges faster on PBMC data
- ✓ Prior clustering refinement produces similar/better results
- ✓ Backward compatible (NULL markers uses current method)

---

### Task 1.3: Adaptive K/L Subcluster Selection
**Priority**: MEDIUM
**Impact**: Better initialization
**Complexity**: LOW
**Estimated Time**: 0.5-1 day

**Implementation Steps**:
1. Create `.adaptiveKSubcluster()` function
   - Sample cells if dataset too large
   - Quick hierarchical clustering
   - Calculate silhouette scores
   - Return data-driven KSubcluster value

2. Create `.adaptiveLSubcluster()` equivalent for genes

3. Modify `.initializeSplitZ()` and `.initializeSplitY()`
   - Replace `ceiling(sqrt(K))` with adaptive calculation
   - Add `adaptiveSubclusters = TRUE` parameter
   - Allow manual override

4. Add unit tests
   - Test adaptive selection on datasets with varying structure
   - Test that it returns reasonable values

**Files to Modify**:
- `R/initialize_clusters.R`
- `tests/testthat/test-initialize_cluster.R`

**Success Criteria**:
- ✓ Adaptive method selects appropriate subclusters
- ✓ Backward compatible with fixed sqrt(K) method

---

## Phase 2: Core Algorithm Enhancements (Week 2)

### Task 2.1: Graph-Based Split Candidate Identification
**Priority**: HIGH
**Impact**: 20-30% better subcluster identification
**Complexity**: HIGH
**Estimated Time**: 2-3 days

**Implementation Steps**:
1. Create `.identifySplitCandidates_GraphBased()` function
   - Calculate statistical heterogeneity (current method)
   - Build kNN graph within each cluster
   - Calculate modularity/graph connectivity score
   - Detect bimodal gene expression patterns
   - Combine scores with weighting

2. Create `.calculateModularity()` helper
   - Accept kNN graph adjacency
   - Calculate modularity score
   - Return value 0-1

3. Create `.findBimodalGenes()` helper
   - Apply Hartigan's dip test
   - Return bimodal gene indices
   - Handle sparse data

4. Create `.detectSubstructure()` for correlation-based approach
   - Build correlation matrix
   - Detect community structure
   - Return substructure score

5. Modify `.cCSplitZ()` and `.cCGSplitZ()`
   - Add `useGraphBased = TRUE` parameter
   - Add optional `reducedDim` parameter for kNN
   - Route to new graph-based method

6. Add comprehensive unit tests
   - Test on clusters with known substructure
   - Test on homogeneous clusters
   - Test graph connectivity calculations
   - Test bimodality detection

**Files to Modify**:
- `R/split_clusters.R` - Extend splitting logic
- `R/graph_split.R` - NEW FILE
- `tests/testthat/test-split_clusters.R`

**Dependencies**:
- Consider adding optional dependency on `igraph` for graph metrics
- Or implement basic graph metrics directly

**Success Criteria**:
- ✓ Identifies subclusters that statistical method misses
- ✓ Doesn't over-split homogeneous clusters
- ✓ Computational overhead acceptable (<20% slowdown)

---

### Task 2.2: Advanced Convergence Detection
**Priority**: MEDIUM
**Impact**: Better resource utilization
**Complexity**: MEDIUM
**Estimated Time**: 1-2 days

**Implementation Steps**:
1. Create `.checkConvergence_Advanced()` function
   - Track log-likelihood history
   - Track cluster assignment history (z, y)
   - Calculate relative tolerance on LL
   - Calculate Adjusted Rand Index on assignments
   - Return convergence status + diagnostics

2. Modify `.celda_CG()` main loop
   - Store z/y history (last 5-10 iterations)
   - Call new convergence check
   - Add `convergenceMethod = "advanced"` parameter
   - Support legacy "simple" method for compatibility

3. Add early stopping for stable clusters
   - If ARI > 0.99 for 5 consecutive iterations, stop
   - Even if LL still slowly improving

4. Add unit tests
   - Test convergence detection on synthetic data
   - Test early stopping doesn't stop too early
   - Test backward compatibility

**Files to Modify**:
- `R/celda_CG.R` - Main loop
- `R/celda_C.R` - Cell clustering
- `R/celda_G.R` - Gene clustering
- `R/convergence.R` - NEW FILE
- `tests/testthat/test-convergence.R` - NEW FILE

**Success Criteria**:
- ✓ Detects convergence more accurately
- ✓ Saves iterations on simple datasets
- ✓ Doesn't stop prematurely on complex datasets

---

### Task 2.3: Internal Cluster Validation Metrics
**Priority**: MEDIUM
**Impact**: Better model selection
**Complexity**: MEDIUM
**Estimated Time**: 1-2 days

**Implementation Steps**:
1. Create `.calculateClusterMetrics()` function
   - Silhouette score (when feasible)
   - Calinski-Harabasz index
   - Davies-Bouldin index
   - Cluster size coefficient of variation
   - Module coherence (gene correlation within modules)
   - Return named list of metrics

2. Create helper functions:
   - `.calculateCH()` - Calinski-Harabasz
   - `.calculateDB()` - Davies-Bouldin
   - `.calculateModuleCoherence()` - Gene module quality
   - `.calculateSilhouette_Sampled()` - Sampled silhouette for large data

3. Modify chain selection in `.celda_CG()`
   - Calculate validation metrics for each chain
   - Add `selectionCriterion = c("logLik", "combined")` parameter
   - "combined" uses weighted score of LL + validation metrics
   - Store metrics in model object metadata

4. Create visualization function
   - `plotClusterMetrics()` to visualize validation scores
   - Compare metrics across different K/L values

5. Add unit tests
   - Test metric calculations on known good/bad clusterings
   - Test metric-based selection vs. LL-based selection
   - Test computational performance on large datasets

**Files to Modify**:
- `R/celda_CG.R` - Chain selection
- `R/cluster_validation.R` - NEW FILE
- `R/plot_validation.R` - NEW FILE
- `tests/testthat/test-cluster_validation.R` - NEW FILE

**Success Criteria**:
- ✓ Metrics correctly identify good vs. poor clusterings
- ✓ Combined selection performs better than LL alone on test data
- ✓ Computational overhead acceptable

---

## Phase 3: Advanced Features (Week 3)

### Task 3.1: Consensus Clustering Across Chains
**Priority**: MEDIUM
**Impact**: More robust results
**Complexity**: MEDIUM
**Estimated Time**: 1-2 days

**Implementation Steps**:
1. Create `.consensusClustering()` function
   - Accept list of chain results
   - Implement co-occurrence matrix method
   - Implement median assignment method
   - Calculate per-cell confidence scores
   - Return consensus clustering + diagnostics

2. Create `.buildCooccurrenceMatrix()` helper
   - Build cell-cell co-clustering matrix
   - Handle large datasets efficiently (sparse representation)
   - Normalize by number of chains

3. Modify `.celda_CG()` chain handling
   - Add `useConsensus = FALSE` parameter
   - When TRUE, run consensus instead of best chain selection
   - Store confidence scores in colData

4. Create visualization functions
   - `plotConsensusConfidence()` - Show cell-level confidence
   - `plotChainAgreement()` - Show agreement across chains

5. Add unit tests
   - Test consensus on chains with known agreement/disagreement
   - Test confidence score calculation
   - Test that consensus is at least as good as best chain

**Files to Modify**:
- `R/celda_CG.R` - Chain handling
- `R/consensus_clustering.R` - NEW FILE
- `R/plot_consensus.R` - NEW FILE
- `tests/testthat/test-consensus.R` - NEW FILE

**Success Criteria**:
- ✓ Consensus clustering more robust than single chain
- ✓ Confidence scores meaningfully reflect uncertainty
- ✓ Low-confidence cells identified correctly

---

### Task 3.2: Stochastic Mini-Batch for Large Datasets
**Priority**: LOW (unless needed)
**Impact**: 3-5x speedup for >50K cells
**Complexity**: HIGH
**Estimated Time**: 2-3 days

**Implementation Steps**:
1. Create `.celda_CG_MiniBatch()` variant
   - Add `useBatching = FALSE` parameter
   - Add `batchSize = NULL` parameter (auto-determined)
   - Add `adaptiveSampling = TRUE` parameter
   - Implement mini-batch sampling strategy

2. Create `.findUncertainCells()` helper
   - Identify cells with low maximum probability
   - Return indices of uncertain cells

3. Modify Gibbs sampling loop
   - Sample mini-batch of cells each iteration
   - Update only sampled cells
   - Periodically do full updates (every 10 iterations)
   - Track global log-likelihood

4. Ensure convergence with mini-batches
   - Modify convergence detection
   - Require full update convergence
   - Handle stochastic fluctuations

5. Add extensive unit tests
   - Test that mini-batch converges to same result
   - Test on datasets of varying sizes
   - Test adaptive sampling vs. uniform sampling

**Files to Modify**:
- `R/celda_CG.R` - Add mini-batch mode
- `R/minibatch.R` - NEW FILE
- `tests/testthat/test-minibatch.R` - NEW FILE

**Success Criteria**:
- ✓ Significant speedup on large datasets (>50K cells)
- ✓ Converges to similar results as full-batch
- ✓ Only activates when beneficial (auto-detection)

**Note**: This is lowest priority - implement only if needed for very large datasets.

---

### Task 3.3: scRNA-seq Biological Priors
**Priority**: LOW
**Impact**: 15-25% better biological interpretability
**Complexity**: HIGH
**Estimated Time**: 3-4 days

**Implementation Steps**:
1. Create `.adjustPriorsByQuality()` function
   - Accept QC metrics (nGenes, mtPercent, etc.)
   - Adjust alpha/beta based on cell quality
   - Return adjusted priors

2. Create `.cCCalcGibbsProbZ_ZeroInflated()` variant
   - Model zero-inflation explicitly
   - Use mixture of Dirichlet + zero component
   - Add `zeroInflationPrior` parameter

3. Create `.cCSplitZ_BatchAware()` variant
   - Accept batch labels
   - Penalize splits that separate batches
   - Allow splits with strong biological evidence

4. Add zero-inflation handling
   - Detect genes with excess zeros
   - Adjust probabilities accordingly
   - Add `handleZeroInflation = TRUE` parameter

5. Add unit tests
   - Test quality-based prior adjustment
   - Test zero-inflation handling
   - Test batch-aware splitting
   - Test on real scRNA-seq data with known biology

**Files to Modify**:
- `R/celda_CG.R` - Add biological priors
- `R/celda_C.R` - Batch-aware splitting
- `R/biological_priors.R` - NEW FILE
- `tests/testthat/test-biological_priors.R` - NEW FILE

**Success Criteria**:
- ✓ Improved clustering on dropout-heavy datasets
- ✓ Respects batch structure when appropriate
- ✓ Higher quality cells have more influence

**Note**: This requires most domain knowledge - implement last or get feedback from biologists.

---

## Phase 4: Testing & Benchmarking (Throughout + Final Week)

### Task 4.1: Comprehensive Unit Testing
**Priority**: CRITICAL
**Ongoing throughout implementation**

**Testing Strategy**:
1. **Correctness Tests**
   - Each new function has dedicated tests
   - Test edge cases (empty clusters, single cell, etc.)
   - Test numerical stability

2. **Backward Compatibility Tests**
   - All new parameters default to legacy behavior
   - Existing tests continue to pass
   - Results identical when features disabled

3. **Integration Tests**
   - Test combinations of improvements
   - Test full workflow with all features enabled
   - Test on multiple dataset types

4. **Regression Tests**
   - Save known-good results
   - Ensure future changes don't break these
   - Use snapshot testing where appropriate

**Files to Create/Modify**:
- All `tests/testthat/test-*.R` files
- `tests/testthat/helper-benchmark.R` - Benchmark utilities

---

### Task 4.2: Performance Benchmarking
**Priority**: HIGH
**Time**: Ongoing + 2-3 days final benchmarking

**Benchmark Datasets**:
1. **PBMC 3K** (small, well-characterized)
   - Test all improvements
   - Measure speedup and quality improvement
   - Compare with ground truth annotations

2. **PBMC 10K** (medium complexity)
   - Test scalability
   - Measure convergence speed
   - Test consensus clustering

3. **Simulated Data** (ground truth known)
   - Various cluster structures
   - Various noise levels
   - Test recovery of true clusters

4. **Mouse Brain 100K** (large, complex) - Optional
   - Test mini-batch if implemented
   - Test performance at scale
   - Memory profiling

**Metrics to Track**:
- **Speed**: Time to convergence (iterations × time/iteration)
- **Quality**: ARI, NMI vs. ground truth
- **Clustering Metrics**: Silhouette, CH index, DB index
- **Biological**: Marker gene enrichment, cell type purity
- **Memory**: Peak memory usage

**Create Benchmarking Scripts**:
- `benchmark_improvements.R` - Comprehensive benchmarking
- `benchmark_comparison.R` - Compare old vs. new
- `benchmark_ablation.R` - Test each improvement individually

**Success Criteria**:
- ✓ 10-20% faster convergence on average
- ✓ 10-20% better clustering quality on average
- ✓ No regression on any dataset
- ✓ Memory usage comparable or better

---

### Task 4.3: Documentation & Examples
**Priority**: HIGH
**Time**: 1-2 days

**Documentation Tasks**:
1. **Update Function Documentation**
   - Add roxygen2 docs for all new functions
   - Update existing docs with new parameters
   - Add examples demonstrating new features

2. **Create Vignettes**
   - "Advanced Clustering Features" vignette
   - Show marker-guided initialization
   - Show feature weighting
   - Show consensus clustering
   - Compare with/without improvements

3. **Update User Guide**
   - Add section on when to use which features
   - Add troubleshooting guide
   - Add performance tuning tips

4. **Create Developer Notes**
   - Document implementation details
   - Explain design decisions
   - Add architecture diagrams

**Files to Create/Modify**:
- `vignettes/advanced_clustering.Rmd` - NEW
- `man/*.Rd` - All updated via roxygen2
- `README.md` - Add feature highlights
- `NEWS.md` - Document all changes

---

## Implementation Schedule

### Week 1: Foundation (High-Impact Quick Wins)
- **Day 1-2**: Task 1.2 - Marker-guided initialization
- **Day 2-3**: Task 1.1 - Feature weighting
- **Day 3-4**: Task 1.3 - Adaptive subclusters
- **Day 4-5**: Task 2.2 - Advanced convergence detection
- **Ongoing**: Unit tests for all implemented features

### Week 2: Core Enhancements
- **Day 1-3**: Task 2.1 - Graph-based split heuristic (most complex)
- **Day 3-4**: Task 2.3 - Cluster validation metrics
- **Day 4-5**: Task 3.1 - Consensus clustering
- **Ongoing**: Integration testing

### Week 3: Polish & Optional Features
- **Day 1-2**: Task 4.2 - Comprehensive benchmarking
- **Day 2-3**: Task 4.3 - Documentation and vignettes
- **Day 3-4**: Bug fixes and refinements
- **Day 4-5**: Task 3.2 or 3.3 if time allows (optional)
- **Final**: Code review and merge preparation

---

## Success Metrics

### Quantitative Goals:
- ✓ 15-25% improvement in ARI vs. ground truth on test datasets
- ✓ 10-20% faster convergence (fewer iterations to convergence)
- ✓ 100% backward compatibility (all existing tests pass)
- ✓ <5% computational overhead when features disabled
- ✓ >90% test coverage for new code

### Qualitative Goals:
- ✓ Improved biological interpretability of clusters
- ✓ More robust clustering across multiple runs
- ✓ Better handling of difficult datasets (high dropout, batch effects)
- ✓ Easier to use for end users
- ✓ Well-documented and maintainable code

---

## Risk Mitigation

### Potential Risks:
1. **Performance Regression**: New features slow down algorithm
   - **Mitigation**: Make all features optional, profile carefully

2. **Backward Incompatibility**: Changes break existing code
   - **Mitigation**: Extensive compatibility testing, default to legacy behavior

3. **Increased Complexity**: Code becomes harder to maintain
   - **Mitigation**: Modular design, comprehensive documentation

4. **Biological Validity**: Improvements don't help real data
   - **Mitigation**: Validate on real datasets with known cell types

5. **Timeline Overrun**: Tasks take longer than expected
   - **Mitigation**: Prioritize high-impact items first, make later items optional

---

## Dependencies

### Required R Packages (already in celda):
- `Matrix` - Sparse matrices
- `Rcpp` - C++ integration
- `SingleCellExperiment` - Data structure

### Optional New Dependencies:
- `mclust` - For ARI calculation (may already be dep)
- `cluster` - For silhouette scores (may already be dep)
- `diptest` - For bimodality testing (NEW - consider vendoring)
- `igraph` - For graph metrics (NEW - optional, can implement basics ourselves)

### Testing Dependencies:
- `testthat` - Already used
- `covr` - Code coverage (already available)

---

## Deliverables

### Code:
- [ ] All tasks implemented and tested
- [ ] All existing tests pass
- [ ] New tests achieve >90% coverage
- [ ] Code reviewed and refined

### Documentation:
- [ ] All functions documented with roxygen2
- [ ] New vignette created
- [ ] README updated
- [ ] NEWS.md updated with all changes

### Benchmarks:
- [ ] Benchmark scripts created
- [ ] Results documented
- [ ] Performance comparison vs. baseline

### Release:
- [ ] All changes merged to development branch
- [ ] Version bumped appropriately
- [ ] Release notes prepared

---

## Next Steps

1. **Get stakeholder buy-in** on this plan
2. **Set up development branch** for these improvements
3. **Begin Phase 1** with marker-guided initialization
4. **Execute tasks** systematically using agent-based development
5. **Review and iterate** based on benchmarking results

---

## Notes

- This plan is aggressive but achievable with focused effort
- Priority should be on high-impact, lower-complexity items first
- Optional features (mini-batch, biological priors) can be deferred
- Continuous integration and testing throughout is critical
- Regular benchmarking will guide which improvements provide most value
