# Module Decision Tree Implementation

## Overview

This implementation adds decision tree functionality for gene modules in celda. The decision tree traces the hierarchical splitting process from `recursiveSplitModule()` and provides tools to analyze and visualize module relationships.

## New Functions

### 1. `buildModuleDecisionTree(celdaList, minL = NULL, maxL = NULL)`

**Purpose**: Constructs a hierarchical decision tree structure from recursive module splitting results.

**Input**:
- `celdaList`: A celdaList object from `recursiveSplitModule()`
- `minL`, `maxL`: Optional range of L values to include

**Output**: A "moduleDecisionTree" object containing:
- `nodes`: Data frame with node information (L, module label, parent, gene count, log-likelihood)
- `edges`: Data frame with parent-child relationships
- `geneAssignments`: Gene-to-module mappings at each L
- `lValues`: L values included in the tree
- `celdaList`: Original celdaList object

**Example**:
```r
data(celdaGSim)
moduleSplit <- recursiveSplitModule(celdaGSim$counts, initialL = 3, maxL = 10)
modTree <- buildModuleDecisionTree(moduleSplit)
```

### 2. `plotModuleDecisionTree(x, labelModules = TRUE, labelGenes = TRUE, plotType = "dendrogram", ...)`

**Purpose**: Visualizes the module decision tree.

**Input**:
- `x`: A moduleDecisionTree object
- `labelModules`: Whether to label modules at terminal nodes
- `labelGenes`: Whether to show gene counts
- `plotType`: "dendrogram" or "network" (network requires igraph)

**Output**: A plot showing the hierarchical module relationships

**Example**:
```r
plotModuleDecisionTree(modTree)
```

### 3. `getModuleLineage(tree, L, module)`

**Purpose**: Traces the ancestry of a specific module back through lower L values.

**Input**:
- `tree`: A moduleDecisionTree object
- `L`: Number of modules where the target module exists
- `module`: Module label to trace

**Output**: Data frame showing the module's lineage from root to current

**Example**:
```r
# Trace module 5 at L=10 back to its origins
lineage <- getModuleLineage(modTree, L = 10, module = 5)
```

### 4. `getModuleDescendants(tree, L, module)`

**Purpose**: Finds all descendant modules that arose from splitting a specific module.

**Input**:
- `tree`: A moduleDecisionTree object
- `L`: Number of modules where the ancestral module exists
- `module`: Module label to find descendants for

**Output**: Data frame showing all descendant modules and their generation distance

**Example**:
```r
# Find all modules that descended from module 1 at L=3
descendants <- getModuleDescendants(modTree, L = 3, module = 1)
```

### 5. `print.moduleDecisionTree(x, ...)`

**Purpose**: Print method for moduleDecisionTree objects that provides a summary.

## Use Cases

1. **Understanding Module Relationships**: Trace how modules are related through the splitting process
2. **Module Stability Analysis**: Identify which modules remain stable vs. those that split frequently
3. **Optimal L Selection**: Visualize the tree to understand module granularity at different L values
4. **Gene Module Assignment**: Understand which genes stay together vs. separate during splitting
5. **Quality Control**: Identify anomalous splits or module behaviors

## Implementation Details

### Tree Construction Algorithm

1. Start with the initial L value and create root nodes for each module
2. For each subsequent L value:
   - Compare gene assignments between L and L-1
   - For each module at L, find its parent module at L-1 (the one sharing most genes)
   - Create edge from parent to child
   - Record node statistics (gene count, log-likelihood)

### Visualization

The dendrogram plot shows:
- Y-axis: L values (number of modules)
- X-axis: Module positions (distributed evenly)
- Node size/color: Number of genes in the module
- Edges: Parent-child relationships showing module splits

## File Location

Implementation: `/home/user/celda/R/moduleDecisionTree.R`

## Dependencies

Only uses base R packages:
- `graphics`: For plotting
- `grDevices`: For color palettes
- `methods`: For S4 class checking

## Next Steps

1. Run `roxygen2::roxygenise()` to generate documentation files
2. Add unit tests in `tests/testthat/`
3. Consider adding igraph integration for network visualization
4. Add examples to vignettes

## Example Workflow

```r
library(celda)

# Load data
data(celdaGSim)

# Perform recursive splitting
moduleSplit <- recursiveSplitModule(
  celdaGSim$counts,
  initialL = 5,
  maxL = 20
)

# Build decision tree
modTree <- buildModuleDecisionTree(moduleSplit)

# Print summary
print(modTree)

# Visualize tree
plotModuleDecisionTree(modTree)

# Analyze specific module
lineage <- getModuleLineage(modTree, L = 15, module = 7)
print(lineage)

# Find descendants of early module
descendants <- getModuleDescendants(modTree, L = 5, module = 2)
print(descendants)
```

## Notes

- The tree structure is built from the gene assignments at each L value
- Modules are linked based on maximum gene overlap between successive L values
- The implementation handles cases where modules may not split cleanly
- Log-likelihood values are stored at each node for quality assessment
