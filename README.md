# Crescendo
Crescendo for single-gene correction in single-cell RNA-sequencing data

Implementation of the Crescendo algorithm, which allows investigators to remove the effects of confounding factors by directly correcting gene expression count data.

## Installation

Install the current version of Crescendo from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("immunogenomics/crescendo")
```

# Usage/Demos

## Quick-start code

The following code uses Crescendo to correct gene expression from a public 10X dataset containing 3 batches (3pV1, 3pV2, and 5p)

``` r
library(crescendo)

# Load dataset with metadata, raw gene counts, and Harmony clusters (result from running the Harmony algorithm)
obj <- readRDS(system.file("extdata", "pbmc_4gene_obj.rds", package = "crescendo"))

# Set which genes to correct and parameters for coorrection
batch_var <- 'batch'
genes_use <- c('TRAC', 'MS4A1')
prop <- 0.05
min_cells <- 50

mc.cores <- NULL
lambda <- NULL
alpha <- 0

# Run Crescendo
corr_counts <- crescendo(
    Ycounts = obj$exprs_raw,
    meta_data = obj$meta_data,
    R = obj$R,
    batch_var = 'batch',
    genes_use = genes_use,
    prop = prop,
    min_cells = min_cells,
    lambda = lambda,
    alpha = alpha,
    mc.cores = mc.cores,
    return_obj = FALSE,
    verbose = FALSE
    # verbose = TRUE
)
```
