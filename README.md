# CDEnhancer

CDEnhancer identifies cell-type-associated differential enhancer RNA (eRNA)
expression between two conditions using bulk eRNA expression and estimated
cell-type fractions.

The framework fits a linear model containing condition, cell-type fractions,
condition–fraction interactions, and optional covariates. Models are fitted
using limma, with empirical Bayes moderation of the test statistics.

## Installation

CDEnhancer is currently provided as an R script.

Install the dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c("limma", "statmod"))
```

Download or clone this repository, then load the function:

```r
source("CDEnhancer.R")
```

The script was tested with R 4.4.1, limma 3.60.4, and statmod 1.5.0.

## Required inputs

| Input | Description |
|---|---|
| `eRNA` | Numeric matrix or data frame with eRNAs in rows and samples in columns. Unique feature and sample names are required. |
| `frac` | Numeric matrix or data frame with samples in rows and cell types in columns. |
| `type` | Numeric vector coded as 0 and 1, containing both conditions. |
| `confounding.factor` | Optional sample-level covariates, such as age, sex, or tumor grade. |

### Input preparation

- Supply filtered and appropriately normalized/transformed eRNA expression.
  The function does not perform these preprocessing steps.
- Each fraction must lie between 0 and 1, and fractions must sum to 1 within
  each sample, allowing for the specified numerical tolerance.
- Row names of `frac` must exactly match column names of `eRNA`, including order.
- Condition labels and covariates must follow the same sample order.
- Missing and nonfinite values are not accepted.
- Use factors or character vectors for categorical covariates.
- Remove constant or linearly dependent predictors before fitting.

Named condition vectors and explicitly named covariate rows are checked
against expression sample IDs. Unnamed covariates are matched by position;
their order must be verified by the user.

## Example

The following synthetic example demonstrates the input format and function call.
It is not a performance benchmark.

```r
source("CDEnhancer.R")

set.seed(123)

n_samples <- 80
n_eRNAs <- 300
sample_ids <- paste0("Sample_", seq_len(n_samples))

# Simulated fractions for three cell populations
frac <- matrix(rexp(n_samples * 3), nrow = n_samples)
frac <- frac / rowSums(frac)

rownames(frac) <- sample_ids
colnames(frac) <- c("Malignant", "Immune", "Stromal")

# Condition 0 is the comparison baseline
type <- rep(c(0, 1), each = n_samples / 2)
names(type) <- sample_ids

# Simulated expression values on an analysis scale
eRNA <- matrix(
    rnorm(n_eRNAs * n_samples),
    nrow = n_eRNAs,
    dimnames = list(
        paste0("eRNA_", seq_len(n_eRNAs)),
        sample_ids
    )
)

# Add a simulated malignant-cell-associated condition effect
eRNA[1:30, ] <- sweep(
    eRNA[1:30, , drop = FALSE],
    2,
    2 * type * frac[, "Malignant"],
    "+"
)

covariates <- data.frame(
    age = rnorm(n_samples, mean = 55, sd = 10),
    sex = factor(rep(c("Female", "Male"), length.out = n_samples)),
    row.names = sample_ids
)

results <- CDEnhancer(
    eRNA = eRNA,
    frac = frac,
    type = type,
    confounding.factor = covariates,
    ref_idx = 3,
    robust = TRUE,
    trend = FALSE,
    adjust.method = "BH"
)

head(results$Malignant)
```

## Main parameters

| Parameter | Default | Description |
|---|---|---|
| `ref_idx` | `NULL` | Reference cell-type column. By default, the population with the largest mean fraction is selected. |
| `robust` | `TRUE` | Use robust empirical Bayes estimation in limma. |
| `trend` | `TRUE` | Allow an abundance-dependent trend in the prior variance. |
| `adjust.method` | `"fdr"` | Multiple-testing adjustment method; `"fdr"` is an alias for `"BH"`. |
| `fraction_sum_tol` | `1e-6` | Allowed deviation of each fraction row sum from 1. |

The reference cell type is distinct from the reference condition.
The reference condition is always `type = 0`.

## Statistical model

For each eRNA, CDEnhancer fits:

```text
expression ~ condition * (K − 1 cell-type fractions) + covariates
```

One cell-type fraction is omitted because the fractions sum to one.

The condition effect is estimated as:

- Reference cell type: the condition coefficient.
- Other cell types: the condition coefficient plus the corresponding
  condition–fraction interaction coefficient.

## Output

The function returns a named list containing one result table per cell type.
Cell types and eRNAs retain their input order.

Key columns include:

| Column | Description |
|---|---|
| `logFC` | Model-estimated condition effect on the input expression scale. |
| `AveExpr` | Mean expression on the input scale. |
| `t` | Moderated t statistic. |
| `P.Value` | Two-sided P value for the cell-type-associated condition effect. |
| `adj.P.Val` | Multiple-testing-adjusted P value. |
| `B` | limma log-odds statistic for differential expression. |

Positive effects indicate higher inferred expression in condition 1 relative
to condition 0; negative effects indicate the reverse.

Despite the limma column name `logFC`, the estimate is not automatically a
conventional cell-type-specific log2 fold change. Its interpretation depends
on the input transformation and model assumptions.

P values are adjusted across eRNAs separately within each cell type.
The function does not apply significance or effect-size thresholds.

## Interpretation and limitations

Cell-type associations are inferred from bulk expression and supplied
fractions; they are not direct measurements of expression in isolated cells.
Results depend on fraction accuracy, sample size, covariate specification,
and the suitability of the linear model.

Interactions fitted to log-transformed bulk expression should be interpreted
on that analysis scale. They do not constitute a physical linear mixture
model of untransformed transcript abundance.

Nearest-gene assignment, pathway enrichment, and functional validation are
separate downstream analyses and are not performed by this function.

## Citation

[Add the CDEnhancer manuscript or preprint citation when available.]

Please also cite limma and the methods used to obtain eRNA expression and
cell-type fractions.

## License

[Specify the repository license and include the corresponding LICENSE file.]

## Contact

[Add the maintainer name and contact information, or direct users to GitHub Issues.]