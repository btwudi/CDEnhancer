# CDEnhancer: cell-type-associated condition effects from bulk eRNA expression.
# Dependencies: limma (Bioconductor); statmod when robust = TRUE.
# Install once, if needed:
#   install.packages("BiocManager")
#   BiocManager::install(c("limma", "statmod"))
#
# Inputs:
#   eRNA: finite numeric features-by-samples matrix/data.frame, already filtered
#         and transformed as appropriate. Unique feature and sample names required.
#   frac: finite numeric samples-by-cell-types matrix/data.frame. Rows must match
#         colnames(eRNA) exactly. Fractions must be in [0,1], with row sums within
#         fraction_sum_tol of 1. No automatic reordering or normalization occurs.
#   type: numeric vector containing both 0 and 1, in expression sample order.
#         Positive effects correspond to condition 1 relative to condition 0.
#         If named, names must match expression sample IDs exactly.
#   confounding.factor: optional vector, matrix or data.frame in sample order.
#         Numeric/logical covariates and categorical factors/characters supported.
#         Explicit row names (or vector names) are checked against sample IDs.
#         Default data.frame row names imply positional matching: the caller must
#         ensure that unnamed covariates are in the correct sample order.
#   ref_idx: cell-type column index. NULL selects the greatest mean fraction.
#   robust, trend, adjust.method: passed to limma; original defaults retained.
#   fraction_sum_tol: validation tolerance only; fractions are never changed.
#
# Model: expression ~ condition * (K-1 cell fractions) + optional covariates.
# Cell-type effects: condition coefficient for the reference cell type;
# condition + corresponding interaction coefficient for each other cell type.
# With log-transformed bulk input, these are model-estimated effects on that
# scale, not necessarily conventional cell-type-specific log2 fold changes.
#
# Returns: a named list of limma topTable data.frames, one per cell type in
# original column order, with features in original row order. `logFC` stores the
# estimated effect. `adj.P.Val` is adjusted separately within each cell type
# across the supplied eRNAs (BH for adjust.method = "fdr" or "BH"). No significance
# or effect-size filtering is performed here. Statistical warnings are preserved.

CDEnhancer <- function(eRNA, frac, type, confounding.factor = NULL,
                       ref_idx = NULL, robust = TRUE, trend = TRUE,
                       adjust.method = "fdr", fraction_sum_tol = 1e-6) {
    fail <- function(message) stop(message, call. = FALSE)
    numeric_matrix <- function(x, label) {
        if (!(is.matrix(x) || is.data.frame(x)))
            fail(paste0(label, " must be a numeric matrix or data.frame."))
        if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1))))
            fail(paste0(label, " must contain only numeric columns."))
        x <- as.matrix(x)
        if (!is.numeric(x) || any(!is.finite(x)))
            fail(paste0(label, " must contain only finite numeric values."))
        x
    }
    valid_names <- function(x, label) {
        if (is.null(x) || anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x))
            fail(paste0(label, " must be present, nonempty and unique."))
    }
    check_ids <- function(ids, label) {
        if (!is.null(ids) && !identical(as.character(ids), sample_ids))
            fail(paste0(label, " must match colnames(eRNA) in the same order."))
    }

    if (!requireNamespace("limma", quietly = TRUE))
        fail("Package 'limma' is required. Install it with BiocManager::install('limma').")
    for (label in c("robust", "trend")) {
        value <- get(label)
        if (!is.logical(value) || length(value) != 1L || is.na(value))
            fail(paste0(label, " must be TRUE or FALSE."))
    }
    if (robust && !requireNamespace("statmod", quietly = TRUE))
        fail("Package 'statmod' is required for robust = TRUE.")
    if (!is.character(adjust.method) || length(adjust.method) != 1L ||
        is.na(adjust.method) || !adjust.method %in% stats::p.adjust.methods)
        fail("adjust.method must be a method supported by stats::p.adjust().")
    if (!is.numeric(fraction_sum_tol) || length(fraction_sum_tol) != 1L ||
        !is.finite(fraction_sum_tol) || fraction_sum_tol < 0)
        fail("fraction_sum_tol must be one finite nonnegative number.")

    # Default data.frame row numbers are not biological feature/sample IDs.
    if (is.data.frame(eRNA) && .row_names_info(eRNA, 1L) < 0L)
        fail("eRNA must have explicit unique feature row names.")
    if (is.data.frame(frac) && .row_names_info(frac, 1L) < 0L)
        fail("frac must have explicit sample row names matching colnames(eRNA).")
    eRNA <- numeric_matrix(eRNA, "eRNA")
    frac <- numeric_matrix(frac, "frac")
    if (nrow(eRNA) < 2L || ncol(eRNA) < 2L)
        fail("eRNA must contain at least two features and two samples.")
    if (ncol(frac) < 2L || nrow(frac) != ncol(eRNA))
        fail("frac must have at least two cell types and one row per expression sample.")
    valid_names(rownames(eRNA), "eRNA feature names")
    valid_names(colnames(eRNA), "eRNA sample names")
    sample_ids <- colnames(eRNA)
    valid_names(rownames(frac), "frac sample names")
    check_ids(rownames(frac), "frac row names")
    K <- ncol(frac)
    if (is.null(colnames(frac))) colnames(frac) <- paste0("CT", seq_len(K))
    valid_names(colnames(frac), "Cell-type names")
    if (any(frac < 0 | frac > 1)) fail("Cell fractions must lie between 0 and 1.")
    if (any(abs(rowSums(frac) - 1) > fraction_sum_tol))
        fail("Cell fractions must sum to 1 per sample within fraction_sum_tol. Check the inputs or explicitly set a tolerance appropriate for rounding.")
    if (!is.numeric(type) || !is.null(dim(type)) || length(type) != ncol(eRNA) ||
        any(!is.finite(type)) || !setequal(unique(type), c(0, 1)))
        fail("type must be a numeric 0/1 vector with both groups and one value per sample; factors are not accepted.")
    check_ids(names(type), "type names")
    if (is.null(ref_idx)) ref_idx <- which.max(colMeans(frac))
    if (!is.numeric(ref_idx) || length(ref_idx) != 1L || !is.finite(ref_idx) ||
        ref_idx != floor(ref_idx) || ref_idx < 1 || ref_idx > K)
        fail("ref_idx must be a single cell-type column index between 1 and ncol(frac).")

    use_idx <- setdiff(seq_len(K), ref_idx)
    frac_use <- as.data.frame(frac[, use_idx, drop = FALSE])
    colnames(frac_use) <- paste0("f", seq_along(use_idx))
    df <- data.frame(type = unname(type), frac_use, row.names = sample_ids)
    conf_names <- character(0)
    if (!is.null(confounding.factor)) {
        z <- confounding.factor
        if (is.data.frame(z)) {
            if (.row_names_info(z, 1L) > 0L) check_ids(rownames(z), "Covariate row names")
        } else if (is.matrix(z)) {
            check_ids(rownames(z), "Covariate row names")
        } else if (is.atomic(z) && is.null(dim(z))) {
            check_ids(names(z), "Covariate names")
            z <- data.frame(covariate = unname(z))
        } else fail("confounding.factor must be a vector, matrix or data.frame.")
        conf_df <- as.data.frame(z)
        if (nrow(conf_df) != ncol(eRNA) || ncol(conf_df) == 0L)
            fail("Covariates must have one row per sample and at least one column, or be NULL.")
        for (j in seq_len(ncol(conf_df))) {
            v <- conf_df[[j]]
            if (!(is.numeric(v) || is.factor(v) || is.character(v) || is.logical(v)) ||
                !is.null(dim(v))) fail("Covariates must be numeric, logical, character or factor vectors.")
            if (anyNA(v) || (is.numeric(v) && any(!is.finite(v))))
                fail("Covariates contain missing or nonfinite values. Resolve these consistently across all inputs before fitting.")
            if (length(unique(v)) < 2L) fail(paste0("Covariate '", names(conf_df)[j], "' is constant; remove it before fitting."))
            if (is.character(v)) v <- factor(v)
            if (is.factor(v)) v <- droplevels(v)
            conf_df[[j]] <- v
        }
        # Internal names prevent collisions with condition/fraction predictors.
        conf_names <- paste0("cov", seq_len(ncol(conf_df)))
        names(conf_df) <- conf_names
        rownames(conf_df) <- sample_ids
        df <- cbind(df, conf_df)
    }
    form <- paste0("~ type * (", paste(names(frac_use), collapse = " + "), ")",
                   if (length(conf_names)) paste0(" + ", paste(conf_names, collapse = " + ")) else "")
    design <- stats::model.matrix(stats::as.formula(form), data = df)
    if (nrow(design) != ncol(eRNA)) fail("The design matrix lost samples during construction.")
    design_rank <- qr(design)$rank
    if (design_rank < ncol(design))
        fail("Design matrix is rank deficient. Check constant/dependent fractions, covariates and condition confounding.")
    if (nrow(design) <= design_rank)
        fail("No residual degrees of freedom remain; reduce model complexity or add samples.")

    # Numeric contrasts avoid restrictions on cell-type labels and factor levels.
    contrasts <- matrix(0, nrow = ncol(design), ncol = K,
                        dimnames = list(colnames(design), paste0("eff", seq_len(K))))
    contrasts["type", ] <- 1
    for (j in seq_along(use_idx)) {
        term <- paste0("type:", names(frac_use)[j])
        if (!term %in% rownames(contrasts)) fail("Expected condition-fraction interaction is missing from the design.")
        contrasts[term, use_idx[j]] <- 1
    }
    fit <- limma::lmFit(eRNA, design)
    fit <- limma::contrasts.fit(fit, contrasts)
    fit <- limma::eBayes(fit, robust = robust, trend = trend)
    results <- lapply(seq_len(K), function(i) {
        tab <- limma::topTable(fit, coef = i, adjust.method = adjust.method,
                              number = Inf, sort.by = "none")
        tab[rownames(eRNA), , drop = FALSE]
    })
    names(results) <- colnames(frac)
    results
}
