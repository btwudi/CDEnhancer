#####
# use eRNAs to predict gene expression using TCGA data
# due to the limitation of the computing power, run one cancer type each time
# run on HPC
#####



#############################
# Robust nested-CV glmnet
# (safe against constant-y folds; numeric confounders)
#############################

suppressPackageStartupMessages({
        library(glmnet)
        library(caret)
        library(rslurm)
        library(rtracklayer)
})

# ---------- Helpers to avoid glmnet crashes on constant y ----------
is_constant <- function(v) {
        s <- stats::sd(v, na.rm = TRUE)
        isTRUE(!is.na(s) && s == 0)
}

safe_glmnet_fit <- function(x, y, lambda = NULL, alpha = 1) {
        x <- as.matrix(x)
        y <- as.numeric(y)
        if (is_constant(y)) {
                const <- unique(y)
                if (length(const) != 1L) const <- y[1]
                return(list(model = NULL, const = const, skipped = TRUE,
                            lambda = if (is.null(lambda)) 1 else lambda))
        }
        m <- glmnet::glmnet(x = x, y = y,
                            family = "gaussian",
                            alpha = alpha,
                            lambda = lambda,
                            standardize.response = FALSE)
        list(model = m, const = NULL, skipped = FALSE,
             lambda = if (is.null(lambda)) m$lambda else lambda)
}

safe_glmnet_predict <- function(fit, newx, s = NULL) {
        newx <- as.matrix(newx)
        if (isTRUE(fit$skipped)) {
                L <- if (is.null(s)) length(fit$lambda) else length(s)
                if (L < 1) L <- 1
                pred <- matrix(fit$const, nrow = nrow(newx), ncol = L)
                colnames(pred) <- if (is.null(s)) fit$lambda else s
                return(pred)
        }
        if (is.null(s)) {
                as.matrix(stats::predict(fit$model, newx = newx))
        } else {
                as.matrix(stats::predict(fit$model, newx = newx, s = s))
        }
}

# ---------- Core Nested CV ----------
nestedCV_glmnet <- function(train.ds = NULL, 
                            confounding.factor = NULL,
                            ncv_folds = c(5, 5),
                            inner_selection_percent = 0.2,
                            seed = NULL,
                            verbose = FALSE)
{
        if (!is.null(seed)) set.seed(seed)
        
        if (verbose) cat("\nRunning nested cross-validation...\n")
        ptm <- proc.time()
        
        # Folds for outer CV based on response (last column of train.ds)
        y_outer <- train.ds[, ncol(train.ds)]
        if (length(unique(y_outer)) <= ncv_folds[1]) {
                outer_folds <- caret::createFolds(factor(y_outer), ncv_folds[1], list = FALSE)
        } else {
                outer_folds <- caret::createFolds(y_outer, ncv_folds[1], list = FALSE)
        }
        
        
        intergrated_pred <- rep(NA_real_, nrow(train.ds))
        tune_params <- vector("list", ncv_folds[1])
        
        for (i in seq_len(ncv_folds[1])) {
                if (verbose) cat(sprintf("\n[Outer fold %d/%d]\n", i, ncv_folds[1]))
                
                # Indices for this outer split
                outer_idx <- which(outer_folds != i)
                test_idx  <- which(outer_folds == i)
                
                # Partition data
                outer_train.data <- as.matrix(train.ds[outer_idx, -ncol(train.ds), drop = FALSE])
                outer_test.data  <- as.matrix(train.ds[test_idx,  -ncol(train.ds), drop = FALSE])
                
                outer_train.cf <- as.matrix(confounding.factor[outer_idx, , drop = FALSE])
                outer_test.cf  <- as.matrix(confounding.factor[test_idx,  , drop = FALSE])
                
                outer_train.pheno <- as.numeric(train.ds[outer_idx, ncol(train.ds)])
                outer_test.pheno  <- as.numeric(train.ds[test_idx,  ncol(train.ds)])
                
                # Rank features by absolute Pearson correlation on outer train
                ranked_vars.train <- suppressWarnings(abs(stats::cor(outer_train.data, outer_train.pheno)))
                if (is.matrix(ranked_vars.train)) ranked_vars.train <- as.vector(ranked_vars.train)
                names(ranked_vars.train) <- colnames(train.ds)[-ncol(train.ds)]
                ranked_vars.train[is.na(ranked_vars.train)] <- 0
                
                wrapper.topN <- max(1, floor(inner_selection_percent * length(ranked_vars.train)))
                top_vars.train <- names(sort(ranked_vars.train, decreasing = TRUE))[seq_len(wrapper.topN)]
                
                # Get lambda grid from a safe outer fit on selected vars (train only)
                glmnet.obj <- safe_glmnet_fit(x = outer_train.data[, top_vars.train, drop = FALSE],
                                              y = outer_train.pheno)
                lambda <- glmnet.obj$lambda
                if (length(lambda) < 1) lambda <- 1
                
                # Prepare inner folds built on the outer training set only
                if (verbose) cat(sprintf(" Create [%d] inner folds of outer fold [%d]\n", ncv_folds[2], i))
                y_inner <- train.ds[outer_idx, ncol(train.ds)]
                if (length(unique(y_inner)) <= ncv_folds[2]) {
                        inner_folds <- caret::createFolds(factor(y_inner), ncv_folds[2], list = TRUE)
                } else {
                        inner_folds <- caret::createFolds(y_inner, ncv_folds[2], list = TRUE)
                }
                
                
                inner_test_acc_sum <- matrix(0, nrow = length(lambda), ncol = ncv_folds[2])
                inner_tune_params  <- vector("list", ncv_folds[2])
                inner_atts         <- vector("list", ncv_folds[2])
                
                for (j in seq_along(inner_folds)) {
                        # Map inner fold indices back to full-index space
                        outer_train_indices <- outer_idx
                        fold_test_local     <- inner_folds[[j]]             # positions within outer training set
                        inner_idx_test      <- outer_train_indices[fold_test_local]
                        inner_idx           <- setdiff(outer_train_indices, inner_idx_test)
                        
                        # Rank features on inner-train only
                        ranked_vars <- suppressWarnings(abs(stats::cor(train.ds[inner_idx, -ncol(train.ds), drop = FALSE],
                                                                       as.numeric(train.ds[inner_idx, ncol(train.ds)]))))
                        if (is.matrix(ranked_vars)) ranked_vars <- as.vector(ranked_vars)
                        names(ranked_vars) <- colnames(train.ds)[-ncol(train.ds)]
                        ranked_vars[is.na(ranked_vars)] <- 0
                        
                        wrapper.topN <- max(1, floor(inner_selection_percent * length(ranked_vars)))
                        top_vars <- names(sort(ranked_vars, decreasing = TRUE))[seq_len(wrapper.topN)]
                        
                        # Build inner design matrices
                        inner_trn.data <- as.matrix(train.ds[inner_idx, top_vars, drop = FALSE])
                        inner_tst.data <- as.matrix(train.ds[inner_idx_test, top_vars, drop = FALSE])
                        
                        inner_trn.cf <- as.matrix(confounding.factor[inner_idx, , drop = FALSE])
                        inner_tst.cf <- as.matrix(confounding.factor[inner_idx_test, , drop = FALSE])
                        
                        inner_trn.pheno <- as.numeric(train.ds[inner_idx, ncol(train.ds)])
                        inner_tst.pheno <- as.numeric(train.ds[inner_idx_test, ncol(train.ds)])
                        
                        # Fit safely on inner-train using the lambda grid
                        inner_train_model <- safe_glmnet_fit(
                                x = cbind(inner_trn.data, inner_trn.cf),
                                y = inner_trn.pheno,
                                lambda = lambda
                        )
                        
                        # Predict on inner-test for all lambdas
                        inner_test_pred <- safe_glmnet_predict(
                                inner_train_model,
                                newx = cbind(inner_tst.data, inner_tst.cf),
                                s = lambda
                        )
                        
                        # Compute correlation per lambda (Pearson, like your original inner metric)
                        if (is_constant(inner_tst.pheno)) {
                                inner_test_acc <- rep(NA_real_, ncol(inner_test_pred))
                        } else {
                                inner_test_acc <- apply(inner_test_pred, 2, function(p) stats::cor(inner_tst.pheno, p))
                        }
                        inner_test_acc[is.na(inner_test_acc)] <- -1
                        
                        inner_tune_params[[j]] <- data.frame(lambda = lambda, inner_test_acc = inner_test_acc)
                        inner_test_acc_sum[, j] <- inner_test_acc
                        inner_atts[[j]] <- top_vars
                }
                
                # Select lambda with the highest summed inner-test correlation
                acc.all <- apply(inner_test_acc_sum, 1, sum, na.rm = TRUE)
                best_idx <- which.max(acc.all)
                lambda.min <- lambda[best_idx]
                
                # Restrict outer train/test to the outer-selected features
                outer_train.data <- outer_train.data[, top_vars.train, drop = FALSE]
                outer_test.data  <- outer_test.data[,  top_vars.train, drop = FALSE]
                
                # Final outer fit & prediction (safe)
                outer_glmnet.model <- safe_glmnet_fit(
                        x = cbind(outer_train.data, outer_train.cf),
                        y = outer_train.pheno,
                        lambda = lambda.min
                )
                outer_test.pred <- as.numeric(safe_glmnet_predict(
                        outer_glmnet.model,
                        newx = cbind(outer_test.data, outer_test.cf),
                        s = lambda.min
                ))
                
                # Outer test accuracy (Spearman, as in your original code)
                outer_test_acc <- if (is_constant(outer_test.pheno) || is_constant(outer_test.pred)) {
                        NA_real_
                } else {
                        suppressWarnings(stats::cor(outer_test.pheno, outer_test.pred, method = "spearman", use = "pairwise.complete.obs"))
                }
                
                intergrated_pred[test_idx] <- outer_test.pred
                
                tune_params[[i]] <- list(accuracy = outer_test_acc,
                                         lambda = lambda.min,
                                         feature.weight = if (!isTRUE(outer_glmnet.model$skipped)) outer_glmnet.model$model$beta else NULL,
                                         feature = top_vars.train,
                                         prediction = outer_test.pred)
                names(tune_params)[i] <- paste0("fold", i)
                
                #print(i)
        }
        
        # Integrated (across all outer folds) Spearman correlation
        intergrated_acc <- if (is_constant(train.ds[, ncol(train.ds)]) || all(is.na(intergrated_pred)) || is_constant(intergrated_pred)) {
                NA_real_
        } else {
                suppressWarnings(stats::cor(train.ds[, ncol(train.ds)], intergrated_pred, method = "spearman", use = "pairwise.complete.obs"))
        }
        
        elapsed.time <- (proc.time() - ptm)[3]
        if (verbose) cat("nestedCV elapsed time (sec):", elapsed.time, "\n")
        
        list(intergrated_acc, tune_params)
}

# ---------- CV driver (serial) ----------
CV_prediction_main <- function(expression.feature, 
                               expression.target,
                               confounding.factor,
                               ncv_folds = c(5, 5),
                               inner_selection_percent = 0.2,
                               seed = NULL,
                               iteration = 1,
                               verbose = FALSE)
{
        if (ncol(expression.feature) != ncol(expression.target))
                stop("The sample sizes are not matched between the two datasets")
        
        average.acc <- matrix(0, nrow = iteration, ncol = nrow(expression.target))
        params <- vector("list", iteration)
        
        for (iter in seq_len(iteration)) {
                temp2 <- vector("list", nrow(expression.target))
                res <- rep(0, nrow(expression.target))
                
                for (i in seq_len(nrow(expression.target))) {
                        train.ds <- t(rbind(expression.feature, expression.target[i, ]))
                        temp <- nestedCV_glmnet(train.ds = train.ds,
                                                confounding.factor = confounding.factor,
                                                ncv_folds = ncv_folds,
                                                inner_selection_percent = inner_selection_percent,
                                                seed = iter,
                                                verbose = verbose)
                        res[i] <- temp[[1]]
                        temp2[[i]] <- temp[[2]]
                        message(sprintf("Finished the prediction for %dth gene!", i))
                }
                
                names(temp2) <- rownames(expression.target)
                average.acc[iter, ] <- res
                params[[iter]] <- temp2
                message(sprintf("Finished the %dth iteration!", iter))
        }
        
        names(params) <- paste0("interation", seq_len(iteration))
        out.res <- list(accuracy = average.acc, parameters = params)
        out.res
}

# ---------- SLURM Utilities (unchanged except for safety of globals) ----------
check_job_status <- function(job) {
        status <- "pending"
        queue_info <- tryCatch(system(paste("squeue -n", job$jobname), intern = TRUE), error = function(e) "")
        if (length(queue_info) <= 1) {
                status <- "done"
        } else {
                st_temp <- character()
                for (i in 2:length(queue_info)) {
                        st_temp <- c(st_temp, scan(text = queue_info[i], what = "character", quiet = TRUE)[5])
                }
                if (sum(st_temp == "R") > 0) status <- "running"
        }
        status
}

CV_prediction_main_parallel <- function(gene_index) {
        i <- gene_index
        average.acc <- rep(0, iteration)
        params <- vector("list", iteration)
        for (iter in seq_len(iteration)) {
                train.ds <- t(rbind(expression.feature, expression.target[i, ]))
                temp <- nestedCV_glmnet(train.ds = train.ds,
                                        confounding.factor = confounding.factor,
                                        ncv_folds = ncv_folds,
                                        inner_selection_percent = inner_selection_percent,
                                        seed = iter,
                                        verbose = verbose)
                average.acc[iter] <- temp[[1]]
                params[[iter]] <- temp[[2]]
        }
        names(params) <- paste0("interation", seq_len(iteration))
        list(average.acc, params, i)
}

CV.glmnet.main <- function(gene.index.vector, filename, time = '2:00:00', mem = '50g', nodes = 2) {
        vec <- data.frame(gene_index = gene.index.vector)
        obj <- c("expression.feature", "expression.target", "confounding.factor",
                 "ncv_folds", "inner_selection_percent", "seed", "iteration", "verbose", "nestedCV_glmnet",
                 "is_constant", "safe_glmnet_fit", "safe_glmnet_predict")  # ensure helpers are shipped
        sopt <- list(time = time, mem = mem)
        sjob <- rslurm::slurm_apply(CV_prediction_main_parallel, vec,
                                    jobname = sprintf("%s.prediction", filename),
                                    add_objects = obj, nodes = nodes, cpus_per_node = 2,
                                    slurm_options = sopt, submit = TRUE)
        while (check_job_status(sjob) != "done") Sys.sleep(120)
        message("Prediction is done!")
        res <- rslurm::get_slurm_out(sjob, outtype = 'table')
        saveRDS(res, sprintf("%s.res.rds", filename))
        
        acc_cor <- matrix(0, nrow = nrow(res), ncol = length(res$V1[[1]]))
        for (i in seq_len(nrow(res))) acc_cor[i, ] <- res$V1[[i]]
        index <- unlist(res$V3)
        acc_cor <- acc_cor[order(index), , drop = FALSE]
        rownames(acc_cor) <- rownames(expression.target)[index[order(index)]]
        acc_cor
}





#############################
# ---------- Data IO & Main ----------
#############################


dt = readRDS("~/TCGA_Pancancer_Gene_Enhancer15k_Metadata.rds")
erna = dt$enhancer
ge = dt$gene
meta = dt$metadata


gencode19 <- rtracklayer::import("~/gencode.v19.annotation_hg19.gtf")
gencode19 <- as.data.frame(gencode19)

# Parameters
ncv_folds <- c(5, 5)
inner_selection_percent <- 0.2
seed <- 123
iteration <- 10
verbose <- TRUE

# Storage for per-cancer-type results
prediction_output <- list()

# Iterate over cancer types
ctypes <- unique(meta$`cancer type abbreviation`)

run_one_cancer <- function(i) {
        
        ct <- ctypes[i]
        message(sprintf("\n==== Cancer type: %s (%d/%d) ====", ct, i, length(ctypes)))
        
        meta_sample <- meta[meta$`cancer type abbreviation` == ct, , drop = FALSE]
        erna_sub <- erna[, colnames(erna) %in% meta_sample$sample_new, drop = FALSE]
        ge_sub   <- ge[,   colnames(ge)   %in% meta_sample$sample_new, drop = FALSE]
        
        ## ---------------- Confounders ----------------
        age    <- as.factor(meta_sample$age)
        gender <- as.factor(meta_sample$gender)
        race   <- as.factor(meta_sample$race)
        group  <- as.factor(meta_sample$group)
        
        varying <- c(
                age    = nlevels(age)    >= 2,
                gender = nlevels(gender) >= 2,
                race   = nlevels(race)   >= 2,
                group  = nlevels(group)  >= 2
        )
        
        vars <- names(varying)[varying]
        
        if (length(vars) == 0) {
                confounding.factor <- matrix(numeric(0),
                                             nrow = nrow(meta_sample),
                                             ncol = 0)
        } else {
                form <- as.formula(paste("~ 0 +", paste(vars, collapse = " + ")))
                confounding.factor <- model.matrix(
                        form,
                        data = data.frame(age = age,
                                          gender = gender,
                                          race = race,
                                          group = group)
                )
        }
        
        rownames(confounding.factor) <- meta_sample$sample_new
        
        ## ---------------- Targets ----------------
        gencode.sub <- gencode19[gencode19$gene_type == "protein_coding", ]
        expression.target <- ge_sub[
                rownames(ge_sub) %in% gencode.sub$gene_name,
                ,
                drop = FALSE
        ]
        
        # sparsity filter
        count0 <- vapply(seq_len(nrow(expression.target)),
                         function(x) sum(expression.target[x, ] == log2(0.001)), # the gene value is log2(tpm + 0.001)
                         numeric(1))
        
        ctoff <- round(ncol(expression.target) * 0.2 + 1)
        expression.target <- expression.target[count0 < ctoff, , drop = FALSE]
        
        # variability filter
        sd.tar <- apply(expression.target, 1, sd)
        expression.target <- expression.target[sd.tar > 1, , drop = FALSE]
        idx <- sample(seq_len(ncol(expression.target)), 500) # subset the target because the server can only run 72 hours job at a time
        expression.target <- expression.target[, idx, drop = FALSE]
        
        
        ## ---------------- Features ----------------
        expression.feature <- erna_sub
        sd.fea <- apply(expression.feature, 1, sd)
        expression.feature <- expression.feature[sd.fea > 1, , drop = FALSE]
        expression.feature <- log2(expression.feature + 0.001)
        expression.feature <- t(scale(t(expression.feature)))
        
        ## ---------------- Align samples ----------------
        common_samples <- Reduce(
                intersect,
                list(colnames(expression.feature),
                     colnames(expression.target),
                     rownames(confounding.factor))
        )
        
        expression.feature <- expression.feature[, common_samples, drop = FALSE]
        expression.target  <- expression.target[,  common_samples, drop = FALSE]
        confounding.factor <- confounding.factor[common_samples, , drop = FALSE]
        
        ## ---------------- Run CV ----------------
        acc_cor <- CV_prediction_main(
                expression.feature,
                expression.target,
                confounding.factor,
                ncv_folds = ncv_folds,
                inner_selection_percent = inner_selection_percent,
                seed = seed,
                iteration = iteration,
                verbose = verbose
        )
        
        acc_cor_mat <- t(acc_cor$accuracy)
        rownames(acc_cor_mat) <- rownames(expression.target)
        colnames(acc_cor_mat) <- paste0("iteration", seq_len(nrow(acc_cor$accuracy)))
        
        ## ---------------- Save ----------------
        saveRDS(
                acc_cor_mat,
                file = sprintf(
                        "~/eRNA15k_predict_gene/%s_eRNA15k.pred.gene1k_iter=10.rds",
                        ct
                )
        )
        
        message(sprintf("Finished cancer type: %s", ct))
        
}


run_one_cancer(16)




















