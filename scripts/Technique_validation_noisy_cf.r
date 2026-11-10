# select N number of eRNAs to alter basing on SNR
# num is better to be even
# direction is "up", "down" or "bi" (bi means some are up, some are down)
# run on a HPC

library(limma)
library(pROC)
library(PRROC)
library(MCMCpack)
library(EpiDISH)
library(TOAST)

cf_generation <- function(
                sample_num = 284,
                var_ctrl   = 1,
                r          = 0.9,   # target correlation 
                out_dir    = NULL
) {
        set.seed(123)
        
        means <- c(0.30, 0.20, 0.10, 0.40)
        alpha <- means * var_ctrl
        true_mat <- rdirichlet(sample_num, alpha)
        
        if (r < 1) {
                achieved_r_for_sd <- function(sd_try) {
                        set.seed(999)  # <- make search stable (minimal but important)
                        noise <- matrix(
                                rnorm(sample_num * length(means), mean = 0, sd = sd_try),
                                nrow = sample_num, ncol = length(means)
                        )
                        noisy <- abs(true_mat + noise)
                        row_sums <- rowSums(noisy)
                        row_sums[row_sums == 0] <- 1
                        pred_tmp <- sweep(noisy, 1, row_sums, "/")
                        cor(as.vector(true_mat), as.vector(pred_tmp))
                }
                
                lo <- 0
                hi <- 1
                while (achieved_r_for_sd(hi) > r && hi < 100) hi <- hi * 2
                
                for (iter in 1:25) {
                        mid <- (lo + hi) / 2
                        rr  <- achieved_r_for_sd(mid)
                        if (rr > r) lo <- mid else hi <- mid
                }
                noise_sd <- hi
        } else {
                noise_sd <- 0
        }
        
        if (noise_sd > 0) {
                set.seed(456)  # <- stable final generation
                noise <- matrix(
                        rnorm(sample_num * length(means), mean = 0, sd = noise_sd),
                        nrow = sample_num,
                        ncol = length(means)
                )
                
                noisy <- abs(true_mat + noise)
                row_sums <- rowSums(noisy)
                row_sums[row_sums == 0] <- 1
                pred_mat <- sweep(noisy, 1, row_sums, "/")
        } else {
                pred_mat <- true_mat
        }
        
        colnames(true_mat) <- colnames(pred_mat) <- c("V1", "V2", "V3", "V4")
        
        ls <- list(
                true       = as.data.frame(true_mat),
                pred       = as.data.frame(pred_mat),
                target_r   = r,
                noise_sd   = noise_sd,
                achieved_r = cor(as.vector(true_mat), as.vector(pred_mat))
        )
        
        # only save if out_dir is provided
        if (!is.null(out_dir)) {
                if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
                file_name <- file.path(
                        out_dir,
                        paste0("generated_cf_and_mimic_predicted_cf_var_ctrl_", var_ctrl,
                               "_target_r_", r, ".RDS")
                )
                saveRDS(ls, file = file_name)
        }
        
        
        return(ls)
}







altered_eRNA_table <- function(data.list, ct, test, num, snr, direction)
{
        dt = list()
        altered.eRNA = list()
        for (i in 1:length(data.list))
        {
                row_sd =  apply(data.list[[i]], 1, sd)
                pos.sd.eRNA = data.list[[i]][row_sd >0, ]
                select.eRNA = sample(x = rownames(pos.sd.eRNA), size = num)
                temp = pos.sd.eRNA[rownames(pos.sd.eRNA) %in% select.eRNA, ]
                
                temp.ct = temp[, colnames(temp) %in% ct]
                temp.test = temp[, colnames(temp) %in% test]
                
                
                if (direction == "up") {
                        for (j in 1:nrow(temp.test))
                        {
                                N = snr * ((sd(temp.ct[j,]) + sd(temp.test[j,]))/2)
                                temp.test[j,] = N + temp.test[j,] + mean(temp.ct[j,]) - mean(temp.test[j,])
                                
                                if (sum(temp.test[j,] < 0) >0 )
                                {
                                        temp.test[j,] = temp.test[j,] + abs(min(temp.test[j,]))
                                        temp.ct[j,] = temp.ct[j,] + abs(min(temp.test[j,]))
                                }
                        }
                        altered.eRNA[[i]] = cbind.data.frame(altered_eRNAs = select.eRNA,
                                                             direction = rep("up", length(select.eRNA)))
                }
                
                
                
                if (direction == "down") {
                        for (j in 1:nrow(temp.test))
                        {
                                N = snr * ((sd(temp.ct[j,]) + sd(temp.test[j,]))/2)
                                temp.ct[j,] = N + temp.ct[j,] + mean(temp.test[j,]) - mean(temp.ct[j,])
                                
                                if (sum(temp.ct[j,] < 0) > 0)
                                {
                                        temp.ct[j,] = temp.ct[j,] + abs(min(temp.ct[j,]))
                                        temp.test[j,] = temp.test[j,] + abs(min(temp.ct[j,]))
                                }
                        }
                        altered.eRNA[[i]] = cbind.data.frame(altered_eRNAs = select.eRNA,
                                                             direction = rep("down", length(select.eRNA)))
                }
                
                
                if (direction == "bi") {
                        
                        select.up = sample(x = select.eRNA, size = (num/2))
                        select.down = select.eRNA[!select.eRNA %in% select.up]
                        
                        temp.up = temp.test[rownames(temp.test) %in% select.up,]
                        temp.ct.up = temp.ct[rownames(temp.ct) %in% select.up,]
                        
                        temp.down = temp.test[rownames(temp.test) %in% select.down,]
                        temp.ct.down = temp.ct[rownames(temp.ct) %in% select.down,]
                        
                        for (j in 1:(num/2))
                        {
                                # half up
                                N.up = snr * ((sd(temp.up[j,]) + sd(temp.ct.up[j,]))/2)
                                temp.up[j,] = N.up + temp.up[j,] + mean(temp.ct.up[j,]) - mean(temp.up[j,])
                                
                                if (sum(temp.up[j,] < 0) > 0)
                                {
                                        temp.up[j,] = temp.up[j,] + abs(min(temp.up[j,]))
                                        temp.ct.up[j,] = temp.ct.up[j,] + abs(min(temp.up[j,]))
                                }
                                
                                # half down
                                N.down = snr * ((sd(temp.down[j,]) + sd(temp.ct.down[j,]))/2)
                                temp.ct.down[j,] = N.down + temp.ct.down[j,] + mean(temp.down[j,]) - mean(temp.ct.down[j,])
                                
                                if (sum(temp.ct.down[j,] < 0) > 0)
                                {
                                        temp.ct.down[j,] = temp.ct.down[j,] + abs(min(temp.ct.down[j,]))
                                        temp.down[j,] = temp.down[j,] + abs(min(temp.ct.down[j,]))
                                }
                                
                        }
                        altered.eRNA[[i]] = cbind.data.frame(altered_eRNAs = c(select.up, select.down),
                                                             direction = c(rep("up", length(select.up)), 
                                                                           rep("down", length(select.down))))
                        ind.alt = match(altered.eRNA[[i]]$altered_eRNAs, select.eRNA)
                        altered.eRNA[[i]] = altered.eRNA[[i]][ind.alt, ]
                        
                        temp.ct.adj = rbind(temp.ct.up, temp.ct.down)
                        temp.test.adj = rbind(temp.up, temp.down)
                        
                        ind.ct = match(rownames(temp.ct), rownames(temp.ct.adj))
                        ind.test = match(rownames(temp.test), rownames(temp.test.adj))
                        
                        temp.ct.adj = temp.ct.adj[ind.ct, ]
                        temp.test.adj = temp.test.adj[ind.test, ]
                        
                        temp.ct = temp.ct.adj
                        temp.test = temp.test.adj
                }
                
                temp.adjusted = cbind(temp.ct, temp.test)
                ind.col = match(colnames(data.list[[i]]), colnames(temp.adjusted))
                temp.adjusted = temp.adjusted[, ind.col]
                
                ind.row = match(rownames(temp.adjusted), rownames(data.list[[i]]))
                data.list[[i]][ind.row, ] = temp.adjusted
                
                dt[[i]] = data.list[[i]]  
        }
        names(dt) = names(data.list)
        names(altered.eRNA) = names(data.list)
        
        altered.table = list(table.list = dt, altered.eRNA = altered.eRNA)
}



#########################
# get stats for CellDMC #
#########################
get.stats_CellDMC <- function(rec, altered.eRNA, num, logFC_cutoff, qv_cutoff, new.table, direction)
{
        if (direction == "up") {
                for(i in seq_along(rec))
                {
                        rec_temp = rec[[i]]
                        colnames(rec_temp)[c(1,4)] = c("logFC", "P.Value")
                        rec_temp$one_side_qv = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC <= 0,]$one_side_qv = 1-rec_temp[rec_temp$logFC <= 0,]$one_side_qv
                        rec_temp$qv = p.adjust(rec_temp$one_side_qv, method = "fdr")
                        rec[[i]] = rec_temp
                }
                
        } else if (direction == "down") {
                for(i in seq_along(rec))
                {
                        rec_temp = rec[[i]]
                        colnames(rec_temp)[c(1,4)] = c("logFC", "P.Value")
                        rec_temp$one_side_qv = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC >= 0,]$one_side_qv = 1-rec_temp[rec_temp$logFC >= 0,]$one_side_qv
                        rec_temp$qv = p.adjust(rec_temp$one_side_qv, method = "fdr")
                        rec[[i]] = rec_temp
                }
                
        } else if (direction == "bi") {
                for(i in seq_along(rec))
                {
                        rec_temp = rec[[i]]
                        colnames(rec_temp)[c(1,4)] = c("logFC", "P.Value")
                        
                        # calculate qv as in upregulated direction
                        rec_temp$one_side_qv.up = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC <= 0,]$one_side_qv.up = 1-rec_temp[rec_temp$logFC <= 0,]$one_side_qv.up
                        rec_temp$qv.up = p.adjust(rec_temp$one_side_qv.up, method = "fdr")
                        
                        # calculate qv as in downregulated direction
                        rec_temp$one_side_qv.down = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC >= 0,]$one_side_qv.down = 1-rec_temp[rec_temp$logFC >= 0,]$one_side_qv.down
                        rec_temp$qv.down = p.adjust(rec_temp$one_side_qv.down, method = "fdr")
                        
                        rec[[i]] = rec_temp
                }
                
        }
        
        
        if (direction == "bi") {
                
                qv.up = cbind(rec[[1]]$qv.up, 
                              rec[[2]]$qv.up,
                              rec[[3]]$qv.up,
                              rec[[4]]$qv.up)
                rownames(qv.up) = rownames(rec[[1]])
                colnames(qv.up) = names(rec)
                
                qv.down = cbind(rec[[1]]$qv.down, 
                                rec[[2]]$qv.down,
                                rec[[3]]$qv.down,
                                rec[[4]]$qv.down)
                rownames(qv.down) = rownames(rec[[1]])
                colnames(qv.down) = names(rec)
                
                estimate = cbind(rec[[1]]$logFC, 
                                 rec[[2]]$logFC,
                                 rec[[3]]$logFC,
                                 rec[[4]]$logFC)
                rownames(estimate) = rownames(rec[[1]])
                colnames(estimate) = names(rec)
                
                dif.eRNA.up = NULL
                dif.eRNA.down = NULL
                for (j in 1: ncol(estimate))
                {
                        temp.up = cbind.data.frame(estimate = estimate[, j], q_value.up = qv.up[, j])
                        rownames(temp.up) = rownames(rec[[j]])
                        temp.up = temp.up[complete.cases(temp.up),]
                        sig.up = temp.up[temp.up$estimate > logFC_cutoff & temp.up$q_value.up < qv_cutoff, ]
                        
                        temp.down = cbind.data.frame(estimate = estimate[, j], q_value.down = qv.down[, j])
                        rownames(temp.down) = rownames(rec[[j]])
                        temp.down = temp.down[complete.cases(temp.down),]
                        sig.down = temp.down[temp.down$estimate < -logFC_cutoff & temp.down$q_value.down < qv_cutoff, ]
                        
                        temp.up = cbind.data.frame(eRNA = rownames(sig.up),
                                                   cell_type = rep(colnames(qv.up)[j], nrow(sig.up)),
                                                   direction = rep("up", nrow(sig.up)),
                                                   qv = sig.up$q_value.up,
                                                   logFC = sig.up$estimate)
                        temp.down = cbind.data.frame(eRNA = rownames(sig.down),
                                                     cell_type = rep(colnames(qv.down)[j], nrow(sig.down)),
                                                     direction = rep("down", nrow(sig.down)),
                                                     qv = sig.down$q_value.down,
                                                     logFC = sig.down$estimate)
                        
        
                        up.eRNA.ind = (temp.up$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "up",]$altered_eRNAs)
                        down.eRNA.ind = (temp.down$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "down",]$altered_eRNAs)
                        
                        
                        
                        # Check if temp.up has rows
                        if (nrow(temp.up) > 0) {
                                temp.up$group = paste("false", temp.up$direction, sep = " ")
                                temp.up$group[up.eRNA.ind] = "true up"
                                dif.eRNA.up = rbind.data.frame(dif.eRNA.up, temp.up)
                        }
                        
                        # Check if temp.down has rows
                        if (nrow(temp.down) > 0) {
                                temp.down$group = paste("false", temp.down$direction, sep = " ")
                                temp.down$group[down.eRNA.ind] = "true down"
                                dif.eRNA.down = rbind.data.frame(dif.eRNA.down, temp.down)
                        }
                        
                        
                }
                
                # Combine the results if any rows exist
                if (exists("dif.eRNA.up") & exists("dif.eRNA.down")) {
                        dif.eRNA = rbind(dif.eRNA.up, dif.eRNA.down)
                } else if (exists("dif.eRNA.up")) {
                        dif.eRNA = dif.eRNA.up
                } else if (exists("dif.eRNA.down")) {
                        dif.eRNA = dif.eRNA.down
                } else {
                        dif.eRNA = data.frame()  # Empty dataframe if no rows exist
                }
                
                
                
                
        } else {
                qv = cbind(rec[[1]]$qv, 
                           rec[[2]]$qv,
                           rec[[3]]$qv,
                           rec[[4]]$qv)
                rownames(qv) = rownames(rec[[1]])
                colnames(qv) = names(rec)
                
                estimate = cbind(rec[[1]]$logFC, 
                                 rec[[2]]$logFC,
                                 rec[[3]]$logFC,
                                 rec[[4]]$logFC)
                rownames(estimate) = rownames(rec[[1]])
                colnames(estimate) = names(rec)
                
                
                
                dif.eRNA.up = NULL
                dif.eRNA.down = NULL
                for (j in 1: ncol(estimate))
                {
                        temp = cbind.data.frame(estimate = estimate[, j], q_value = qv[, j])
                        rownames(temp) = rownames(rec[[j]])
                        temp = temp[complete.cases(temp),]
                        sig.up = temp[temp$estimate > logFC_cutoff & temp$q_value < qv_cutoff, ]
                        sig.down = temp[temp$estimate < -logFC_cutoff & temp$q_value < qv_cutoff, ]
                        
                        temp.up = cbind.data.frame(eRNA = rownames(sig.up),
                                                   cell_type = rep(names(rec)[j], nrow(sig.up)),
                                                   direction = rep("up", nrow(sig.up)),
                                                   qv = sig.up$q_value,
                                                   logFC = sig.up$estimate)
                        temp.down = cbind.data.frame(eRNA = rownames(sig.down),
                                                     cell_type = rep(names(rec)[j], nrow(sig.down)),
                                                     direction = rep("down", nrow(sig.down)),
                                                     qv = sig.down$q_value,
                                                     logFC = sig.down$estimate)
                        
                        up.eRNA.ind = (temp.up$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "up",]$altered_eRNAs)
                        down.eRNA.ind = (temp.down$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "down",]$altered_eRNAs)
                        
                        
                        
                        # Check if temp.up has rows
                        if (nrow(temp.up) > 0) {
                                temp.up$group = paste("false", temp.up$direction, sep = " ")
                                temp.up$group[up.eRNA.ind] = "true up"
                                dif.eRNA.up = rbind.data.frame(dif.eRNA.up, temp.up)
                        }
                        
                        # Check if temp.down has rows
                        if (nrow(temp.down) > 0) {
                                temp.down$group = paste("false", temp.down$direction, sep = " ")
                                temp.down$group[down.eRNA.ind] = "true down"
                                dif.eRNA.down = rbind.data.frame(dif.eRNA.down, temp.down)
                        }
                        
                        
                }
                
                # Combine the results if any rows exist
                if (exists("dif.eRNA.up") & exists("dif.eRNA.down")) {
                        dif.eRNA = rbind(dif.eRNA.up, dif.eRNA.down)
                } else if (exists("dif.eRNA.up")) {
                        dif.eRNA = dif.eRNA.up
                } else if (exists("dif.eRNA.down")) {
                        dif.eRNA = dif.eRNA.down
                } else {
                        dif.eRNA = data.frame()  # Empty dataframe if no rows exist
                }
        }
        
        
        sens.single.run = NULL
        speci.single.run = NULL
        preci.single.run = NULL
        auc_results = list()
        pr_auc_results = list()
        
        for (n in seq_along(unique(dif.eRNA$cell_type)))
        {
                sub = dif.eRNA[dif.eRNA$cell_type == unique(dif.eRNA$cell_type)[n], ]
                sens = nrow(sub[sub$group == "true up" | sub$group == "true down", ]) / num
                sens.single.run = c(sens.single.run, sens)
                
                speci = (nrow(new.table[[n]]) - num - nrow(sub[sub$group == "false up" | sub$group == "false down", ])) / (nrow(new.table[[n]]) - num)
                speci.single.run = c(speci.single.run, speci)
                
                preci = nrow(sub[sub$group == "true up" | sub$group == "true down", ]) / (nrow(sub[sub$group == "true up" | sub$group == "true down", ]) + nrow(sub[sub$group == "false up" | sub$group == "false down", ]))
                preci.single.run = c(preci.single.run, preci)
                
                
                
                
                # AUC
                
                if (direction == "bi") {
                        
                        # ROC_AUC
                        predictions.up = qv.up[, n]
                        names(predictions.up) = rownames(rec[[n]])
                        
                        predictions.down = qv.down[, n]
                        names(predictions.down) = rownames(rec[[n]])
                        
                        temp = new.table[[n]]
                        
                        ind.temp.up = (rownames(temp) %in% names(predictions.up))
                        temp.up = temp[ind.temp.up, ]
                        
                        ind.temp.down = (rownames(temp) %in% names(predictions.down))
                        temp.down = temp[ind.temp.down, ]
                        
                        actual.up = rep(0, nrow(temp.up))
                        ind.up = (rownames(temp.up) %in% altered.eRNA[[n]][altered.eRNA[[n]]$direction == "up",]$altered_eRNAs)
                        actual.up[ind.up] = 1
                        
                        roc_result.up = roc(actual.up, predictions.up, direction = ">", quiet = TRUE)
                        auc_results[["AUC_up"]][n] = auc(roc_result.up)
                        
                        actual.down = rep(0, nrow(temp.down))
                        ind.down = (rownames(temp.down) %in% altered.eRNA[[n]][altered.eRNA[[n]]$direction == "down",]$altered_eRNAs)
                        actual.down[ind.down] = 1
                        
                        roc_result.down = roc(actual.down, predictions.down, direction = ">", quiet = TRUE)
                        auc_results[["AUC_down"]][n] = auc(roc_result.down)
                        
                        
                        
                        # PR-AUC
                        pr_predictions.up = -predictions.up
                        pr.up <- pr.curve(scores.class0 = pr_predictions.up[actual.up == 1],
                                          scores.class1 = pr_predictions.up[actual.up == 0],
                                          curve = TRUE)
                        pr_auc_results[["PR_AUC_up"]][n] = pr.up$auc.integral 
                        
                        pr_predictions.down = -predictions.down
                        pr.down <- pr.curve(scores.class0 = pr_predictions.down[actual.down == 1],
                                            scores.class1 = pr_predictions.down[actual.down == 0],
                                            curve = TRUE)
                        pr_auc_results[["PR_AUC_down"]][n] = pr.down$auc.integral 
                        
                        
                        
                        
                        
                } else {
                        # ROC-AUC
                        predictions = qv[, n]
                        names(predictions) = rownames(rec[[n]])
                        temp = new.table[[n]]
                        ind.temp = (rownames(temp) %in% names(predictions))
                        temp = temp[ind.temp, ]
                        
                        actual = rep(0, nrow(temp))
                        ind = (rownames(temp) %in% altered.eRNA[[n]]$altered_eRNAs)
                        actual[ind] = 1
                        
                        roc_result = roc(actual, predictions, direction = ">", quiet = TRUE)
                        auc_results[["AUC"]][n] = auc(roc_result)
                        
                        
                        # PR-AUC
                        pr_predictions = -predictions
                        pr <- pr.curve(scores.class0 = pr_predictions[actual == 1],
                                       scores.class1 = pr_predictions[actual == 0],
                                       curve = TRUE)
                        pr_auc_results[["PR_AUC"]][n] = pr$auc.integral
                        
                        
                }
                
                
                
                
        }
        
        
        stat = cbind(sensitivity = sens.single.run, specificity = speci.single.run, precision = preci.single.run)
        rownames(stat) = unique(dif.eRNA$cell_type)
        stat_final = list(stat = stat, ROC_AUC = auc_results, PR_AUC = pr_auc_results)
        
        return(list(differential_eRNA = dif.eRNA, performance = stat_final))
        
}



#############################
# get stats for TOAST model #
#############################
get.stats_TOAST <- function(rec, altered.eRNA, num, logFC_cutoff, qv_cutoff, new.table, direction)
{
        if (direction == "up") {
                for(i in seq_along(rec))
                {
                        rec_temp = rec[[i]]
                        colnames(rec_temp)[c(1,6)] = c("logFC", "P.Value")
                        rec_temp$one_side_qv = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC <= 0,]$one_side_qv = 1-rec_temp[rec_temp$logFC <= 0,]$one_side_qv
                        rec_temp$qv = p.adjust(rec_temp$one_side_qv, method = "fdr")
                        rec[[i]] = rec_temp
                }
                
        } else if (direction == "down") {
                for(i in seq_along(rec))
                {
                        rec_temp = rec[[i]]
                        colnames(rec_temp)[c(1,6)] = c("logFC", "P.Value")
                        rec_temp$one_side_qv = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC >= 0,]$one_side_qv = 1-rec_temp[rec_temp$logFC >= 0,]$one_side_qv
                        rec_temp$qv = p.adjust(rec_temp$one_side_qv, method = "fdr")
                        rec[[i]] = rec_temp
                }
                
        } else if (direction == "bi") {
                for(i in seq_along(rec))
                {
                        rec_temp = rec[[i]]
                        colnames(rec_temp)[c(1,6)] = c("logFC", "P.Value")
                        
                        # calculate qv as in upregulated direction
                        rec_temp$one_side_qv.up = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC <= 0,]$one_side_qv.up = 1-rec_temp[rec_temp$logFC <= 0,]$one_side_qv.up
                        rec_temp$qv.up = p.adjust(rec_temp$one_side_qv.up, method = "fdr")
                        
                        # calculate qv as in downregulated direction
                        rec_temp$one_side_qv.down = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC >= 0,]$one_side_qv.down = 1-rec_temp[rec_temp$logFC >= 0,]$one_side_qv.down
                        rec_temp$qv.down = p.adjust(rec_temp$one_side_qv.down, method = "fdr")
                        
                        rec[[i]] = rec_temp
                }
                
        }
        
        
        if (direction == "bi") {
                
                qv.up = cbind(rec[[1]]$qv.up, 
                              rec[[2]]$qv.up,
                              rec[[3]]$qv.up,
                              rec[[4]]$qv.up)
                rownames(qv.up) = rownames(rec[[1]])
                colnames(qv.up) = names(rec)
                
                qv.down = cbind(rec[[1]]$qv.down, 
                                rec[[2]]$qv.down,
                                rec[[3]]$qv.down,
                                rec[[4]]$qv.down)
                rownames(qv.down) = rownames(rec[[1]])
                colnames(qv.down) = names(rec)
                
                estimate = cbind(rec[[1]]$logFC, 
                                 rec[[2]]$logFC,
                                 rec[[3]]$logFC,
                                 rec[[4]]$logFC)
                rownames(estimate) = rownames(rec[[1]])
                colnames(estimate) = names(rec)
                
                dif.eRNA.up = NULL
                dif.eRNA.down = NULL
                for (j in 1: ncol(estimate))
                {
                        temp.up = cbind.data.frame(estimate = estimate[, j], q_value.up = qv.up[, j])
                        temp.up = temp.up[complete.cases(temp.up),]
                        sig.up = temp.up[temp.up$estimate > logFC_cutoff & temp.up$q_value.up < qv_cutoff, ]
                        
                        temp.down = cbind.data.frame(estimate = estimate[, j], q_value.down = qv.down[, j])
                        temp.down = temp.down[complete.cases(temp.down),]
                        sig.down = temp.down[temp.down$estimate < -logFC_cutoff & temp.down$q_value.down < qv_cutoff, ]
                        
                        temp.up = cbind.data.frame(eRNA = rownames(sig.up),
                                                   cell_type = rep(colnames(qv.up)[j], nrow(sig.up)),
                                                   direction = rep("up", nrow(sig.up)),
                                                   qv = sig.up$q_value.up,
                                                   logFC = sig.up$estimate)
                        temp.down = cbind.data.frame(eRNA = rownames(sig.down),
                                                     cell_type = rep(colnames(qv.down)[j], nrow(sig.down)),
                                                     direction = rep("down", nrow(sig.down)),
                                                     qv = sig.down$q_value.down,
                                                     logFC = sig.down$estimate)
                        
                        
                        up.eRNA.ind = (temp.up$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "up",]$altered_eRNAs)
                        down.eRNA.ind = (temp.down$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "down",]$altered_eRNAs)
                        
                        
                        
                        # Check if temp.up has rows
                        if (nrow(temp.up) > 0) {
                                temp.up$group = paste("false", temp.up$direction, sep = " ")
                                temp.up$group[up.eRNA.ind] = "true up"
                                dif.eRNA.up = rbind.data.frame(dif.eRNA.up, temp.up)
                        }
                        
                        # Check if temp.down has rows
                        if (nrow(temp.down) > 0) {
                                temp.down$group = paste("false", temp.down$direction, sep = " ")
                                temp.down$group[down.eRNA.ind] = "true down"
                                dif.eRNA.down = rbind.data.frame(dif.eRNA.down, temp.down)
                        }
                        
                        
                }
                
                # Combine the results if any rows exist
                if (exists("dif.eRNA.up") & exists("dif.eRNA.down")) {
                        dif.eRNA = rbind(dif.eRNA.up, dif.eRNA.down)
                } else if (exists("dif.eRNA.up")) {
                        dif.eRNA = dif.eRNA.up
                } else if (exists("dif.eRNA.down")) {
                        dif.eRNA = dif.eRNA.down
                } else {
                        dif.eRNA = data.frame()  # Empty dataframe if no rows exist
                }
                
                
                
                
        } else {
                qv = cbind(rec[[1]]$qv, 
                           rec[[2]]$qv,
                           rec[[3]]$qv,
                           rec[[4]]$qv)
                rownames(qv) = rownames(rec[[1]])
                colnames(qv) = names(rec)
                
                estimate = cbind(rec[[1]]$logFC, 
                                 rec[[2]]$logFC,
                                 rec[[3]]$logFC,
                                 rec[[4]]$logFC)
                rownames(estimate) = rownames(rec[[1]])
                colnames(estimate) = names(rec)
                
                
                
                dif.eRNA.up = NULL
                dif.eRNA.down = NULL
                for (j in 1: ncol(estimate))
                {
                        temp = cbind.data.frame(estimate = estimate[, j], q_value = qv[, j])
                        rownames(temp) = rownames(rec[[j]])
                        temp = temp[complete.cases(temp),]
                        sig.up = temp[temp$estimate > logFC_cutoff & temp$q_value < qv_cutoff, ]
                        sig.down = temp[temp$estimate < -logFC_cutoff & temp$q_value < qv_cutoff, ]
                        
                        temp.up = cbind.data.frame(eRNA = rownames(sig.up),
                                                   cell_type = rep(names(rec)[j], nrow(sig.up)),
                                                   direction = rep("up", nrow(sig.up)),
                                                   qv = sig.up$q_value,
                                                   logFC = sig.up$estimate)
                        temp.down = cbind.data.frame(eRNA = rownames(sig.down),
                                                     cell_type = rep(names(rec)[j], nrow(sig.down)),
                                                     direction = rep("down", nrow(sig.down)),
                                                     qv = sig.down$q_value,
                                                     logFC = sig.down$estimate)
                        
                        up.eRNA.ind = (temp.up$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "up",]$altered_eRNAs)
                        down.eRNA.ind = (temp.down$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "down",]$altered_eRNAs)
                        
                        
                        
                        # Check if temp.up has rows
                        if (nrow(temp.up) > 0) {
                                temp.up$group = paste("false", temp.up$direction, sep = " ")
                                temp.up$group[up.eRNA.ind] = "true up"
                                dif.eRNA.up = rbind.data.frame(dif.eRNA.up, temp.up)
                        }
                        
                        # Check if temp.down has rows
                        if (nrow(temp.down) > 0) {
                                temp.down$group = paste("false", temp.down$direction, sep = " ")
                                temp.down$group[down.eRNA.ind] = "true down"
                                dif.eRNA.down = rbind.data.frame(dif.eRNA.down, temp.down)
                        }
                        
                        
                }
                
                # Combine the results if any rows exist
                if (exists("dif.eRNA.up") & exists("dif.eRNA.down")) {
                        dif.eRNA = rbind(dif.eRNA.up, dif.eRNA.down)
                } else if (exists("dif.eRNA.up")) {
                        dif.eRNA = dif.eRNA.up
                } else if (exists("dif.eRNA.down")) {
                        dif.eRNA = dif.eRNA.down
                } else {
                        dif.eRNA = data.frame()  # Empty dataframe if no rows exist
                }
        }
        
        
        sens.single.run = NULL
        speci.single.run = NULL
        preci.single.run = NULL
        auc_results = list()
        pr_auc_results = list()
        
        for (n in seq_along(unique(dif.eRNA$cell_type)))
        {
                sub = dif.eRNA[dif.eRNA$cell_type == unique(dif.eRNA$cell_type)[n], ]
                sens = nrow(sub[sub$group == "true up" | sub$group == "true down", ]) / num
                sens.single.run = c(sens.single.run, sens)
                
                speci = (nrow(new.table[[n]]) - num - nrow(sub[sub$group == "false up" | sub$group == "false down", ])) / (nrow(new.table[[n]]) - num)
                speci.single.run = c(speci.single.run, speci)
                
                preci = nrow(sub[sub$group == "true up" | sub$group == "true down", ]) / (nrow(sub[sub$group == "true up" | sub$group == "true down", ]) + nrow(sub[sub$group == "false up" | sub$group == "false down", ]))
                preci.single.run = c(preci.single.run, preci)
                
                
                
                
                # AUC
                
                if (direction == "bi") {
                        
                        # ROC_AUC
                        predictions.up = qv.up[, n]
                        names(predictions.up) = rownames(rec[[n]])
                        
                        predictions.down = qv.down[, n]
                        names(predictions.down) = rownames(rec[[n]])
                        
                        temp = new.table[[n]]
                        
                        ind.temp.up = (rownames(temp) %in% names(predictions.up))
                        temp.up = temp[ind.temp.up, ]
                        
                        ind.temp.down = (rownames(temp) %in% names(predictions.down))
                        temp.down = temp[ind.temp.down, ]
                        
                        actual.up = rep(0, nrow(temp.up))
                        ind.up = (rownames(temp.up) %in% altered.eRNA[[n]][altered.eRNA[[n]]$direction == "up",]$altered_eRNAs)
                        actual.up[ind.up] = 1
                        
                        roc_result.up = roc(actual.up, predictions.up, direction = ">", quiet = TRUE)
                        auc_results[["AUC_up"]][n] = auc(roc_result.up)
                        
                        actual.down = rep(0, nrow(temp.down))
                        ind.down = (rownames(temp.down) %in% altered.eRNA[[n]][altered.eRNA[[n]]$direction == "down",]$altered_eRNAs)
                        actual.down[ind.down] = 1
                        
                        roc_result.down = roc(actual.down, predictions.down, direction = ">", quiet = TRUE)
                        auc_results[["AUC_down"]][n] = auc(roc_result.down)
                        
                        
                        
                        # PR-AUC
                        pr_predictions.up = -predictions.up
                        pr.up <- pr.curve(scores.class0 = pr_predictions.up[actual.up == 1],
                                          scores.class1 = pr_predictions.up[actual.up == 0],
                                          curve = TRUE)
                        pr_auc_results[["PR_AUC_up"]][n] = pr.up$auc.integral 
                        
                        pr_predictions.down = -predictions.down
                        pr.down <- pr.curve(scores.class0 = pr_predictions.down[actual.down == 1],
                                            scores.class1 = pr_predictions.down[actual.down == 0],
                                            curve = TRUE)
                        pr_auc_results[["PR_AUC_down"]][n] = pr.down$auc.integral 
                        
                        
                        
                        
                        
                } else {
                        # ROC-AUC
                        predictions = qv[, n]
                        names(predictions) = rownames(rec[[n]])
                        temp = new.table[[n]]
                        ind.temp = (rownames(temp) %in% names(predictions))
                        temp = temp[ind.temp, ]
                        
                        actual = rep(0, nrow(temp))
                        ind = (rownames(temp) %in% altered.eRNA[[n]]$altered_eRNAs)
                        actual[ind] = 1
                        
                        roc_result = roc(actual, predictions, direction = ">", quiet = TRUE)
                        auc_results[["AUC"]][n] = auc(roc_result)
                        
                        
                        # PR-AUC
                        pr_predictions = -predictions
                        pr <- pr.curve(scores.class0 = pr_predictions[actual == 1],
                                       scores.class1 = pr_predictions[actual == 0],
                                       curve = TRUE)
                        pr_auc_results[["PR_AUC"]][n] = pr$auc.integral
                        
                        
                }
                
                
                
                
        }
        
        
        stat = cbind(sensitivity = sens.single.run, specificity = speci.single.run, precision = preci.single.run)
        rownames(stat) = unique(dif.eRNA$cell_type)
        stat_final = list(stat = stat, ROC_AUC = auc_results, PR_AUC = pr_auc_results)
        
        return(list(differential_eRNA = dif.eRNA, performance = stat_final))
        
}





######################################
# get stats for our CDEnhancer model #
######################################
get.stats <- function(rec, altered.eRNA, num, logFC_cutoff, qv_cutoff, new.table, direction, cv)
{
        if (direction == "up") {
                for(i in 1:(length(rec)/2))
                {
                        rec_temp = cbind(rec[[i]], rec[[i + (length(rec)/2)]])
                        colnames(rec_temp)[7:12] = paste0(colnames(rec_temp)[7:12], "_int")
                        
                        rec_temp$one_side_pv = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC <= 0,]$one_side_pv = 1-rec_temp[rec_temp$logFC <= 0,]$one_side_pv
                        rec_temp$qv = p.adjust(rec_temp$one_side_pv, method = "fdr")
                        rec_temp$dir = rec_temp$logFC * rec_temp$logFC_int
                        rec_temp$cv = cv[,i]
                        rec[[i]] = rec_temp
                        
                }
                rec = rec[1:4]
                
        } else if (direction == "down") {
                for(i in 1:(length(rec)/2))
                {
                        rec_temp = cbind(rec[[i]], rec[[i + (length(rec)/2)]])
                        colnames(rec_temp)[7:12] = paste0(colnames(rec_temp)[7:12], "_int")
                        
                        rec_temp$one_side_pv = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC >= 0,]$one_side_pv = 1-rec_temp[rec_temp$logFC >= 0,]$one_side_pv
                        rec_temp$qv = p.adjust(rec_temp$one_side_pv, method = "fdr")
                        rec_temp$dir = rec_temp$logFC * rec_temp$logFC_int
                        rec_temp$cv = cv[,i]
                        rec[[i]] = rec_temp
                }
                rec = rec[1:4]
                
        } else if (direction == "bi") {
                for(i in 1:(length(rec)/2))
                {
                        rec_temp = cbind(rec[[i]], rec[[i + (length(rec)/2)]])
                        colnames(rec_temp)[7:12] = paste0(colnames(rec_temp)[7:12], "_int")
                        
                        # calculate qv as in upregulated direction
                        rec_temp$one_side_pv.up = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC <= 0,]$one_side_pv.up = 1-rec_temp[rec_temp$logFC <= 0,]$one_side_pv.up
                        rec_temp$qv.up = p.adjust(rec_temp$one_side_pv.up, method = "fdr")
                        
                        # calculate qv as in downregulated direction
                        rec_temp$one_side_pv.down = rec_temp$P.Value/2
                        rec_temp[rec_temp$logFC >= 0,]$one_side_pv.down = 1-rec_temp[rec_temp$logFC >= 0,]$one_side_pv.down
                        rec_temp$qv.down = p.adjust(rec_temp$one_side_pv.down, method = "fdr")
                        rec_temp$dir = rec_temp$logFC * rec_temp$logFC_int
                        rec_temp$cv = cv[,i]
                        rec[[i]] = rec_temp
                }
                rec = rec[1:4]
                
        }
        
        
        if (direction == "bi") {
                
                qv.up = cbind(rec[[1]]$qv.up, 
                              rec[[2]]$qv.up,
                              rec[[3]]$qv.up,
                              rec[[4]]$qv.up)
                rownames(qv.up) = rownames(rec[[1]])
                colnames(qv.up) = names(rec)[1:4]
                
                qv.down = cbind(rec[[1]]$qv.down, 
                                rec[[2]]$qv.down,
                                rec[[3]]$qv.down,
                                rec[[4]]$qv.down)
                rownames(qv.down) = rownames(rec[[1]])
                colnames(qv.down) = names(rec)[1:4]
                
                estimate = cbind(rec[[1]]$logFC, 
                                 rec[[2]]$logFC,
                                 rec[[3]]$logFC,
                                 rec[[4]]$logFC)
                rownames(estimate) = rownames(rec[[1]])
                colnames(estimate) = names(rec)[1:4]
                
                dir = cbind(rec[[1]]$dir,
                            rec[[2]]$dir,
                            rec[[3]]$dir,
                            rec[[4]]$dir)
                rownames(dir) = rownames(rec[[1]])
                colnames(dir) = names(rec)[1:4]
                
                one_side_pv.up = cbind(rec[[1]]$one_side_pv.up,
                                       rec[[2]]$one_side_pv.up,
                                       rec[[3]]$one_side_pv.up,
                                       rec[[4]]$one_side_pv.up)
                rownames(one_side_pv.up) = rownames(rec[[1]])
                colnames(one_side_pv.up) = names(rec)[1:4]
                
                one_side_pv.down = cbind(rec[[1]]$one_side_pv.down,
                                         rec[[2]]$one_side_pv.down,
                                         rec[[3]]$one_side_pv.down,
                                         rec[[4]]$one_side_pv.down)
                rownames(one_side_pv.down) = rownames(rec[[1]])
                colnames(one_side_pv.down) = names(rec)[1:4]
                
                cv = cbind(rec[[1]]$cv,
                           rec[[2]]$cv,
                           rec[[3]]$cv,
                           rec[[4]]$cv)
                rownames(cv) = rownames(rec[[1]])
                colnames(cv) = names(rec)[1:4]
                
                
                
                dif.eRNA.up = NULL
                dif.eRNA.down = NULL
                for (j in 1: ncol(estimate))
                {
                        temp.up = cbind.data.frame(estimate = estimate[,j], one_side_pv.up =one_side_pv.up[,j], dir = dir[,j], cv = cv[,j])
                        rownames(temp.up) = rownames(rec[[j]])
                        temp.up$q_value.up = p.adjust(temp.up$one_side_pv.up, method = "fdr")
                        temp.up = temp.up[complete.cases(temp.up),]
                        
                        sig.up = temp.up[temp.up$estimate > logFC_cutoff & temp.up$q_value.up < qv_cutoff, ]
                        
                        temp.down = cbind.data.frame(estimate = estimate[,j], one_side_pv.down =one_side_pv.down[,j], dir = dir[,j], cv = cv[,j])
                        rownames(temp.down) = rownames(rec[[j]])
                        temp.down$q_value.down = p.adjust(temp.down$one_side_pv.down, method = "fdr")
                        temp.down = temp.down[complete.cases(temp.down),]
                        
                        sig.down = temp.down[temp.down$estimate < -logFC_cutoff & temp.down$q_value.down < qv_cutoff, ]
                        
                        temp.up = cbind.data.frame(eRNA = rownames(sig.up),
                                                   cell_type = rep(colnames(qv.up)[j], nrow(sig.up)),
                                                   direction = rep("up", nrow(sig.up)),
                                                   qv = sig.up$q_value.up,
                                                   logFC = sig.up$estimate)
                        temp.down = cbind.data.frame(eRNA = rownames(sig.down),
                                                     cell_type = rep(colnames(qv.down)[j], nrow(sig.down)),
                                                     direction = rep("down", nrow(sig.down)),
                                                     qv = sig.down$q_value.down,
                                                     logFC = sig.down$estimate)
                        
                          
                        up.eRNA.ind = (temp.up$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "up",]$altered_eRNAs)
                        down.eRNA.ind = (temp.down$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "down",]$altered_eRNAs)
                        
                        
                        
                        # Check if temp.up has rows
                        if (nrow(temp.up) > 0) {
                                temp.up$group = paste("false", temp.up$direction, sep = " ")
                                temp.up$group[up.eRNA.ind] = "true up"
                                dif.eRNA.up = rbind.data.frame(dif.eRNA.up, temp.up)
                        }
                        
                        # Check if temp.down has rows
                        if (nrow(temp.down) > 0) {
                                temp.down$group = paste("false", temp.down$direction, sep = " ")
                                temp.down$group[down.eRNA.ind] = "true down"
                                dif.eRNA.down = rbind.data.frame(dif.eRNA.down, temp.down)
                        }
                        
                        
                }
                
                # Combine the results if any rows exist
                if (exists("dif.eRNA.up") & exists("dif.eRNA.down")) {
                        dif.eRNA = rbind(dif.eRNA.up, dif.eRNA.down)
                } else if (exists("dif.eRNA.up")) {
                        dif.eRNA = dif.eRNA.up
                } else if (exists("dif.eRNA.down")) {
                        dif.eRNA = dif.eRNA.down
                } else {
                        dif.eRNA = data.frame()  # Empty dataframe if no rows exist
                }
                
                
                
                
        } else {
                
                qv = cbind(rec[[1]]$qv, 
                           rec[[2]]$qv,
                           rec[[3]]$qv,
                           rec[[4]]$qv)
                rownames(qv) = rownames(rec[[1]])
                colnames(qv) = names(rec)[1:4]
                
                estimate = cbind(rec[[1]]$logFC, 
                                 rec[[2]]$logFC,
                                 rec[[3]]$logFC,
                                 rec[[4]]$logFC)
                rownames(estimate) = rownames(rec[[1]])
                colnames(estimate) = names(rec)[1:4]
                
                dir = cbind(rec[[1]]$dir,
                            rec[[2]]$dir,
                            rec[[3]]$dir,
                            rec[[4]]$dir)
                rownames(dir) = rownames(rec[[1]])
                colnames(dir) = names(rec)[1:4]
                
                one_side_pv = cbind(rec[[1]]$one_side_pv,
                                    rec[[2]]$one_side_pv,
                                    rec[[3]]$one_side_pv,
                                    rec[[4]]$one_side_pv)
                rownames(one_side_pv) = rownames(rec[[1]])
                colnames(one_side_pv) = names(dir)[1:4]
                
                cv = cbind(rec[[1]]$cv,
                           rec[[2]]$cv,
                           rec[[3]]$cv,
                           rec[[4]]$cv)
                rownames(cv) = rownames(rec[[1]])
                colnames(cv) = names(rec)[1:4]
                
                
                
                
                
                dif.eRNA.up = NULL
                dif.eRNA.down = NULL
                for (j in 1: length(rec))
                {
                        temp = cbind.data.frame(estimate = estimate[,j], one_side_pv =one_side_pv[,j], dir = dir[,j], cv = cv[,j])
                        rownames(temp) = rownames(rec[[j]])
                        temp$q_value = p.adjust(temp$one_side_pv, method = "fdr")
                        temp = temp[complete.cases(temp),]
                        sig.up = temp[temp$estimate > logFC_cutoff & temp$q_value < qv_cutoff, ]
                        sig.down = temp[temp$estimate < -logFC_cutoff & temp$q_value < qv_cutoff, ]
                        
                        temp.up = cbind.data.frame(eRNA = rownames(sig.up),
                                                   cell_type = rep(names(rec)[j], nrow(sig.up)),
                                                   direction = rep("up", nrow(sig.up)),
                                                   qv = sig.up$q_value,
                                                   logFC = sig.up$estimate)
                        temp.down = cbind.data.frame(eRNA = rownames(sig.down),
                                                     cell_type = rep(names(rec)[j], nrow(sig.down)),
                                                     direction = rep("down", nrow(sig.down)),
                                                     qv = sig.down$q_value,
                                                     logFC = sig.down$estimate)
                        
                        up.eRNA.ind = (temp.up$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "up",]$altered_eRNAs)
                        down.eRNA.ind = (temp.down$eRNA %in% altered.eRNA[[j]][altered.eRNA[[1]]$direction == "down",]$altered_eRNAs)
                        
                        
                        
                        # Check if temp.up has rows
                        if (nrow(temp.up) > 0) {
                                temp.up$group = paste("false", temp.up$direction, sep = " ")
                                temp.up$group[up.eRNA.ind] = "true up"
                                dif.eRNA.up = rbind.data.frame(dif.eRNA.up, temp.up)
                        }
                        
                        # Check if temp.down has rows
                        if (nrow(temp.down) > 0) {
                                temp.down$group = paste("false", temp.down$direction, sep = " ")
                                temp.down$group[down.eRNA.ind] = "true down"
                                dif.eRNA.down = rbind.data.frame(dif.eRNA.down, temp.down)
                        }
                        
                        
                }
                
                # Combine the results if any rows exist
                if (exists("dif.eRNA.up") & exists("dif.eRNA.down")) {
                        dif.eRNA = rbind(dif.eRNA.up, dif.eRNA.down)
                } else if (exists("dif.eRNA.up")) {
                        dif.eRNA = dif.eRNA.up
                } else if (exists("dif.eRNA.down")) {
                        dif.eRNA = dif.eRNA.down
                } else {
                        dif.eRNA = data.frame()  # Empty dataframe if no rows exist
                }
        }
        
        
        sens.single.run = NULL
        speci.single.run = NULL
        preci.single.run = NULL
        auc_results = list()
        pr_auc_results = list()
        
        for (n in seq_along(unique(dif.eRNA$cell_type)))
        {
                sub = dif.eRNA[dif.eRNA$cell_type == unique(dif.eRNA$cell_type)[n], ]
                sens = nrow(sub[sub$group == "true up" | sub$group == "true down", ]) / num
                sens.single.run = c(sens.single.run, sens)
                
                speci = (nrow(new.table[[n]]) - num - nrow(sub[sub$group == "false up" | sub$group == "false down", ])) / (nrow(new.table[[n]]) - num)
                speci.single.run = c(speci.single.run, speci)
                
                preci = nrow(sub[sub$group == "true up" | sub$group == "true down", ]) / (nrow(sub[sub$group == "true up" | sub$group == "true down", ]) + nrow(sub[sub$group == "false up" | sub$group == "false down", ]))
                preci.single.run = c(preci.single.run, preci)
                
                
                
                
                # AUC
                
                if (direction == "bi") {
                        
                        # ROC_AUC
                        predictions.up = qv.up[, n]
                        names(predictions.up) = rownames(rec[[n]])
                        
                        predictions.down = qv.down[, n]
                        names(predictions.down) = rownames(rec[[n]])
                        
                        temp = new.table[[n]]
                        
                        ind.temp.up = (rownames(temp) %in% names(predictions.up))
                        temp.up = temp[ind.temp.up, ]
                        
                        ind.temp.down = (rownames(temp) %in% names(predictions.down))
                        temp.down = temp[ind.temp.down, ]
                        
                        actual.up = rep(0, nrow(temp.up))
                        ind.up = (rownames(temp.up) %in% altered.eRNA[[n]][altered.eRNA[[n]]$direction == "up",]$altered_eRNAs)
                        actual.up[ind.up] = 1
                        
                        roc_result.up = roc(actual.up, predictions.up, direction = ">", quiet = TRUE)
                        auc_results[["AUC_up"]][n] = auc(roc_result.up)
                        
                        actual.down = rep(0, nrow(temp.down))
                        ind.down = (rownames(temp.down) %in% altered.eRNA[[n]][altered.eRNA[[n]]$direction == "down",]$altered_eRNAs)
                        actual.down[ind.down] = 1
                        
                        roc_result.down = roc(actual.down, predictions.down, direction = ">", quiet = TRUE)
                        auc_results[["AUC_down"]][n] = auc(roc_result.down)
                        
                        
                        
                        # PR-AUC
                        pr_predictions.up = -predictions.up
                        pr.up <- pr.curve(scores.class0 = pr_predictions.up[actual.up == 1],
                                          scores.class1 = pr_predictions.up[actual.up == 0],
                                          curve = TRUE)
                        pr_auc_results[["PR_AUC_up"]][n] = pr.up$auc.integral 
                        
                        pr_predictions.down = -predictions.down
                        pr.down <- pr.curve(scores.class0 = pr_predictions.down[actual.down == 1],
                                            scores.class1 = pr_predictions.down[actual.down == 0],
                                            curve = TRUE)
                        pr_auc_results[["PR_AUC_down"]][n] = pr.down$auc.integral 
                        
                        
                        
                        
                        
                } else {
                        # ROC-AUC
                        predictions = qv[, n]
                        names(predictions) = rownames(rec[[n]])
                        temp = new.table[[n]]
                        ind.temp = (rownames(temp) %in% names(predictions))
                        temp = temp[ind.temp, ]
                        
                        actual = rep(0, nrow(temp))
                        ind = (rownames(temp) %in% altered.eRNA[[n]]$altered_eRNAs)
                        actual[ind] = 1
                        
                        roc_result = roc(actual, predictions, direction = ">", quiet = TRUE)
                        auc_results[["AUC"]][n] = auc(roc_result)
                        
                        
                        # PR-AUC
                        pr_predictions = -predictions
                        pr <- pr.curve(scores.class0 = pr_predictions[actual == 1],
                                       scores.class1 = pr_predictions[actual == 0],
                                       curve = TRUE)
                        pr_auc_results[["PR_AUC"]][n] = pr$auc.integral
                        
                        
                }
                
                
                
                
        }
        
        
        stat = cbind(sensitivity = sens.single.run, specificity = speci.single.run, precision = preci.single.run)
        rownames(stat) = unique(dif.eRNA$cell_type)
        stat_final = list(stat = stat, ROC_AUC = auc_results, PR_AUC = pr_auc_results)
        
        return(list(differential_eRNA = dif.eRNA, performance = stat_final))
        
}






eRNA_model_test <- function(data.list,
                            num,
                            snr,
                            direction,
                            cf,
                            var_ctrl = 1,
                            logFC_cutoff,
                            qv_cutoff,
                            ref_idx = 4 # <--- which cell type is reference (1–4)
)
{
        # randomly split data into control and test groups
        sample_size = floor(0.5*ncol(data.list[[1]]))
        
        #set.seed(123)
        inc.ct = sample(seq_len(ncol(data.list[[1]])),size = sample_size)
        ct = colnames(data.list[[1]])[inc.ct]
        test = colnames(data.list[[1]])[-inc.ct]
        label = cbind(samples = c(ct, test),
                      type = c(rep("CT", length(ct)), rep("test", length(test))))
        
        # select N number of eRNAs to alter basing on SNR
        altered = altered_eRNA_table(data.list, ct, test, num, snr, direction) 
        new.table = altered$table.list
        altered.eRNA = altered$altered.eRNA
        
        
        # creat the bulk eRNA file
        bulk.dt = sweep(altered[[1]][[1]], 2, cf$true$V1, `*`) +
                sweep(altered[[1]][[2]], 2, cf$true$V2, `*`) + 
                sweep(altered[[1]][[3]], 2, cf$true$V3, `*`) + 
                sweep(altered[[1]][[4]], 2, cf$true$V4, `*`)
        
        bulk.dt = bulk.dt[apply(bulk.dt, 1, mean) > 1, ]
        
        
        
        
        # take log of the bulk data
        log.bulk.dt = log2(bulk.dt + 0.01)
        log.bulk.dt = log.bulk.dt[complete.cases(log.bulk.dt), ]
        
        eRNA = log.bulk.dt
        frac = cf$pred
        
        ## --- NEW BLOCK: reorder so that ref_idx is last (4th) ---
        K <- ncol(frac)  # should be 4
        perm <- c(setdiff(seq_len(K), ref_idx), ref_idx)
        frac <- frac[, perm, drop = FALSE]
        altered.eRNA <- altered.eRNA[perm]
        new.table <- new.table[perm]
        ## --------------------------------------------------------
        
        
        
        colnames(frac) = names(new.table)
        rownames(frac) = colnames(bulk.dt)
        
        all = cbind(c(ct, test), c(rep("CT", length(ct)), rep("Test", length(test))))
        ind.all = match(colnames(bulk.dt), all[,1])
        all = all[ind.all,]
        type = as.numeric(all[,2] == "Test") # 1 is test, 0 is ct
        
        
        ## build data.frame for model.matrix
        df <- data.frame(
                type = type,
                f1   = frac[, 1],
                f2   = frac[, 2],
                f3   = frac[, 3]
        )
        
        ## design matrix for the full model:
        ## y = δ + β¹*C + γ¹*f1 + γ²*f2 + γ³*f3 + θ¹*(f1*C) + θ²*(f2*C) + θ³*(f3*C)
        design <- model.matrix(~ type * (f1 + f2 + f3), data = df)
        colnames(design)[c(6,7,8)] = c("interaction1", "interaction2", "interaction3")
        colnames(design)[colnames(design) == "(Intercept)"] <- "Intercept"
        
        ## fit limma
        fit <- lmFit(eRNA, design)
        
        ## contrasts:
        ## β¹ alone
        ## β¹ + θ¹, β¹ + θ², β¹ + θ³
        contrasts.mat <- makeContrasts(
                eff_f1 = type + interaction1,    # β¹ + θ¹ summation coeficient for cell type 1
                eff_f2 = type + interaction2,    # β¹ + θ² summation coeficient for cell type 2
                eff_f3 = type + interaction3,    # β¹ + θ³ summation coeficient for cell type 3
                eff_f4 = type,                   # β¹  coeficient for cell type 4
                int1   = interaction1,
                int2   = interaction2,
                int3   = interaction3,
                int4   = type,
                levels = design
        )
        
        fit2 <- contrasts.fit(fit, contrasts.mat)
        fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
        
        
        SE_all <- fit2$stdev.unscaled * sqrt(fit2$s2.post)
        SE_sum = SE_all[, 1:4]
        coef_sum = fit2$coefficients[,1:4]
        cv = SE_sum / abs(coef_sum)
        
        
        
        ## output: rec$beta_main, rec$eff_f1, rec$eff_f2, rec$eff_f3
        rec <- list()
        
        for (coef_name in colnames(contrasts.mat)) {
                tmp <- topTable(
                        fit2, coef = coef_name,
                        adjust.method = "fdr",
                        number = Inf
                )
                tmp <- tmp[rownames(eRNA), ]  # reorder to match eRNA
                rec[[coef_name]] <- tmp
        }
        
        ## rec now contains:
        ## rec$beta_main  -> β¹ (main condition effect)
        ## rec$eff_f1     -> β¹ + θ¹
        ## rec$eff_f2     -> β¹ + θ²
        ## rec$eff_f3     -> β¹ + θ³
        
        names(rec) = c(names(new.table), paste0("int_", names(new.table)))
        performance = get.stats(rec, altered.eRNA, num, logFC_cutoff, qv_cutoff, new.table, direction, cv)
        
        
        #Test in CellDMC
        cdmc.o <- CellDMC(eRNA, pheno.v=type+1, frac.m=as.matrix(frac), 
                          adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, sort = FALSE, mc.cores = 4);
        rec.CellDMC = cdmc.o$coe
        performance.CellDMC = get.stats_CellDMC(rec.CellDMC, altered.eRNA, num, logFC_cutoff, qv_cutoff, new.table, direction)
        
        
        #Test in TOAST
        design.TOAST <- data.frame(type = as.factor(type))
        Design_out <- makeDesign(design.TOAST, frac)
        fitted_model <- fitModel(Design_out, eRNA)
        res_cstest <- csTest(fitted_model, coef = "type", sort = FALSE)
        rec.TOAST = res_cstest[1:4]
        performance.TOAST = get.stats_TOAST(rec.TOAST, altered.eRNA, num, logFC_cutoff, qv_cutoff, new.table, direction)
        
        
        
        
        
        
        
        return(list(
                performance = performance,
                altered_eRNAs_in_each_tissue = altered.eRNA,
                performance.CellDMC = performance.CellDMC,
                performance.TOAST = performance.TOAST
        ))
        
}





run_param_grid <- function(
                data.list,
                var_ctrl,
                r,
                logFC_cutoff = 1,
                qv_cutoff    = 0.05,
                ref_idx      = 4,
                out_dir,
                iterations   = 100
                
) {
        # generate cf
        cf <- cf_generation(
                var_ctrl = var_ctrl,
                r = r,
                out_dir = out_dir
        )
        
        
        # Build parameter grid
        nums       <- c(1000, 2000, 3000)
        snrs       <- c(0.75, 1, 2, 3)
        directions <- c("up", "down", "bi")
        
        param_grid <- expand.grid(
                num       = nums,
                snr       = snrs,
                direction = directions,
                stringsAsFactors = FALSE
        )
        
        # Loop over combinations
        for (i in seq_len(nrow(param_grid))) {
                
                current_num <- param_grid$num[i]
                current_snr <- param_grid$snr[i]
                current_dir <- param_grid$direction[i]
                
                cat("\nRunning combination:",
                    "num =", current_num,
                    "snr =", current_snr,
                    "direction =", current_dir, "\n")
                
                # Storage for iterations
                combination_results <- vector("list", iterations)
                
                for (j in seq_len(iterations)) {
                        time_taken <- system.time({
                                combination_results[[j]] <- eRNA_model_test(
                                        data.list,
                                        current_num,
                                        current_snr,
                                        current_dir,
                                        cf,
                                        var_ctrl,
                                        logFC_cutoff,
                                        qv_cutoff,
                                        ref_idx = ref_idx
                                )
                        })
                        
                        cat("  Iteration:", j, 
                            "Time:", time_taken["elapsed"], "sec\n")
                }
                
                # output filename
                combination_name <- paste0(
                        "num_", current_num,
                        "_snr_", current_snr,
                        "_direction_", current_dir,
                        "_model_comparing"
                )
                
                file_name <- paste0(out_dir, combination_name, ".RDS")
                
                saveRDS(combination_results, file = file_name)
                
                cat("Saved to:", file_name, "\n")
        }
        
        invisible(TRUE)
}




extract_stats <- function(out_dir, var_ctrl, r, iterations)
{
        file_paths <- list.files(out_dir, pattern = "_model_comparing\\.RDS$", full.names = TRUE)
        
        stats_table_bi_all = NULL
        stats_table_up_down_all = NULL
        stats_table_bi_all_CellDMC = NULL
        stats_table_up_down_all_CellDMC = NULL
        stats_table_bi_all_TOAST = NULL
        stats_table_up_down_all_TOAST = NULL
        
        for (i in seq_along(file_paths))
        {
                temp = readRDS(file_paths[i])
                spl_name = unlist(strsplit(file_paths[i], "num_|_snr_|_direction_|_model_comparing.RDS"))
                
                stats_table = NULL
                stats_table_CellDMC = NULL
                stats_table_TOAST = NULL
                stats_table2 = NULL
                stats_table2_CellDMC = NULL
                stats_table2_TOAST = NULL
                if (spl_name[4] == "bi") {
                        for (n in seq_along(temp))
                        {
                                stats_single = temp[[n]]$performance$performance
                                stats_single_CellDMC = temp[[n]]$performance.CellDMC$performance
                                stats_single_TOAST = temp[[n]]$performance.TOAST$performance
                                
                                if (!is.null(stats_single$stat) && nrow(stats_single$stat) > 0) {
                                        stats_table_single <- cbind.data.frame(
                                                sensitivity   = stats_single$stat[, 1],
                                                specificity   = stats_single$stat[, 2],
                                                precision     = stats_single$stat[, 3],
                                                AUC_up        = stats_single$ROC_AUC$AUC_up,
                                                AUC_down      = stats_single$ROC_AUC$AUC_down,
                                                PR_AUC_up     = stats_single$PR_AUC$PR_AUC_up,
                                                PR_AUC_down   = stats_single$PR_AUC$PR_AUC_down,
                                                cell_type     = rownames(stats_single$stat),
                                                snr           = spl_name[3],
                                                num           = spl_name[2],
                                                direction     = spl_name[4]
                                        )
                                        
                                        stats_table <- rbind.data.frame(stats_table, stats_table_single)
                                }
                                
                                if (!is.null(stats_single_CellDMC$stat) && nrow(stats_single_CellDMC$stat) > 0) {
                                        stats_table_single_CellDMC <- cbind.data.frame(
                                                sensitivity   = stats_single_CellDMC$stat[, 1],
                                                specificity   = stats_single_CellDMC$stat[, 2],
                                                precision     = stats_single_CellDMC$stat[, 3],
                                                AUC_up        = stats_single_CellDMC$ROC_AUC$AUC_up,
                                                AUC_down      = stats_single_CellDMC$ROC_AUC$AUC_down,
                                                PR_AUC_up     = stats_single_CellDMC$PR_AUC$PR_AUC_up,
                                                PR_AUC_down   = stats_single_CellDMC$PR_AUC$PR_AUC_down,
                                                cell_type     = rownames(stats_single_CellDMC$stat),
                                                snr           = spl_name[3],
                                                num           = spl_name[2],
                                                direction     = spl_name[4]
                                        )
                                        
                                        stats_table_CellDMC <- rbind.data.frame(stats_table_CellDMC, stats_table_single_CellDMC)
                                }
                                
                                if (!is.null(stats_single_TOAST$stat) && nrow(stats_single_TOAST$stat) > 0) {
                                        stats_table_single_TOAST <- cbind.data.frame(
                                                sensitivity   = stats_single_TOAST$stat[, 1],
                                                specificity   = stats_single_TOAST$stat[, 2],
                                                precision     = stats_single_TOAST$stat[, 3],
                                                AUC_up        = stats_single_TOAST$ROC_AUC$AUC_up,
                                                AUC_down      = stats_single_TOAST$ROC_AUC$AUC_down,
                                                PR_AUC_up     = stats_single_TOAST$PR_AUC$PR_AUC_up,
                                                PR_AUC_down   = stats_single_TOAST$PR_AUC$PR_AUC_down,
                                                cell_type     = rownames(stats_single_TOAST$stat),
                                                snr           = spl_name[3],
                                                num           = spl_name[2],
                                                direction     = spl_name[4]
                                        )
                                        
                                        stats_table_TOAST <- rbind.data.frame(stats_table_TOAST, stats_table_single_TOAST)
                                }
                                
                        }
                        
                } else {
                        
                        for (n in seq_along(temp))
                        {
                                stats_single = temp[[n]]$performance$performance
                                stats_single_CellDMC = temp[[n]]$performance.CellDMC$performance
                                stats_single_TOAST = temp[[n]]$performance.TOAST$performance
                                
                                
                                if (!is.null(stats_single$stat) && nrow(stats_single$stat) > 0) {
                                        stats_table_single = cbind.data.frame(sensitivity = stats_single$stat[,1],
                                                                              specificity = stats_single$stat[,2],
                                                                              precision = stats_single$stat[,3],
                                                                              AUC = stats_single$ROC_AUC$AUC,
                                                                              PR_AUC = stats_single$PR_AUC$PR_AUC,
                                                                              cell_type = rownames(stats_single$stat),
                                                                              snr = spl_name[3],
                                                                              num = spl_name[2],
                                                                              direction = spl_name[4])
                                        
                                        
                                        stats_table2 = rbind.data.frame(stats_table2, stats_table_single)
                                }
                                
                                if (!is.null(stats_single_CellDMC$stat) && nrow(stats_single_CellDMC$stat) > 0) {
                                        stats_table_single_CellDMC = cbind.data.frame(sensitivity = stats_single_CellDMC$stat[,1],
                                                                                      specificity = stats_single_CellDMC$stat[,2],
                                                                                      precision = stats_single_CellDMC$stat[,3],
                                                                                      AUC = stats_single_CellDMC$ROC_AUC$AUC,
                                                                                      PR_AUC = stats_single_CellDMC$PR_AUC$PR_AUC,
                                                                                      cell_type = rownames(stats_single_CellDMC$stat),
                                                                                      snr = spl_name[3],
                                                                                      num = spl_name[2],
                                                                                      direction = spl_name[4])
                                        
                                        
                                        stats_table2_CellDMC = rbind.data.frame(stats_table2_CellDMC, stats_table_single_CellDMC)
                                }
                                
                                if (!is.null(stats_single_TOAST$stat) && nrow(stats_single_TOAST$stat) > 0) {
                                        stats_table_single_TOAST = cbind.data.frame(sensitivity = stats_single_TOAST$stat[,1],
                                                                                    specificity = stats_single_TOAST$stat[,2],
                                                                                    precision = stats_single_TOAST$stat[,3],
                                                                                    AUC = stats_single_TOAST$ROC_AUC$AUC,
                                                                                    PR_AUC = stats_single_TOAST$PR_AUC$PR_AUC,
                                                                                    cell_type = rownames(stats_single_TOAST$stat),
                                                                                    snr = spl_name[3],
                                                                                    num = spl_name[2],
                                                                                    direction = spl_name[4])
                                        
                                        
                                        stats_table2_TOAST = rbind.data.frame(stats_table2_TOAST, stats_table_single_TOAST)
                                }
                                
                                
                               
                        }
                }
                
                stats_table_bi_all = rbind.data.frame(stats_table_bi_all, stats_table)
                stats_table_up_down_all = rbind.data.frame(stats_table_up_down_all, stats_table2)
                
                stats_table_bi_all_CellDMC = rbind.data.frame(stats_table_bi_all_CellDMC, stats_table_CellDMC)
                stats_table_up_down_all_CellDMC = rbind.data.frame(stats_table_up_down_all_CellDMC, stats_table2_CellDMC)
                
                stats_table_bi_all_TOAST = rbind.data.frame(stats_table_bi_all_TOAST, stats_table_TOAST)
                stats_table_up_down_all_TOAST = rbind.data.frame(stats_table_up_down_all_TOAST, stats_table2_TOAST)
                
                
        }
        
        ls = list()
        ls[["stats_bi_direction"]] = stats_table_bi_all
        ls[["stats_up_down_directions"]] = stats_table_up_down_all
        
        ls[["stats_bi_direction_CellDMC"]] = stats_table_bi_all_CellDMC
        ls[["stats_up_down_directions_CellDMC"]] = stats_table_up_down_all_CellDMC
        
        ls[["stats_bi_direction_TOAST"]] = stats_table_bi_all_TOAST
        ls[["stats_up_down_directions_TOAST"]] = stats_table_up_down_all_TOAST
        
        
        file_name <- paste0(out_dir, "extracted_stats_var_ctrl_", var_ctrl, "target_r", r, "it", iterations, ".RDS")
        saveRDS(ls, file = file_name)
        cat("Saved to:", file_name, "\n")
        
}




run_and_extract <- function(
                data.list,
                var_ctrl,
                r,
                out_dir,
                logFC_cutoff = 1,
                qv_cutoff    = 0.05,
                ref_idx      = 4,
                iterations   = 100
) {
        # 1) run the grid and save all .RDS files
        run_param_grid(
                data.list    = data.list,
                var_ctrl     = var_ctrl,
                r            = r,
                logFC_cutoff = logFC_cutoff,
                qv_cutoff    = qv_cutoff,
                ref_idx      = ref_idx,
                out_dir      = out_dir,
                iterations   = iterations
        )
        
        # 2) extract stats from everything just created
        stats_list <- extract_stats(out_dir = out_dir, var_ctrl = var_ctrl, r = r, iterations = iterations)
        
}






sample.sub_tpm <- readRDS("~/GTEx_4tissue_lap_samples284_TPM.rds")

# remove low variance eRNAs from all 4 tissues simultaneously
sd_table <- matrix(nrow = nrow(sample.sub_tpm[[1]]), ncol = length(sample.sub_tpm))
for (i in seq_along(sample.sub_tpm)) {
        row_sd <- apply(sample.sub_tpm[[i]], 1, sd)
        sd_table[, i] <- row_sd
}
rownames(sd_table) <- rownames(sample.sub_tpm[[1]])
colnames(sd_table) <- names(sample.sub_tpm)

sum_sd  <- apply(sd_table, 1, sum)
row_keep <- sum_sd > 5

filtered_data <- lapply(sample.sub_tpm, function(tissue) {
        tissue[row_keep, , drop = FALSE]
})





run_and_extract(data.list = filtered_data,
                var_ctrl = 2, 
                r  = 0.9,
                out_dir = "~tech_val/noisy_cf/",
                iterations = 100)
























