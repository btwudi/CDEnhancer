## Plotting code for Figure panels A, B, D, and E.
## A/B use the true-fraction results; D/E use the noisy-fraction results.
## Each panel shows bidirectional simulations above unidirectional simulations.

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(grid)

##--- Figure 3ABDE ---##
stats_true_rds <- path.expand("~/extracted_stats_var_ctrl_2it100.RDS")
stats_noisy_rds <- path.expand("~/extracted_stats_var_ctrl_2target_r0.9it100.RDS")

out_dir <- path.expand("~/Figure3")

# Helpers retained from the source plotting script.
mean_every_4_rows <- function(df) {
        grp <- (seq_len(nrow(df)) - 1) %/% 4
        blocks <- split(df, grp)
        
        out <- lapply(blocks, function(x) {
                means <- colMeans(
                        x[, vapply(x, is.numeric, logical(1)) &
                                  !(names(x) %in% c("snr", "num")), drop = FALSE],
                        na.rm = TRUE
                )
                
                data.frame(
                        snr = x$snr[1],
                        num = x$num[1],
                        direction = x$direction[1],
                        as.list(means),
                        check.names = FALSE
                )
        })
        
        out <- do.call(rbind, out)
        rownames(out) <- NULL
        out
}

add_group <- function(df,
                      group_levels = c("CDEnhancer", "CellDMC", "TOAST")) {
        df$group <- ifelse(
                grepl("CellDMC", rownames(df)), "CellDMC",
                ifelse(grepl("TOAST", rownames(df)), "TOAST", "CDEnhancer")
        )
        df$group <- factor(df$group, levels = group_levels)
        df
}

collect_metric <- function(lst, cols_keep, complete = FALSE) {
        out <- do.call(rbind, lapply(lst, function(df) {
                missing_cols <- setdiff(cols_keep, colnames(df))
                if (length(missing_cols) > 0) df[missing_cols] <- NA
                df[, cols_keep, drop = FALSE]
        }))
        if (complete) out <- out[complete.cases(out), , drop = FALSE]
        add_group(out)
}

make_bi_metric <- function(df, up_col, down_col, value_name, dir_name) {
        out <- rbind.data.frame(df, df)
        out[[dir_name]] <- factor(
                c(rep("up", nrow(df)), rep("down", nrow(df))),
                levels = c("down", "up")
        )
        out[[value_name]] <- c(df[[up_col]], df[[down_col]])
        out
}

base_plot_theme <- function(axis_text_size = 8, axis_title_size = 10,
                            strip_size = 8, axis_line_size = 0.3) {
        theme_minimal() +
                theme(
                        axis.text.x = element_text(
                                face = "bold", size = axis_text_size, color = "black"
                        ),
                        axis.title.x = element_text(face = "bold", size = axis_title_size),
                        axis.text.y = element_text(
                                face = "bold", size = axis_text_size, color = "black"
                        ),
                        axis.title.y = element_text(face = "bold", size = axis_title_size),
                        axis.ticks.y = element_line(color = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.spacing.y = unit(1.2, "lines"),
                        axis.line.x = element_blank(),
                        axis.line = element_line(color = "black", linewidth = axis_line_size),
                        plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
                        legend.position = "right",
                        legend.text = element_text(size = 8, color = "black"),
                        legend.title = element_text(size = 8, face = "bold", color = "black"),
                        strip.text = element_text(size = strip_size, face = "bold", vjust = 2)
                )
}

prepare_mean_metrics <- function(stats_rds) {
        stats_list <- readRDS(stats_rds)
        mean_list <- lapply(stats_list, mean_every_4_rows)
        
        mean_AUC_bi <- collect_metric(
                mean_list,
                c("AUC_up", "AUC_down", "snr", "num", "direction"),
                complete = TRUE
        )
        mean_AUC_uni <- collect_metric(
                mean_list,
                c("AUC", "snr", "num", "direction"),
                complete = TRUE
        )
        mean_PRAUC_bi <- collect_metric(
                mean_list,
                c("PR_AUC_up", "PR_AUC_down", "snr", "num", "direction"),
                complete = TRUE
        )
        mean_PRAUC_uni <- collect_metric(
                mean_list,
                c("PR_AUC", "snr", "num", "direction"),
                complete = TRUE
        )
        
        list(
                AUC_bi = make_bi_metric(
                        mean_AUC_bi, "AUC_up", "AUC_down", "AUC", "AUC_dir"
                ),
                AUC_uni = mean_AUC_uni,
                PRAUC_bi = make_bi_metric(
                        mean_PRAUC_bi, "PR_AUC_up", "PR_AUC_down", "PR_AUC", "PR_AUC_dir"
                ),
                PRAUC_uni = mean_PRAUC_uni
        )
}

# Top plot: bidirectional simulations, faceted by number of altered eRNAs.
plot_bi_metric <- function(df, value_col, fill_col, y_lab, title) {
        ggplot(df, aes(x = factor(snr), y = .data[[value_col]],
                       fill = .data[[fill_col]])) +
                stat_boxplot(
                        geom = "errorbar", width = 0.6,
                        position = position_dodge(width = 0.7),
                        show.legend = FALSE
                ) +
                geom_jitter(
                        position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2),
                        shape = 21, stroke = 0.1, size = 1, alpha = 0.7,
                        show.legend = TRUE
                ) +
                geom_boxplot(
                        width = 0.6, position = position_dodge(width = 0.7),
                        outlier.shape = NA, linewidth = 0.3, show.legend = FALSE
                ) +
                facet_wrap(~ num) +
                scale_fill_brewer(palette = "Accent") +
                labs(title = title, x = "SNR", y = y_lab) +
                scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
                geom_hline(yintercept = 0, color = "black", linewidth = 1) +
                base_plot_theme(10, 12, 10, 0.5) +
                geom_vline(
                        xintercept = c(1.5, 2.5, 3.5), linetype = "dotted",
                        color = "gray", linewidth = 0.25
                ) +
                guides(fill = guide_legend(
                        title = NULL, override.aes = list(shape = 21, size = 3, alpha = 1)
                ))
}

# Bottom plot: unidirectional simulations, faceted by number of altered eRNAs.
plot_uni_metric <- function(df, value_col, y_lab, title) {
        df$direction <- factor(df$direction, levels = c("down", "up"))
        
        ggplot(df, aes(x = factor(snr), y = .data[[value_col]], fill = direction)) +
                stat_boxplot(
                        geom = "errorbar", width = 0.6,
                        position = position_dodge(width = 0.7),
                        show.legend = FALSE
                ) +
                geom_jitter(
                        position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2),
                        shape = 21, stroke = 0.1, size = 1, alpha = 0.7,
                        show.legend = TRUE
                ) +
                geom_boxplot(
                        width = 0.6, position = position_dodge(width = 0.7),
                        outlier.shape = NA, linewidth = 0.3, show.legend = FALSE
                ) +
                facet_grid(. ~ num) +
                scale_fill_brewer(palette = "Accent") +
                labs(title = title, x = "SNR", y = y_lab) +
                scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
                geom_hline(yintercept = 0, color = "black", linewidth = 1) +
                base_plot_theme(10, 12, 10, 0.5) +
                geom_vline(
                        xintercept = c(1.5, 2.5, 3.5), linetype = "dotted",
                        color = "gray", linewidth = 0.25
                ) +
                guides(fill = guide_legend(
                        title = NULL, override.aes = list(shape = 21, size = 3, alpha = 1)
                ))
}

make_auc_panel <- function(metrics, panel_tag) {
        bi <- metrics$AUC_bi[metrics$AUC_bi$group == "CDEnhancer", , drop = FALSE]
        uni <- metrics$AUC_uni[metrics$AUC_uni$group == "CDEnhancer", , drop = FALSE]
        
        top <- plot_bi_metric(bi, "AUC", "AUC_dir", "AUC", "bi direction") +
                labs(tag = panel_tag)
        bottom <- plot_uni_metric(uni, "AUC", "AUC", "uni direction")
        top / bottom + plot_layout(heights = c(1, 1))
}

make_prauc_panel <- function(metrics, panel_tag) {
        bi <- metrics$PRAUC_bi[
                metrics$PRAUC_bi$group == "CDEnhancer", , drop = FALSE
        ]
        uni <- metrics$PRAUC_uni[
                metrics$PRAUC_uni$group == "CDEnhancer", , drop = FALSE
        ]
        
        top <- plot_bi_metric(
                bi, "PR_AUC", "PR_AUC_dir", "PR_AUC", "bi direction"
        ) + labs(tag = panel_tag)
        bottom <- plot_uni_metric(uni, "PR_AUC", "PR_AUC", "uni direction")
        top / bottom + plot_layout(heights = c(1, 1))
}

# Panels A/B: true cell fractions.
true_metrics <- prepare_mean_metrics(stats_true_rds)
panel_A <- make_auc_panel(true_metrics, "A")
panel_B <- make_prauc_panel(true_metrics, "B")

# Panels D/E: noisy estimated cell fractions.
noisy_metrics <- prepare_mean_metrics(stats_noisy_rds)
panel_D <- make_auc_panel(noisy_metrics, "D")
panel_E <- make_prauc_panel(noisy_metrics, "E")

# Save each panel and the combined A/B/D/E layout.
ggsave(file.path(out_dir, "Figure3A_ROC_AUC_true_fractions.png"), panel_A,
       width = 5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "Figure3B_PR_AUC_true_fractions.png"), panel_B,
       width = 5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "Figure3D_ROC_AUC_noisy_fractions.png"), panel_D,
       width = 5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "Figure3E_PR_AUC_noisy_fractions.png"), panel_E,
       width = 5, height = 6, dpi = 300)





##--- Figure 3C ---##
cfs = readRDS("~/generated_cf_and_mimic_predicted_cf_var_ctrl_2_target_r_0.9.RDS")

for (j in 1:4)
{
        spearman_correlation <- cor(cfs$true[,paste0("V", j)], cfs$pred[,paste0("V", j)], method = "spearman")
        
        x <- cfs$true[,paste0("V", j)]
        y <- cfs$pred[,paste0("V", j)]
        data <- data.frame(x, y)
        
        p <- ggplot(data, aes(x = x, y = y)) +
                geom_point(color = "#bba6e3", size = 2, alpha = 0.8) +  # Add scatter points
                geom_smooth(method = "lm", color = "#404040", fill = "#efcde3", se = TRUE) +  # Add the slope line with shaded confidence interval
                
                # Add slope and p-value text at the top left corner
                annotate("text", x = min(data$x), y = max(data$y) * 0.8, hjust = 0, vjust = -1.5, 
                         label = paste("Spearman rho =", round(spearman_correlation, 2)), size = 5, color = "black") +
                labs(title = paste0("\nCorrelation: cell_type ", j),
                     x = "Ground truth fraction",
                     y = "Imputed cell cell fraction") +
                # Customize theme and position of elements
                theme_minimal() +
                theme(
                        plot.title = element_text(hjust = 0.5, size = 16),
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 14),
                        axis.line = element_line(color = "black", size = 1),  # Keep x and y axis lines
                        panel.grid.major = element_blank(),  # Remove major grid lines
                        panel.grid.minor = element_blank(),  # Remove minor grid lines
                        axis.ticks.x = element_line(color = "black"),  
                        axis.ticks.y = element_line(color = "black"),  
                        
                        legend.position = "none"
                )
        
        path <- paste(paste("~", paste0("Figure3C_cell_type_", j), sep = "/"), ".png", sep = "")
        ggsave(path, plot = p, width = 5, height = 4, dpi = 300)
        
}






##--- Figure 3E ---##
data <- readRDS("~/prAUC_across_cell_types_num_3000_snr_3_corrected.rds")

cell_types <- paste("cell type", 1:4)
items_to_plot <- c("up", "down", "bi")

df_list <- lapply(items_to_plot, function(item) {
        if (length(data[[item]]) != 16) {
                stop("Expected 16 values for list item: ", item)
        }
        
        mat <- matrix(data[[item]], nrow = 4, ncol = 4)
        dimnames(mat) <- list(cell_types, cell_types)
        
        df <- as.data.frame(as.table(mat))
        names(df) <- c("Predicted", "Ground_truth", "value")
        
        df$Predicted <- factor(df$Predicted, levels = rev(cell_types))
        df$Ground_truth <- factor(df$Ground_truth, levels = cell_types)
        df$item <- item
        df$percent_value <- sprintf("%.2f%%", df$value * 100)
        df
})

plot_df <- do.call(rbind, df_list)
plot_df$item <- factor(plot_df$item, levels = items_to_plot)

p <- ggplot(plot_df, aes(x = Ground_truth, y = Predicted, fill = value)) +
        geom_tile() +
        geom_text(aes(label = percent_value), color = "black", size = 4) +
        facet_wrap(~item, nrow = 1) +
        scale_fill_gradient(
                low = "white",
                high = "red",
                labels = scales::label_percent(accuracy = 1),
                limits = c(0, 1)
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        labs(
                title = "Cross-cell-type PR-AUC",
                x = "Ground truth",
                y = "Predicted",
                fill = "PR-AUC"
        ) +
        theme_minimal() +
        theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                axis.text.y = element_text(size = 12, color = "black"),
                axis.title.x = element_text(face = "bold", size = 14),
                axis.title.y = element_text(face = "bold", size = 14),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 10),
                strip.text = element_text(face = "bold", size = 12),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                panel.spacing.x = grid::unit(1.2, "cm"),
                plot.title = element_text(hjust = 0.5)
        )

ggsave("~/Figure3E.png", plot = p, width = 12, height = 4.5, dpi = 300)


