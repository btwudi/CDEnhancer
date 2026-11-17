library(ggplot2)
library(reshape2)
library(paletteer)
library(dplyr)
library(tibble)
library(purrr)
library(readxl)
library(tidyr)
library(patchwork)


##--- Figure 6A ---##
de_eRNA_table <- read.csv("~/DE_eRNAs_TCGA_pancancer_hot_cold_tumor_rerun_q01_logFC05.csv", row.names = 1)
de_eRNA_table2 <- data.frame(type=factor(de_eRNA_table$cell_type), group=factor(de_eRNA_table$direction))

dat <- reshape2::dcast(de_eRNA_table2, type ~ group, fun.aggregate = length)
dat.melt <- reshape2::melt(dat, id.vars = "type", measure.vars = c("Hot Tumors", "Cold Tumors"))
dat.melt

p <- ggplot(dat.melt, aes(x = type,y = value, fill = variable)) + 
        geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
        
        geom_text(aes(label = value), position = position_dodge(width = .8), vjust = -0.5, size = 4) +
        #scale_fill_manual(values = c("#fdc086", "#a6cee3")) + 
        scale_fill_paletteer_d("fishualize::Acanthurus_olivaceus") +
        labs(x = NULL,
             y = "Differentially expressed\neRNAs numbers") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14, color = "black"),
              axis.text.y = element_text(face = "bold", size = 14, color = "black"),
              axis.title.y = element_text(face = "bold", size = 14),
              axis.ticks.y = element_line(color = "black"),  # Ensure y-axis ticks are displayed
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              panel.border = element_blank(),      # Remove panel border
              axis.line = element_line(color = "black", size = 1),  # Keep x and y axis lines
              plot.title = element_text(hjust = 0.5),
              legend.text = element_text(size = 10),
        ) +
        guides(fill = guide_legend(title = NULL))

ggsave("~/Figure6A.png", plot = p, width = 6, height = 5, dpi = 300)




##--- Figure 6B ---##
cor_list <- readRDS("~/CD8T_fraction_association_with_DE_eRNAs_and_closest_gene_signature_TCGA_pancancer_hot_cold_tumor_q01_logFC05.RDS")

cor_df <- purrr::imap_dfr(cor_list, ~ {
        .x %>%
                rownames_to_column("dataset") %>%
                mutate(panel = .y)
})
cor_df$panel <- factor(
        cor_df$panel,
        levels = c("CT eRNA", "CT closest genes", "Bulk eRNA", "Bulk closest genes")
)

cor_df$stars <- dplyr::case_when(
        cor_df$pv <= 0.0001 ~ "****",
        cor_df$pv <= 0.001  ~ "***",
        cor_df$pv <= 0.01   ~ "**",
        cor_df$pv <= 0.05   ~ "*",
        TRUE ~ "ns"
)

p <- ggplot(cor_df, aes(x = dataset, y = cor, fill = panel)) +
        geom_col(width = 0.7, show.legend = FALSE) +
        geom_text(aes(label = stars), vjust = -0.2, size = 3) +
        facet_wrap(~ panel, ncol = 2, scales = "fixed") +
        labs(x = NULL, y = "Correlation") +
        scale_fill_paletteer_d("rcartocolor::Antique") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold", color = "black"),
              axis.text.y = element_text(size = 9, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold"), 
              axis.title = element_text(size = 12, face = "bold"), 
              legend.text = element_text(size = 10, face = "bold"),
              legend.title = element_blank(),   # Hides the legend title
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.text = element_text(face = "bold", size = 12)
        )

ggsave("~/Figure6B.png", plot = p, width = 4, height = 4, dpi = 300)




##--- Figure 6C ---##
auc_df_eRNA <- readRDS("/Users/wud8/Documents/GPS Lab/deconvolution/eRNA_deconvolution/CDEnhancer_summation/CDEnhancer_summation_k-1_cf_Application_TCGA_hot_cold_tumors/rerun/DE_eRNAs_TCGA32_hot_cold_tumor_eRNA_score_in_ICB_separated_cohorts_timeline.rds")
auc_df_eRNA_sub_cancer = subset(auc_df_eRNA,timeline == "PRE" & signature == "Cancer_Hot Tumors")
p_eRNA <- ggplot(auc_df_eRNA_sub_cancer, aes(x = cohort, y = AUC, fill = cohort)) +
        geom_col(width = 0.7, show.legend = FALSE) +
        geom_text(aes(label = round(AUC, 2)), vjust = -0.3, size = 3) +
        #facet_wrap(~ cohort, ncol = 2, scales = "fixed") +
        labs(title = "eRNAs",x = NULL, y = "AUC") +
        scale_fill_paletteer_d("rcartocolor::Antique") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.05))) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold", color = "black"),
              axis.text.y = element_text(size = 9, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold"), 
              axis.title = element_text(size = 12, face = "bold"), 
              legend.text = element_text(size = 10, face = "bold"),
              legend.title = element_blank(),   # Hides the legend title
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              #panel.border = element_rect(color = "black", fill = NA, size = 1),
              axis.ticks.y = element_line(color = "black"),
              plot.title = element_text(face = "bold", hjust = 0.5),
              panel.border = element_blank(),      # Remove panel border
              axis.line = element_line(color = "black", size = 1)
        )

auc_df <- readRDS("/Users/wud8/Documents/GPS Lab/deconvolution/eRNA_deconvolution/CDEnhancer_summation/CDEnhancer_summation_k-1_cf_Application_TCGA_hot_cold_tumors/rerun/DE_eRNAs_TCGA32_hot_cold_tumor_closest_gene_score_in_ICB_separated_cohorts_timeline.rds")
auc_df_sub_cancer <- subset(auc_df,timeline == "PRE" & signature == "Cancer_Hot Tumors")
p_ge <- ggplot(auc_df_sub_cancer, aes(x = cohort, y = AUC, fill = cohort)) +
        geom_col(width = 0.7, show.legend = FALSE) +
        geom_text(aes(label = round(AUC, 2)), vjust = -0.3, size = 3) +
        #facet_wrap(~ cohort, ncol = 2, scales = "fixed") +
        labs(title = "genes",x = NULL, y = "AUC") +
        scale_fill_paletteer_d("rcartocolor::Antique") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.05))) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold", color = "black"),
              axis.text.y = element_text(size = 9, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold"), 
              axis.title = element_text(size = 12, face = "bold"), 
              legend.text = element_text(size = 10, face = "bold"),
              legend.title = element_blank(),   # Hides the legend title
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              #panel.border = element_rect(color = "black", fill = NA, size = 1),
              axis.ticks.y = element_line(color = "black"),
              plot.title = element_text(face = "bold", hjust = 0.5),
              panel.border = element_blank(),      # Remove panel border
              axis.line = element_line(color = "black", size = 1)
        )

combined_plot <- (p_eRNA | p_ge) +
        plot_layout(guides = "collect") 

ggsave("~/Figure6C.png", plot = combined_plot, width = 4.5, height = 4, dpi = 300)




##--- Figure 6D ---##
auc_df_eRNA_bulk <- readRDS("~/DE_eRNAs_TCGA32_hot_cold_tumor_BULK_eRNA_score_predict_ICB_response_AUC.rds")
auc_df_eRNA_bulk_hot <- auc_df_eRNA_bulk[auc_df_eRNA_bulk$timeline == "PRE" & auc_df_eRNA_bulk$signature == "Hot", ]

p_eRNA_bulk <- ggplot(auc_df_eRNA_bulk_hot, aes(x = cohort, y = AUC, fill = cohort)) +
        geom_col(width = 0.7, show.legend = FALSE) +
        geom_text(aes(label = round(AUC, 2)), vjust = -0.3, size = 3) +
        #facet_wrap(~ cohort, ncol = 2, scales = "fixed") +
        labs(title = "eRNAs",x = NULL, y = "AUC") +
        scale_fill_paletteer_d("rcartocolor::Antique") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.05))) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold", color = "black"),
              axis.text.y = element_text(size = 9, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold"), 
              axis.title = element_text(size = 12, face = "bold"), 
              legend.text = element_text(size = 10, face = "bold"),
              legend.title = element_blank(),   # Hides the legend title
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              #panel.border = element_rect(color = "black", fill = NA, size = 1),
              axis.ticks.y = element_line(color = "black"),
              plot.title = element_text(face = "bold", hjust = 0.5),
              panel.border = element_blank(),      # Remove panel border
              axis.line = element_line(color = "black", size = 1)
        )
        
auc_df_ge_bulk <- readRDS("~/DE_eRNAs_TCGA32_hot_cold_tumor_BULK_gene_score_predict_ICB_response_AUC.rds")
auc_df_ge_bulk_hot <- auc_df_ge_bulk[auc_df_ge_bulk$timeline == "PRE" & auc_df_ge_bulk$signature == "Hot", ]

p_ge_bulk <- ggplot(auc_df_ge_bulk_hot, aes(x = cohort, y = AUC, fill = cohort)) +
        geom_col(width = 0.7, show.legend = FALSE) +
        geom_text(aes(label = round(AUC, 2)), vjust = -0.3, size = 3) +
        #facet_wrap(~ cohort, ncol = 2, scales = "fixed") +
        labs(title = "genes",x = NULL, y = "AUC") +
        scale_fill_paletteer_d("rcartocolor::Antique") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.05))) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold", color = "black"),
              axis.text.y = element_text(size = 9, face = "bold", color = "black"),
              axis.text = element_text(size = 12, face = "bold"), 
              axis.title = element_text(size = 12, face = "bold"), 
              legend.text = element_text(size = 10, face = "bold"),
              legend.title = element_blank(),   # Hides the legend title
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              #panel.border = element_rect(color = "black", fill = NA, size = 1),
              axis.ticks.y = element_line(color = "black"),
              plot.title = element_text(face = "bold", hjust = 0.5),
              panel.border = element_blank(),      # Remove panel border
              axis.line = element_line(color = "black", size = 1)
        )

combined_plot <- (p_eRNA_bulk | p_ge_bulk) +
        plot_layout(guides = "collect") 

ggsave("~/DE_eRNAs_TCGA32_hot_cold_tumor_BULK_eRNA_and_gene_score_predict_ICB_response_AUC_Hot_only.png", plot = combined_plot, width = 4.5, height = 4, dpi = 300)




##--- Figure 6E ---##
summary_table <- read_excel(
        "~/GO_BP_category_assignments_hot_tumor_bulk_vs_cancer_cell.xlsx",
        sheet = "Category summary",
        range = "A5:E16"
)
category_order <- as.character(summary_table[[1]])
plot_table <- data.frame(
        Category = category_order,
        Bulk = as.numeric(summary_table[[3]]) * 100,
        Cancer_cell_specific = as.numeric(summary_table[[5]]) * 100,
        check.names = FALSE
)

plot_long <- pivot_longer(
        plot_table,
        cols = c(Bulk, Cancer_cell_specific),
        names_to = "Analysis",
        values_to = "Percent"
)
plot_long$Category <- factor(plot_long$Category, levels = rev(category_order))
plot_long$Analysis <- factor(
        plot_long$Analysis,
        levels = c("Bulk", "Cancer_cell_specific")
)
plot_long$Label <- ifelse(
        plot_long$Percent == 0,
        "0",
        sprintf("%.1f", plot_long$Percent)
)

p <- ggplot(plot_long, aes(x = Analysis, y = Category, fill = Percent)) +
        geom_tile(color = "white", linewidth = 0.4) +
        geom_text(aes(label = Label), size = 3.3, color = "black") +
        scale_fill_gradient(
                low = "#fff3b0",
                high = "#fb6f92",
                limits = c(0, 40),
                breaks = c(0, 10, 20, 30, 40),
                name = "%"
        ) +
        scale_x_discrete(
                labels = c(
                        Bulk = "Bulk",
                        Cancer_cell_specific = "Cancer-cell-specific"
                ),
                expand = expansion(add = 0)
        ) +
        scale_y_discrete(expand = expansion(add = 0)) +
        labs(x = NULL, y = NULL) +
        theme_classic(base_size = 11) +
        theme(
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_text(angle = 50, hjust = 1, color = "black"),
                axis.text.y = element_text(color = "black"),
                legend.position = "right",
                legend.title = element_text(size = 11),
                legend.text = element_text(size = 10),
                legend.key.height = grid::unit(0.7, "cm"),
                plot.margin = margin(10, 10, 10, 10)
        )

ggsave("~/Figure6E.png", plot = p, width = 4.5, height = 3.5, dpi = 300)




##--- Figure 6F ---##
plot_tb <- readRDS("~/correlation_DEeRNA_and_TRAF1.rds")

ct <- cor.test(plot_tb$eRNA, plot_tb$gene, method = "spearman")
cor_value <- ct$estimate
p_value <- ct$p.value

p <- ggplot(plot_tb, aes(x = eRNA, y = gene)) +
        geom_point(color = "#c2a5cf",size = 1, alpha = 0.6) +
        geom_smooth(method = "lm", color = "#404040", fill = "#69b3a2", se = TRUE) +
        annotate("text", x = 0.3, y = Inf, hjust = 0, vjust = 1.5, 
                 label = paste0("Correlation: ", round(cor_value, 2), "\nP-value = ", formatC(p_value)),
                 size = 4, color = "black") +
        labs(
                x = erna_id,
                y = gene_id
        ) +
        theme_classic()

ggsave("~/Figure6F.png", plot = p, width = 4, height = 4, dpi = 300)
















