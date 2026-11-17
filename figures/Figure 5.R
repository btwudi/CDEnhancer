library(survival)
library(survminer)
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(reshape2)


##--- Figure 5A ---##
tcga <- readRDS("~/tcga.brain.cancer.data_514_updated_123022.rds")

myData <- data.frame(
        time   = tcga$ph$PFI.time,
        status = tcga$ph$PFI,
        group  = tcga$metadata$Label,
        stringsAsFactors = FALSE
)

# Set group order to match the legend and colors
myData$group <- factor(
        myData$group,
        levels = c("GBM", "IDH-A", "IDH-O")
)

# Retain samples with complete survival and group information
myData <- myData[
        complete.cases(myData[, c("time", "status", "group")]),
]

fit <- survfit(
        Surv(time, status) ~ group,
        data = myData
)

p <- ggsurvplot(
        fit,
        data = myData,
        surv.median.line = "none",
        
        legend.title = "",
        legend.labs = c(
                "Glioblastoma",
                "Astrocytoma",
                "Oligodendroglioma"
        ),
        legend = c(0.80, 0.70),
        
        xlab = "Time (days)",
        ylab = "PFI probability",
        
        pval = TRUE,
        pval.method = FALSE,
        pval.size = 5,
        pval.method.size = 4,
        pval.coord = c(50, 0.05),
        
        risk.table = TRUE,
        risk.table.title = "Number at risk",
        risk.table.height = 0.25,
        risk.table.y.text.col = TRUE,
        tables.theme = theme_cleantable(),
        
        palette = c("#E7B800", "#2E9FDF", "#b800e7"),
        ggtheme = theme_bw(base_size = 14),
        
        font.x = c(16, "plain", "black"),
        font.y = c(16, "plain", "black"),
        font.tickslab = c(14, "plain", "black"),
        font.legend = c(12, "plain", "black")
)

p$plot <- p$plot +
        ggplot2::theme(
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                
                legend.position = "inside",
                legend.position.inside = c(0.98, 0.98),
                legend.justification = c(1, 1),
                
                legend.background = ggplot2::element_rect(
                        fill = "white", colour = NA
                )
        ) +
        ggplot2::guides(
                colour = ggplot2::guide_legend(ncol = 1)
        )

ggsave("~/Figure5A_left.png", plot = p$plot, width = 4.5, height = 4, dpi = 300)



# plot survival between IDH WT and IDH Mut
# Define IDH groups
labels <- as.character(tcga$metadata$Label)
idh_group <- rep(NA_character_, length(labels))
idh_group[labels == "GBM"] <- "IDH-wild-type"
idh_group[labels %in% c("IDH-A", "IDH-O")] <- "IDH-mutant"

# Assumes metadata and survival records have matching sample order
myData_IDH <- data.frame(
        time   = tcga$ph$PFI.time,
        status = tcga$ph$PFI,
        group  = factor(
                idh_group,
                levels = c("IDH-wild-type", "IDH-mutant")
        )
)

myData_IDH <- myData_IDH[
        complete.cases(myData_IDH),
        ,
        drop = FALSE
]

fit_IDH <- survfit(
        Surv(time, status) ~ group,
        data = myData_IDH
)

p_IDH <- ggsurvplot(
        fit_IDH,
        data = myData_IDH,
        surv.median.line = "none",
        
        legend.title = "",
        legend.labs = c("IDH-wild-type", "IDH-mutant"),
        
        xlab = "Time (days)",
        ylab = "PFI probability",
        
        pval = TRUE,
        pval.method = FALSE,
        pval.coord = c(100, 0.15),
        pval.size = 5,
        
        risk.table = TRUE,
        risk.table.title = "Number at risk",
        risk.table.height = 0.25,
        risk.table.y.text.col = TRUE,
        tables.theme = theme_cleantable(),
        
        palette = c("#E7B800", "#2E9FDF"),
        ggtheme = theme_bw(base_size = 14),
        
        font.x = c(16, "plain", "black"),
        font.y = c(16, "plain", "black"),
        font.tickslab = c(14, "plain", "black"),
        font.legend = c(12, "plain", "black")
)

p_IDH$plot <- p_IDH$plot +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "inside",
                legend.position.inside = c(0.98, 0.98),
                legend.justification = c(1, 1),
                legend.background = element_rect(
                        fill = "white", colour = NA
                )
        ) +
        guides(colour = guide_legend(ncol = 1))

ggsave("~/Figure5A_right.png", plot = p_IDH$plot, width = 4.5, height = 4, dpi = 300)





##--- Figure 5C ---##
glio <- readRDS("~/inhouse_idh_wt_mut_merged_n_processed.RDS")
glio.mtx <- GetAssayData(object = glio, assay = "RNA", layer = "data")
glio.mtx <- as.matrix(glio.mtx)
glio.meta <- glio@meta.data

meta <- glio.meta[glio.meta$Anno != "Unknown", ]
meta$Anno_new <- str_replace_all(meta$Anno, "B cell|T cell", "Immune")
meta$Anno_new <- str_replace_all(meta$Anno_new, "Endothelial|Stromal", "Endo/Stromal")
meta$Anno_new <- str_replace_all(meta$Anno_new, "Immune", "Lymphocyte")
meta$Anno_new <- str_replace_all(meta$Anno_new, "Oligos", "Oligo")

glio_sub <- subset(glio, cells = rownames(meta))
glio_sub@meta.data$Anno_new <- meta$Anno_new

Fig5C <- DimPlot(glio_sub, reduction = "umap", label = TRUE, label.size = 6, group.by = "Anno_new", alpha = 0.2) +
        labs(title = "Clustering of in-house glioma single cells") +
        scale_color_brewer(palette = "Set2") +
        theme(legend.position = "none")

ggsave("~/Figure5C.png", plot = Fig5C, width = 6, height = 6, dpi = 600)





##--- Figure 5D ---##
de_eRNA_table <- read.csv("~/DE_eRNA_TCGA558new_glioma453_log2FC_0.5_qv_0.2_3confounding.csv", row.names = 1)
de_eRNA_table2 <- data.frame(type=factor(de_eRNA_table$cell_type), group=factor(de_eRNA_table$direction))

dat <- dcast(de_eRNA_table2, type ~ group, fun.aggregate = length)
dat.melt <- reshape2::melt(dat, id.vars = "type", measure.vars = c("IDH-WT", "IDH-Mut"))
dat.melt

Fig5D <- ggplot(dat.melt, aes(x = type,y = value, fill = variable)) + 
        geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
        
        geom_text(aes(label = value), position = position_dodge(width = .8), vjust = -0.5, size = 4) +
        scale_fill_manual(values = c("#fdc086", "#a6cee3")) + 
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

ggsave(filename = "~/Figure5D.png",
       plot = Fig5D, width = 6, height = 4, dpi = 300)




##--- Figure 5E ---##
plot_long_ih_eRNA <- readRDS("/Users/wud8/Documents/GPS Lab/deconvolution/eRNA_deconvolution/CDEnhancer_summation/CDEnhancer_summation_k-1_cf_Application_Gliomas/inhouse/TCGA558new_glioma453_rerun_eRNA_score_in_inhouse_all.rds")

p_df_ih_eRNA <- plot_long_ih_eRNA %>%
        group_by(signature) %>%
        summarise(
                p_value = case_when(
                        str_detect(signature[1], regex("WT", ignore_case = TRUE)) ~
                                wilcox.test(
                                        score[IDH == "IDH_WT"],
                                        score[IDH == "IDH_Mut"],
                                        alternative = "greater"
                                )$p.value,
                        str_detect(signature[1], regex("Mut", ignore_case = TRUE)) ~
                                wilcox.test(
                                        score[IDH == "IDH_Mut"],
                                        score[IDH == "IDH_WT"],
                                        alternative = "greater"
                                )$p.value,
                        TRUE ~ NA_real_
                ),
                y_pos = max(score, na.rm = TRUE) + 0.2,
                .groups = "drop"
        ) %>%
        mutate(
                label = case_when(
                        p_value <= 0.0001 ~ "****",
                        p_value <= 0.001 ~ "***",
                        p_value <= 0.01 ~ "**",
                        p_value <= 0.05 ~ "*",
                        TRUE ~ "ns"
                )
        )

p <- ggplot(plot_long_ih_eRNA, aes(x = signature, y = score, fill = IDH)) +
        geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
        geom_text(data = p_df_ih_eRNA,
                  aes(x = signature, y = y_pos, label = label),
                  inherit.aes = FALSE,
                  size = 4) +
        #ylim(0, 3) +
        scale_fill_manual(values = c("#fdc086", "#a6cee3")) +
        #scale_fill_paletteer_d("MetBrewer::Austria") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold", size = 10, color = "black"),
              axis.text.y = element_text(face = "bold", size = 14, color = "black"),
              axis.title.y = element_text(face = "bold", size = 14),
              axis.ticks.y = element_line(color = "black"),  # Ensure y-axis ticks are displayed
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              panel.border = element_blank(),      # Remove panel border
              axis.line = element_line(color = "black", size = 1),  # Keep x and y axis lines
              plot.title = element_text(hjust = 0.5)
        ) +
        guides(fill = guide_legend(title = NULL))

ggsave("~/Figure5D.png", plot = p, width = 6, height = 4, dpi = 300)






##--- Figure 5F ---##
random_AUC_long <- readRDS("/Users/wud8/Documents/GPS Lab/deconvolution/eRNA_deconvolution/CDEnhancer_summation/CDEnhancer_summation_k-1_cf_Application_Gliomas/inhouse/TCGA558new_glioma453_rerun_eRNA_score_in_inhouse_random_AUC.rds")
observed_AUC_df <- readRDS("/Users/wud8/Documents/GPS Lab/deconvolution/eRNA_deconvolution/CDEnhancer_summation/CDEnhancer_summation_k-1_cf_Application_Gliomas/inhouse/TCGA558new_glioma453_rerun_eRNA_score_in_inhouse_observed_AUC_df.rds")

p <- ggplot(random_AUC_long, aes(x = signature, y = random_AUC)) +
        geom_violin(aes(fill = "Random AUC"), color = "grey50", scale = "width") +
        geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, show.legend = FALSE) +
        geom_point(data = observed_AUC_df, aes(x = signature, y = observed_AUC, color = "Observed AUC"), inherit.aes = FALSE, size = 3) +
        geom_text(data = observed_AUC_df, aes(x = signature, y = 1.03, label = significance), inherit.aes = FALSE, color = "black",  size = 4, show.legend = FALSE) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "steelblue", show.legend = FALSE) +
        coord_flip(clip = "off") +
        scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, by = 0.2)) +
        scale_fill_manual(name = NULL, values = c("Random AUC" = "grey85")) +
        scale_color_manual(name = NULL, values = c("Observed AUC" = "red")) +
        labs(x = NULL, y = "AUC", title = "Observed AUC versus random eRNA signatures") +
        theme_classic() +
        theme(legend.position = "right",
              legend.text = element_text(size = 8),
              plot.title = element_text(hjust = 0.5)) +
        guides(color = guide_legend(order = 1), fill = guide_legend(order = 2))

ggsave("~/Figure5E.png", plot = p, width = 5, height = 4, dpi = 300)




##--- Figure 5G ---##
tsne_plot <- readRDS(
        "/Users/wud8/Documents/GPS Lab/deconvolution/eRNA_deconvolution/CDEnhancer_summation/CDEnhancer_summation_k-1_cf_Application_Gliomas/inhouse/TCGA558new_glioma453_rerun_eRNA_malignant_stratify_inhouse_TSNE.rds"
)

tsne.r_plot <- readRDS(
        "/Users/wud8/Documents/GPS Lab/deconvolution/eRNA_deconvolution/CDEnhancer_summation/CDEnhancer_summation_k-1_cf_Application_Gliomas/inhouse/TCGA558new_glioma453_rerun_eRNA_malignant_random_stratify_inhouse_TSNE.rds"
)

colors <- c("IDH_WT" = "#fc993b", "IDH_Mut" = "#6caed1")

p1 <- ggplot(tsne_plot) + 
        geom_point(aes(x = tSNE_1, y = tSNE_2, fill = cancer_type, color = cancer_type), 
                   shape = 21, size = 3, stroke = 1, alpha = 0.7) +
        scale_color_manual(values = c("#fc993b", "#6caed1")) +
        scale_fill_manual(values = c("#fc993b", "#6caed1")) +
        labs(title = "Malignant") +
        theme_minimal() +
        theme(
                legend.title = element_blank(),   # Hides the legend title
                legend.text = element_text(size = 12),
                panel.grid.major = element_blank(),  # Remove major grid lines
                panel.grid.minor = element_blank(),  # Remove minor grid lines
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
        )

p2 <- ggplot(tsne.r_plot) + 
        geom_point(aes(x = tSNE_1, y = tSNE_2, fill = cancer_type, color = cancer_type), 
                   shape = 21, size = 3, stroke = 1, alpha = 0.7) +
        scale_color_manual(values = c("#fc993b", "#6caed1")) +
        scale_fill_manual(values = c("#fc993b", "#6caed1")) +
        labs(title = "Random") +
        theme_minimal() +
        theme(
                legend.title = element_blank(),   # Hides the legend title
                legend.text = element_text(size = 12),
                panel.grid.major = element_blank(),  # Remove major grid lines
                panel.grid.minor = element_blank(),  # Remove minor grid lines
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
        )

combined_plot <- (p1 | p2) +
        plot_layout(guides = "collect") 

ggsave("~/Figure5G.png", plot = combined_plot, width = 10, height = 4.5, dpi = 300)








