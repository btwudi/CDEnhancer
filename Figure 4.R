library(ggplot2)


##--- Figure 4B ---##
results <- readRDS("~/DACR_biological_validation_results.rds")
stat_lap_tb <- results$fisher_enrichment
plot.table <- results$qvalue_plot_data
wilcox_results <- results$wilcoxon_results

fig4B <- ggplot(stat_lap_tb, aes(fill=cancer_type, x = reorder(cancer_type, -odds_ratio), y = odds_ratio)) + 
        geom_bar(position = position_dodge(width = 0.8), fill = "#E5EBFB", stat = "identity", color = "black", width = 0.7) +
        geom_label(aes(label = round(odds_ratio, 2)),  # Add correlation values as text on the bars
                   vjust = 2, size = 3, color = "black", fontface = "bold", 
                   fill = "white", label.size = NA) +  # Set label background to white, remove label border
        
        geom_text(stat='identity', aes(label = label), vjust=-0.4, 
                  position = position_dodge(width = 0.8), size = 5, fontface = "bold") +
        labs(x = NULL,
             y = "odds ratio") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
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
              legend.position = "none"
        ) +
        guides(fill = guide_legend(title = NULL))

ggsave("~/Figure4B.png", plot = fig4B, width = 6, height = 4, dpi = 300)





##--- Figure 4C ---##
fig4C <- ggplot(plot.table, aes(x = status, y = -log10(pmax(qv, 1e-300)), fill = status)) +
        geom_boxplot(width = 0.6, outliers = FALSE, staplewidth = 0.5, show.legend = TRUE) +
        facet_wrap(~ cancer, scales = "free_y", ncol = 7, axes = "all_x", axis.labels = "margins") +
        scale_fill_brewer(palette = "Accent") +
        scale_x_discrete(limits = c("no-overlap", "overlap"), labels = c("No overlap", "Overlap")) +
        scale_y_continuous(expand = expansion(mult = c(0.03, 0.20))) +
        labs(x = NULL, y = expression(-log[10](q)), fill = "status") +
        theme_classic(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
              axis.text.y = element_text(size = 12, color = "black"),
              axis.title.y = element_text(face = "bold"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold", size = 14),
              legend.position = "bottom",
              panel.spacing = grid::unit(1, "lines")) +
        geom_text(data = wilcox_results,
                  aes(x = 1.5, y = Inf, label = label),
                  vjust = 1.3,
                  size = 5,
                  fontface = "bold",
                  inherit.aes = FALSE)

ggsave("~/Figure4C.png", plot = fig4C, width = 6, height = 4, dpi = 300)
