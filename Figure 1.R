library(ggplot2)
library(dplyr)

##--- Figure 1C ---##
files <- list.files("~/one_caner_type_function", pattern = "\\.rds$", full.names = TRUE)
files.spl <- as.data.frame(t(sapply(1:length(files), function(i) unlist(strsplit(files[i], "function/|_eRNA15k")))))    

result <- lapply(files, function(f) {
        dat = readRDS(f)
        row_means = rowMeans(dat)
        
        spl = unlist(strsplit(f, "function/|_eRNA15k"))
        cor = cbind.data.frame(gene = rownames(dat), correlation = row_means, type = spl[2])
        return(cor)
})   
names(result) <- files.spl$V2
final.cor <- do.call(rbind, result)        
final.cor$type = reorder(final.cor$type, final.cor$correlation, median, decreasing = TRUE) # sort types by median correlation

Fig1C <- ggplot(final.cor, aes(x = type, y = correlation)) +
        geom_violin(trim=FALSE, fill = "lightblue") +
        geom_boxplot(width = 0.1, outlier.shape = 16, outlier.size = 1) +
        geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
        theme_classic() +
        theme(
                axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
                axis.title.y = element_text(face = "bold")
        ) +
        labs(x = NULL, y = "Correlation")
            
ggsave("~/Figure1C.png", plot = Fig1C, width = 10, height = 3, dpi = 300)      
   




##--- Figure 1D ---##
final.cor_good <- final.cor[final.cor$correlation >= 0.4, ]
good_pred_gene_nums <- table(final.cor_good$type)
vari_gene_nums <- readRDS("~/vari_gene_nums.rds")
vari_gene_nums <- vari_gene_nums[match(names(good_pred_gene_nums), names(vari_gene_nums))] 
good_pct <- (round(good_pred_gene_nums/vari_gene_nums, 3))*100

pct_df <- data.frame(
        type = factor(names(good_pct), levels = names(good_pct)),
        percentage = as.numeric(good_pct)
)

Fig1D <- ggplot(pct_df, aes(x = type, y = percentage)) +
        geom_col(fill = "lightyellow", color = "black", width = 0.7) +
        geom_text(aes(label = percentage), vjust = -0.5, size = 2.5, fontface = "bold") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
              axis.title.y = element_text(face = "bold")) +
        labs(x = NULL, y = "Percentage") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("~/Figure1D.png", plot = Fig1D, width = 10, height = 2, dpi = 300)  




