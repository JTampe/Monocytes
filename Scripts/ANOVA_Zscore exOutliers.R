folder <- "ANOVA_Zscore exOutliers"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric(),
                      TP1_P_Value = numeric(), TP2_P_Value = numeric(), TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (j in 1:length(unique(data_rE_mean$Gene))) {
    gene <- data_rE_mean$Gene[j]
    df_gene <- filter(data_rE_mean, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Category), !is.na(mean_Z), !is.na(Subpopulation), !is.na(Timepoint))
    
    for (i in 1:length(unique(df_gene$Subpopulation))) {
        sup <- df_gene$Subpopulation[i]
        df <- filter(df_gene, Subpopulation == sup)
        if (nrow(df) <= 0 || gene %in% exclude_genes) {
            cat("Skipping:", gene, "for subpopulation:", sup, "due to insufficient data or in excluded genes list\n")
            next  # Skip to the next iteration of the loop
        }
        
        df <- df[!df$mean_Z %in% boxplot.stats(df$mean_Z)$out, ]
        df$Timepoint <- factor(as.vector(df$Timepoint))
        baseline <- filter(df, Timepoint == "TP0")
        median_TP0 <- median(baseline$mean_Z)
        
        # Perform ANOVA
        anova_result <- aov(mean_Z ~ Timepoint + Error(SampleID/Timepoint), data = df)
        
        # Perform Tukey's post-hoc test
        tukey_result <- TukeyHSD(anova_result, "Timepoint")
        
        # Extract p-values for each timepoint comparison against TP0
        tp1_p_value <- tukey_result$Timepoint["TP1-TP0", "p adj"]
        tp2_p_value <- tukey_result$Timepoint["TP2-TP0", "p adj"]
        tp3_p_value <- tukey_result$Timepoint["TP3-TP0", "p adj"]
        tp4_p_value <- tukey_result$Timepoint["TP4-TP0", "p adj"]
        
        # Save results to data frame
        results <- rbind(results, data.frame(Gene = gene, Subpopulation = sup, 
                                             P_Value = summary(anova_result)[[1]]$"Pr(>F)"[1],
                                             TP1_P_Value = tp1_p_value, TP2_P_Value = tp2_p_value, 
                                             TP3_P_Value = tp3_p_value, TP4_P_Value = tp4_p_value))
        
        if (plot_save == "y") {
            my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
            
            ggplot(df, aes(x = Timepoint, y = mean_Z, color = Timepoint)) + 
                geom_boxplot(outlier.shape = NA, width = 0.8, position = position_dodge(0.8)) + 
                geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
                geom_hline(yintercept = median_TP0, linetype = 2) +
                scale_color_manual(values = my_colors) + 
                ylab(paste(gene, "expression relative to B2M")) +
                theme_minimal() +
                theme(legend.position = "none")
            
            # Add ANOVA label
            max_y <- max(df$mean_Z, na.rm = TRUE)
            p_value <- summary(anova_result)[[1]]$"Pr(>F)"[1]
            p_label <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))
            annotate("text", x = 2.5, y = max_y + 0.1 * max_y, label = paste("ANOVA p =", round(p_value, 3), p_label), size = 5)
            
            # Add Tukey's post-hoc results
            for (tp in c("TP1", "TP2", "TP3", "TP4")) {
                tp_p_value <- get(paste0(tp, "_p_value"))
                if (!is.na(tp_p_value)) {
                    tp_label <- ifelse(tp_p_value < 0.001, "***", ifelse(tp_p_value < 0.01, "**", ifelse(tp_p_value < 0.05, "*", "ns")))
                    annotate("text", x = which(levels(df$Timepoint) == tp), 
                             y = max_y + 0.1 * max_y, 
                             label = paste(tp, "vs TP0:", round(tp_p_value, 3), tp_label), 
                             size = 4, vjust = -0.5)
                }
            }
            
            file_name <- file.path(folder, paste(gene, "_", sup, "_Zscore_ANOVA.png", sep = ""))
            ggsave(filename = file_name, plot = last_plot())
        }
    }
}
ANOVA_Zscore <- results
ANOVA_Zscore$Significance <- sapply(ANOVA_Zscore$P_Value, get_significance)
ANOVA_Zscore$TP1_Significance <- sapply(ANOVA_Zscore$TP1_P_Value, get_significance)
ANOVA_Zscore$TP2_Significance <- sapply(ANOVA_Zscore$TP2_P_Value, get_significance)
ANOVA_Zscore$TP3_Significance <- sapply(ANOVA_Zscore$TP3_P_Value, get_significance)
ANOVA_Zscore$TP4_Significance <- sapply(ANOVA_Zscore$TP4_P_Value, get_significance)
write.csv(ANOVA_Zscore, file = "ANOVA_results/ANOVA_Results_Zscore.csv", row.names = FALSE)
