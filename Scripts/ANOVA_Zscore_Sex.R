#### *Male vs Female Analysis* --------------------------------------------------
## you can try to adjust the loop that it also only plots the facet comparisons upon request
folder <- "ANOVA_Zscore_Sex_Plots"
createFolder(folder)  # Ensure the folder exists
results <- data.frame(Sex = character(), Gene = character(), Subpopulation = character(), 
                      P_Value = numeric(), TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                      TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (i in 1:length(unique(data_rE_mean_matched$Gene))) {
    gene <- unique(data_rE_mean_matched$Gene)[i]
    df_gene <- filter(data_rE_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_Z), !is.na(Subpopulation), !is.na(Sex))
    
    for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- unique(df_gene$Subpopulation)[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
        
        # Remove outliers
        df_subpop_clean <- df_subpop[!df_subpop$mean_Z %in% boxplot.stats(df_subpop$mean_Z)$out, ]
        
        if (nrow(df_subpop_clean) <= 0) next  # Skip if no rows remain
        
        tryCatch({
            title_facet <- paste(gene, "expression in", sup, "monocytes (LinReg, rE)", sep=" ")
            
            # Faceted boxplot: comparing Male vs Female across Timepoints
            ggboxplot(df_subpop_clean, 
                      x = "Sex", 
                      y = "mean_Z", 
                      color = "Sex",  
                      title = title_facet,
                      add = "jitter", 
                      facet.by = "Timepoint",
                      legend = "none", 
                      ylab = paste(gene, "expression relative to B2M"), 
                      width = 0.8, 
                      add.params = list(size = 1, alpha = 0.5)) +  
                stat_compare_means(method = "anova", label.y = max(df_subpop_clean$mean_Z)) +        
                stat_compare_means(label = "p.signif", 
                                   method = "wilcox", 
                                   ref.group = "Male", 
                                   hide.ns = TRUE, 
                                   label.y = max(df_subpop_clean$mean_Z)) +
                scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2"))
            
            if (plot_save == "y") {
                file_name_facet <- paste0(folder, "/", title_facet, "_Facets.png")
                ggsave(filename = file_name_facet)
            }
            
            # Male analysis
            male <- filter(df_subpop_clean, Sex == "Male")
            title_male <- paste(gene, "expression in", sup, "(Male) monocytes (LinReg, rE)", sep = "")
            male_clean <- male[!male$mean_Z %in% boxplot.stats(male$mean_Z)$out, ]
            
            baseline_male <- filter(male_clean, Timepoint == "TP0")
            median_TP0_male <- median(baseline_male$mean_Z)
            
            # ANOVA for Male
            anova_m <- aov(mean_Z ~ Timepoint, data = male_clean)
            p_value_m <- summary(anova_m)[[1]]$"Pr(>F)"[1]
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            timepoints <- c("TP1", "TP2", "TP3", "TP4")
            p_values_tp_m <- sapply(timepoints, function(tp) {
                if (sum(male_clean$Timepoint == tp) > 0) {
                    tryCatch({
                        wilcox.test(male_clean$mean_Z[male_clean$Timepoint == "TP0"],
                                    male_clean$mean_Z[male_clean$Timepoint == tp])$p.value
                    }, error = function(e) NA)
                } else {
                    NA
                }
            })
            
            results <- rbind(results, data.frame(Sex = "Male", Gene = gene, Subpopulation = sup, 
                                                 P_Value = p_value_m, 
                                                 TP1_P_Value = p_values_tp_m[1], 
                                                 TP2_P_Value = p_values_tp_m[2], 
                                                 TP3_P_Value = p_values_tp_m[3], 
                                                 TP4_P_Value = p_values_tp_m[4]))
            
            if (plot_save == "y") {
                ggboxplot(male_clean, 
                          x = "Timepoint", 
                          y = "mean_Z", 
                          color = "Timepoint", 
                          title = title_male,
                          add = "jitter", 
                          legend = "none", 
                          ylab = paste(gene, "expression relative to B2M"), 
                          width = 0.6, 
                          add.params = list(size = 1, alpha = 0.5)) +  
                    geom_hline(yintercept = median_TP0_male, linetype = 2) + 
                    stat_compare_means(method = "anova", label.y = max(male_clean$mean_Z)) +        
                    stat_compare_means(label = "p.signif", 
                                       method = "wilcox", 
                                       ref.group = "TP0", 
                                       hide.ns = TRUE, 
                                       label.y = max(male_clean$mean_Z)) +
                    scale_color_manual(values = my_colors)
                
                file_name_male <- paste0(folder, "/", title_male, ".png")
                ggsave(filename = file_name_male)
            }
            
            # Female analysis
            female <- filter(df_subpop_clean, Sex == "Female")
            title_female <- paste(gene, "expression in", sup, "(Female) monocytes (LinReg, rE)", sep = "")
            female_clean <- female[!female$mean_Z %in% boxplot.stats(female$mean_Z)$out, ]
            
            baseline_female <- filter(female_clean, Timepoint == "TP0")
            median_TP0_female <- median(baseline_female$mean_Z)
            
            # ANOVA for Female
            anova_f <- aov(mean_Z ~ Timepoint, data = female_clean)
            p_value_f <- summary(anova_f)[[1]]$"Pr(>F)"[1]
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            p_values_tp_f <- sapply(timepoints, function(tp) {
                if (sum(female_clean$Timepoint == tp) > 0) {
                    tryCatch({
                        wilcox.test(female_clean$mean_Z[female_clean$Timepoint == "TP0"],
                                    female_clean$mean_Z[female_clean$Timepoint == tp])$p.value
                    }, error = function(e) NA)
                } else {
                    NA
                }
            })
            
            results <- rbind(results, data.frame(Sex = "Female", Gene = gene, Subpopulation = sup, 
                                                 P_Value = p_value_f, 
                                                 TP1_P_Value = p_values_tp_f[1], 
                                                 TP2_P_Value = p_values_tp_f[2], 
                                                 TP3_P_Value = p_values_tp_f[3], 
                                                 TP4_P_Value = p_values_tp_f[4]))
            
            if (plot_save == "y") {
                ggboxplot(female_clean, 
                          x = "Timepoint", 
                          y = "mean_Z", 
                          color = "Timepoint", 
                          title = title_female,
                          add = "jitter", 
                          legend = "none", 
                          ylab = paste(gene, "expression relative to B2M"), 
                          width = 0.6, 
                          add.params = list(size = 1, alpha = 0.5)) +  
                    geom_hline(yintercept = median_TP0_female, linetype = 2) + 
                    stat_compare_means(method = "anova", label.y = max(female_clean$mean_Z)) +        
                    stat_compare_means(label = "p.signif", 
                                       method = "wilcox", 
                                       ref.group = "TP0", 
                                       hide.ns = TRUE, 
                                       label.y = max(female_clean$mean_Z)) +
                    scale_color_manual(values = my_colors)
                
                file_name_female <- paste0(folder, "/", title_female, ".png")
                ggsave(filename = file_name_female)
            }
        }, error = function(e) {
            cat("Error processing", gene, ":", e$message, "\n")
        })
    }
}
results$Significance <- sapply(results$P_Value, get_significance)
results$TP1_Significance <- sapply(results$TP1_P_Value, get_significance)
results$TP2_Significance <- sapply(results$TP2_P_Value, get_significance)
results$TP3_Significance <- sapply(results$TP3_P_Value, get_significance)
results$TP4_Significance <- sapply(results$TP4_P_Value, get_significance)
ANOVA_Sex_Zscore <- results 
# Save ANOVA results with significance levels
write.csv(ANOVA_Sex_Zscore, file = "ANOVA_Zscore_Sex.csv", row.names = FALSE)
