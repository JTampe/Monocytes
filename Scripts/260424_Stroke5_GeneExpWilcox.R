# *Gene Expr. Statistics* ----------------------------------------------------------------------------------------------------
# As they are non-normal disributed, dependent with similar variance
# Adjust to yes if you want plots to be saved or n if not
plot_save <- "n"

# remove unmatched CTR for ANOVA too 
data_mean_matched <- dataFINALmean %>% filter(!(Timepoint == "TP0" & !SampleID %in% Matched_TP0_Gene))
# Check if all groups have the same size:
check_sample_counts(data_mean_matched)

TP_colors <- c("TP0" = "darkgreen", "TP1" = "orange", "TP2" = "red", "TP3" = "magenta", "TP4" = "purple")

#### Wilcox test for Z of Genes --------------------------------------------------
folder <- "Wilcox_capZ_Plots_unpaired"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric(),
                      TP1_P_Value = numeric(), TP2_P_Value = numeric(), TP3_P_Value = numeric(), TP4_P_Value = numeric())
check_sample_counts(data_mean_matched)

j = 1
for (j in 1:length(unique(data_mean_matched$Gene))) {
    gene <- data_mean_matched$Gene[j]
    df_gene <- filter(data_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Category), !is.na(mean_capZ), !is.na(Subpopulation), !is.na(Timepoint))
    i = 1
    for (i in 1:length(unique(df_gene$Subpopulation))) {
        sup <- df_gene$Subpopulation[i]
        df <- filter(df_gene, Subpopulation == sup)
        if (nrow(df) <= 0 || gene %in% exclude_genes) {
            cat("Skipping:", gene, "for subpopulation:", sup, "due to insufficient data or in excluded genes list\n")
            next  # Skip to the next iteration of the loop
        }
        
        # you can only exclud outliers if you are not using paired
        df <- df[!df$mean_capZ %in% boxplot.stats(df$mean_capZ)$out, ]
        df$Timepoint <- factor(as.vector(df$Timepoint))
        baseline <- filter(df, Timepoint == "TP0")
        median_TP0 <- median(baseline$mean_capZ)
        
        # Perform Wilcoxon test for each comparison with TP0
        timepoints <- c("TP1", "TP2", "TP3", "TP4")
        tp_p_values <- c(NA, NA, NA, NA)  # Initialize p-values for TP1 to TP4
        
        # Loop through each timepoint and compare to TP0
        for (tp in seq_along(timepoints)) {
            comparison_tp <- timepoints[tp]
            
            # Check if there are enough data for both TP0 and the current timepoint
            if (nrow(filter(df, Timepoint == comparison_tp)) > 1 && nrow(baseline) > 1) {
                # Perform Wilcoxon test between TP0 and the current timepoint
                test_result <- wilcox.test(df$mean_capZ[df$Timepoint == "TP0"], 
                                           df$mean_capZ[df$Timepoint == comparison_tp],
                                           paired = FALSE, exact = TRUE)
                tp_p_values[tp] <- test_result$p.value  # Store the p-value
            }
        }
        
        # Save results to data frame
        p_value_min <- min(tp_p_values, na.rm = TRUE)
        results <- rbind(results, data.frame(Gene = gene, Subpopulation = sup, 
                                             P_Value = min(tp_p_values, na.rm = TRUE),  # Store the smallest p-value as general
                                             TP1_P_Value = tp_p_values[1], 
                                             TP2_P_Value = tp_p_values[2], 
                                             TP3_P_Value = tp_p_values[3], 
                                             TP4_P_Value = tp_p_values[4]))
        # Plot only if there is a significant p-value (<= 0.05)
        # plot_save <- "n"
        # if (p_value_min <= 0.05) {
        #     plot_save <- "y"
        # }
        # 
        if (plot_save == "y") {
            #my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
            plot_name <- paste(gene, " expression of ", sup, " Monocytes (Control vs post-stroke", sep = "")
            
            # Create a boxplot with significance levels
            p <- ggboxplot(df, 
                           x = "Timepoint", 
                           y = "mean_capZ", 
                           color = "Timepoint", 
                           add = "jitter", 
                           legend = "none", 
                           ylab = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                           title = plot_name,
                           width = 0.8, 
                           add.params = list(size = 1, alpha = 0.5)) +  
                geom_hline(yintercept = median_TP0, linetype = 2) +
                scale_color_manual(values = TP_colors) +
                scale_x_discrete(labels = c(
                    "TP0" = "Control",
                    "TP1" = "24 hours",
                    "TP2" = "3-5 days",
                    "TP3" = "1 month",
                    "TP4" = "3 months"))
            
            # Add significance levels based on Wilcoxon tests
            for (tp in seq_along(timepoints)) {
                comparison_tp <- timepoints[tp]
                p_value <- tp_p_values[tp]
                
                if (!is.na(p_value)) {
                    # Use 'p' for raw p-value or 'p.signif' for significance symbols
                    p <- p + stat_compare_means(method = "wilcox", 
                                                ref.group = "TP0", 
                                                hide.ns = TRUE, 
                                                label = "p.signif",  # You can use 'p.signif' if you prefer symbols
                                                label.y = max(df$mean_capZ))  # Adjust label.y as needed
                }
            }
            
            # Save the plot
            file_name <- file.path(folder, paste(gene, "_", sup, "_capZ_Wilcox.png", sep = ""))
            ggsave(filename = file_name, plot = p)
        }
        
    }
}
results$Significance <- sapply(results$P_Value, get_significance)
results$TP1_Significance <- sapply(results$TP1_P_Value, get_significance)
results$TP2_Significance <- sapply(results$TP2_P_Value, get_significance)
results$TP3_Significance <- sapply(results$TP3_P_Value, get_significance)
results$TP4_Significance <- sapply(results$TP4_P_Value, get_significance)

Wilcox_capZ <- results
write.csv(Wilcox_capZ, file = "Wilcox_Results_capZ_unpaired.csv", row.names = FALSE)

#### Wilcox test for Z of Genes -SubtypesCombined --------------------------------------------------
folder <- "Wilcox_capZ_Plots_SubtypesCombined_unpaired"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric(),
                      TP1_P_Value = numeric(), TP2_P_Value = numeric(), TP3_P_Value = numeric(), TP4_P_Value = numeric())
check_sample_counts(data_mean_matched)

SP_colors <- c("all" = "black", "classical" = "red", "intermediate" = "purple", "nonclassical" = "blue")

# Define new x-axis labels
timepoint_labels <- c(
    "TP0" = "Healthy",
    "TP1" = "24 hours post-stroke",
    "TP2" = "3-5 days post-stroke",
    "TP3" = "1 month post-stroke",
    "TP4" = "3 months post-stroke"
)

for (j in 1:length(unique(data_mean_matched$Gene))) {
    gene <- data_mean_matched$Gene[j]
    df_gene <- filter(data_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Category), !is.na(mean_capZ), !is.na(Subpopulation), !is.na(Timepoint))
    
    if (nrow(df_gene) <= 0 || gene %in% exclude_genes) {
        cat("Skipping:", gene, "due to insufficient data or in excluded genes list\n")
        next
    }
    
    # Exclude outliers
    df_gene <- df_gene[!df_gene$mean_capZ %in% boxplot.stats(df_gene$mean_capZ)$out, ]
    df_gene$Timepoint <- factor(as.vector(df_gene$Timepoint))
    df_gene$Subpopulation <- factor(df_gene$Subpopulation, levels = names(SP_colors))
    baseline <- filter(df_gene, Timepoint == "TP0")
    median_TP0 <- median(baseline$mean_capZ)
    
    # Perform Wilcoxon test for each comparison with TP0, for each Subpopulation
    timepoints <- c("TP1", "TP2", "TP3", "TP4")
    tp_p_values <- data.frame(Subpopulation = unique(df_gene$Subpopulation),
                              TP1 = NA, TP2 = NA, TP3 = NA, TP4 = NA)
    
    for (subpop in unique(df_gene$Subpopulation)) {
        for (tp in seq_along(timepoints)) {
            comparison_tp <- timepoints[tp]
            df_subpop <- filter(df_gene, Subpopulation == subpop)
            baseline_subpop <- filter(df_subpop, Timepoint == "TP0")
            
            if (nrow(filter(df_subpop, Timepoint == comparison_tp)) > 1 && nrow(baseline_subpop) > 1) {
                test_result <- wilcox.test(df_subpop$mean_capZ[df_subpop$Timepoint == "TP0"], 
                                           df_subpop$mean_capZ[df_subpop$Timepoint == comparison_tp],
                                           paired = FALSE, exact = TRUE)
                tp_p_values[tp_p_values$Subpopulation == subpop, timepoints[tp]] <- test_result$p.value
            }
        }
    }
    
    # Save results to data frame
    results <- rbind(results, data.frame(Gene = gene, tp_p_values))
    
    # Calculate the y-axis limit dynamically
    y_max <- max(df_gene$mean_capZ, na.rm = TRUE)
    y_extend <- 0.5 + 0.2 * length(unique(df_gene$Subpopulation))  # Additional space for annotations
    y_axis_limit <- y_max + y_extend
    
    # Create a boxplot for all subpopulations
    p <- ggboxplot(df_gene, 
                   x = "Timepoint", 
                   y = "mean_capZ", 
                   color = "Subpopulation", 
                   add = "jitter", 
                   legend.title = "Subpopulation", 
                   ylab = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                   width = 0.5, 
                   add.params = list(size = 1, alpha = 0.5)) +
        scale_color_manual(values = SP_colors) +
        scale_x_discrete(labels = timepoint_labels) +
        ylim(-3, 3) +  # Set the y-axis limit dynamically
        labs(title = paste("Gene:", gene), 
             subtitle = "Expression across timepoints and subpopulations")
    
    # Add asterisks for significant comparisons
    for (tp in seq_along(timepoints)) {
        comparison_tp <- timepoints[tp]
        y_base <- 2.6
        y_step <- 0.1  # Vertical step between asterisks
        
        for (subpop_idx in seq_along(unique(df_gene$Subpopulation))) {
            subpop <- unique(df_gene$Subpopulation)[subpop_idx]
            p_value <- tp_p_values[tp_p_values$Subpopulation == subpop, comparison_tp]
            
            if (!is.na(p_value) && p_value <= 0.05) {
                asterisk <- ifelse(p_value <= 0.001, "***",
                                   ifelse(p_value <= 0.01, "**",
                                          ifelse(p_value <= 0.05, "*", "")))
                
                if (asterisk != "") {
                    p <- p + annotate("text", 
                                      x = which(levels(df_gene$Timepoint) == comparison_tp), 
                                      y = y_base + y_step * (subpop_idx - 0.1),  # Increment y for each subpopulation
                                      label = asterisk, 
                                      color = SP_colors[subpop], 
                                      size = 5)  # Adjust size as needed
                }
            }
        }
    }
    # Add geom_hline for each subpopulation
    # Add geom_hline for each subpopulation using median TP0 values
    for (subpop in unique(df_gene$Subpopulation)) {
        subpop_data <- filter(df_gene, Subpopulation == subpop & Timepoint == "TP0")
        median_TP0 <- median(subpop_data$mean_capZ, na.rm = TRUE)  # Calculate the median TP0 value for the subpopulation
        
        if (!is.na(median_TP0)) {
            p <- p + geom_hline(yintercept = median_TP0, 
                                linetype = 2, 
                                color = SP_colors[subpop], 
                                size = 0.5)  # Adjust line thickness as needed
        }
    }
    
    
    
    # Save the plot
    file_name <- file.path(folder, paste(gene, "_capZ_Wilcoxon_Subpopulations.png", sep = ""))
    ggsave(filename = file_name, plot = p)
}

# results$Significance <- sapply(results$P_Value, get_significance)
# results$TP1_Significance <- sapply(results$TP1_P_Value, get_significance)
# results$TP2_Significance <- sapply(results$TP2_P_Value, get_significance)
# results$TP3_Significance <- sapply(results$TP3_P_Value, get_significance)
# results$TP4_Significance <- sapply(results$TP4_P_Value, get_significance)
# 
# Wilcox_capZ_SP <- results
# write.csv(Wilcox_capZ, file = "Wilcox_Results_capZ_unpaired.csv", row.names = FALSE)


# *Wilcox for Comparison at individual TPs* -------------------------------------
plot_save <- "n"
#### Male vs Female Analysis (Mean ± SEM) -------------------------------------
folder <- "Mean_SEM_capZ_Sex_Plots"
createFolder(folder)  # Ensure the folder exists
results <- data.frame(Sex = character(), Gene = character(), Subpopulation = character(), 
                      P_Value = numeric(), TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                      TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (i in 1:length(unique(data_mean_matched$Gene))) {
    gene <- unique(data_mean_matched$Gene)[i]
    df_gene <- filter(data_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_capZ), !is.na(Subpopulation), !is.na(Sex))
    
    for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- unique(df_gene$Subpopulation)[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
        
        # Remove outliers
        df_subpop_clean <- df_subpop[!df_subpop$mean_capZ %in% boxplot.stats(df_subpop$mean_capZ)$out, ]
        
        if (nrow(df_subpop_clean) <= 0) next  # Skip if no rows remain
        
        tryCatch({
            title_facet <- paste(gene, "expression in", sup, "monocytes (Zscored dCT)", sep=" ")
            
            # Line plot comparing Male vs Female with mean ± SEM across Timepoints
            ggline(subset(df_subpop_clean, !is.na(Sex)), 
                   x = "Timepoint", 
                   y = "mean_capZ", 
                   color = "Sex", 
                   add = "mean_se",  # Add mean and standard error
                   size = 1.2) +  
                labs(y = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                     title = title_facet) +
                scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2")) +
                theme_bw() +  # Apply white background
                theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
                      panel.grid.minor = element_blank(), 
                      panel.border = element_blank(),
                      axis.line = element_line(color = "black")) +
                stat_compare_means(method = "wilcox.test",  # or "t.test" 
                                   aes(group = Sex),
                                   label = "p.signif",  
                                   size = 5) 
            
            if (plot_save == "y") {
                file_name_facet <- paste0(folder, "/", title_facet, "_Mean_SEM.png")
                ggsave(filename = file_name_facet)
            }
            
            # Initialize a list to store p-values for timepoints
            p_values_list <- list(TP0 = NA, TP1 = NA, TP2 = NA, TP3 = NA, TP4 = NA)
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            timepoints <- c("TP0", "TP1", "TP2", "TP3", "TP4")
            for (tp in timepoints) {
                male_values <- df_subpop_clean$mean_capZ[df_subpop_clean$Timepoint == tp & df_subpop_clean$Sex == "Male"]
                female_values <- df_subpop_clean$mean_capZ[df_subpop_clean$Timepoint == tp & df_subpop_clean$Sex == "Female"]
                
                # Perform Wilcoxon test for the comparison of Male vs Female at each timepoint
                if (length(male_values) > 0 && length(female_values) > 0) {
                    wilcox_result <- wilcox.test(male_values, female_values)
                    p_values_list[[tp]] <- wilcox_result$p.value
                }
            }
            
            # Add results for Male vs Female for all timepoints as a single row
            results <- rbind(results, data.frame(WilcoxComparison = "Male vs Female", Gene = gene, Subpopulation = sup, 
                                                 TP0_P_Value = p_values_list[["TP0"]],
                                                 TP1_P_Value = p_values_list[["TP1"]],
                                                 TP2_P_Value = p_values_list[["TP2"]],
                                                 TP3_P_Value = p_values_list[["TP3"]],
                                                 TP4_P_Value = p_values_list[["TP4"]]))
            
            
        }, error = function(e) {
            cat("Error processing", gene, ":", e$message, "\n")
        })
    }
}
Timecourse_Wilcox_Sex_Z <- results
Timecourse_Wilcox_Sex_Z$TP0_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP0_P_Value, get_significance)
Timecourse_Wilcox_Sex_Z$TP1_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP1_P_Value, get_significance)
Timecourse_Wilcox_Sex_Z$TP2_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP2_P_Value, get_significance)
Timecourse_Wilcox_Sex_Z$TP3_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP3_P_Value, get_significance)
Timecourse_Wilcox_Sex_Z$TP4_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP4_P_Value, get_significance)
write.csv(Timecourse_Wilcox_Sex_Z , file = "Timecourse_Wilcox_Sex_Z_Results.csv", row.names = FALSE)

#### Healthy vs MILD vs MODERATE  (Mean ± SEM)* -------------------------------------
folder <- "Mean_SEM_capZ_Category_Plots"
createFolder(folder)  # Ensure the folder exists
results <- data.frame(Category = character(), Gene = character(), Subpopulation = character(), 
                      P_Value = numeric(), TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                      TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (i in 1:length(unique(data_mean_matched$Gene))) {
    gene <- unique(data_mean_matched$Gene)[i]
    df_gene <- filter(data_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_capZ), !is.na(Subpopulation), !is.na(Category))
    
    for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- unique(df_gene$Subpopulation)[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
        
        # Remove outliers
        df_subpop_clean <- df_subpop[!df_subpop$mean_capZ %in% boxplot.stats(df_subpop$mean_capZ)$out, ]
        
        if (nrow(df_subpop_clean) <= 0) next  # Skip if no rows remain
        
        tryCatch({
            title_facet <- paste(gene, "expression in", sup, "monocytes (Zscored dCT)", sep=" ")
            
            # Line plot comparing Male vs Female with mean ± SEM across Timepoints
            ggline(subset(df_subpop_clean, !is.na(Category)), 
                   x = "Timepoint", 
                   y = "mean_capZ", 
                   color = "Category", 
                   add = "mean_se",  # Add mean and standard error
                   size = 1.2) +  
                labs(y = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                     title = title_facet) +
                scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606")) +
                theme_bw() +  # Apply white background
                theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
                      panel.grid.minor = element_blank(), 
                      panel.border = element_blank(),
                      axis.line = element_line(color = "black")) +
                stat_compare_means(method = "wilcox.test",  # or "t.test" 
                                   aes(group = Category),
                                   label = "p.signif",  
                                   size = 5) 
            
            if (plot_save == "y") {
                file_name_facet <- paste0(folder, "/", title_facet, "_Mean_SEM.png")
                ggsave(filename = file_name_facet)
            }
            
            # Initialize a list to store p-values for timepoints
            p_values_list <- list(TP0 = NA, TP1 = NA, TP2 = NA, TP3 = NA, TP4 = NA)
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            timepoints <- c("TP0", "TP1", "TP2", "TP3", "TP4")
            for (tp in timepoints) {
                minor_values <- df_subpop_clean$mean_capZ[df_subpop_clean$Timepoint == tp & df_subpop_clean$Category == "MINOR"]
                moderate_values <- df_subpop_clean$mean_capZ[df_subpop_clean$Timepoint == tp & df_subpop_clean$Category == "MODERATE"]
                
                # Perform Wilcoxon test for the comparison of Minor vs Moderate at each timepoint
                if (length(minor_values) > 0 && length(moderate_values) > 0) {
                    wilcox_result <- wilcox.test(minor_values, moderate_values)
                    p_values_list[[tp]] <- wilcox_result$p.value
                }
            }
            
            # Add results for Male vs Female for all timepoints as a single row
            results <- rbind(results, data.frame(WilcoxComparison = "Minor vs Moderate", Gene = gene, Subpopulation = sup, 
                                                 TP0_P_Value = p_values_list[["TP0"]],
                                                 TP1_P_Value = p_values_list[["TP1"]],
                                                 TP2_P_Value = p_values_list[["TP2"]],
                                                 TP3_P_Value = p_values_list[["TP3"]],
                                                 TP4_P_Value = p_values_list[["TP4"]]))
            
        }, error = function(e) {
            cat("Error processing", gene, ":", e$message, "\n")
        })
    }
}

Timecourse_Wilcox_Category_Z <- results
Timecourse_Wilcox_Category_Z$TP0_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP0_P_Value, get_significance)
Timecourse_Wilcox_Category_Z$TP1_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP1_P_Value, get_significance)
Timecourse_Wilcox_Category_Z$TP2_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP2_P_Value, get_significance)
Timecourse_Wilcox_Category_Z$TP3_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP3_P_Value, get_significance)
Timecourse_Wilcox_Category_Z$TP4_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP4_P_Value, get_significance)
write.csv(Timecourse_Wilcox_Category_Z , file = "Timecourse_Wilcox_Category_Z_Results.csv", row.names = FALSE)

#### Good vs Bad Recovery (Mean ± SEM)* -------------------------------------
folder <- "Mean_SEM_capZ_Recovery_Plots"
plot_save <- "n"
createFolder(folder)  # Ensure the folder exists
results <- data.frame(Recovery = character(), Gene = character(), Subpopulation = character(), 
                      P_Value = numeric(), TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                      TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (i in 1:length(unique(data_mean_matched$Gene))) {
    gene <- unique(data_mean_matched$Gene)[i]
    df_gene <- filter(data_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_capZ), !is.na(Subpopulation), !is.na(Recovery))
    
    for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- unique(df_gene$Subpopulation)[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
        
        # Remove outliers
        df_subpop_clean <- df_subpop[!df_subpop$mean_capZ %in% boxplot.stats(df_subpop$mean_capZ)$out, ]
        
        if (nrow(df_subpop_clean) <= 0) next  # Skip if no rows remain
        
        tryCatch({
            title_facet <- paste(gene, "expression in", sup, "monocytes (Zscored dCT)", sep=" ")
            
            # Line plot comparing Male vs Female with mean ± SEM across Timepoints
            ggline(subset(df_subpop_clean, !is.na(Recovery)), 
                   x = "Timepoint", 
                   y = "mean_capZ", 
                   color = "Recovery", 
                   add = "mean_se",  # Add mean and standard error
                   size = 1.2) +  
                labs(y = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                     title = title_facet) +
                scale_color_manual(values = c("Good" = "darkgreen", "Bad" = "#700606")) +
                theme_bw() +  # Apply white background
                theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
                      panel.grid.minor = element_blank(), 
                      panel.border = element_blank(),
                      axis.line = element_line(color = "black")) +
                stat_compare_means(method = "wilcox.test",  # or "t.test" 
                                   aes(group = Recovery),
                                   label = "p.signif",  
                                   size = 5) 
            
            if (plot_save == "y") {
                file_name_facet <- paste0(folder, "/", title_facet, "_Mean_SEM.png")
                ggsave(filename = file_name_facet, width = 9, height = 5)
            }
            
            # Initialize a list to store p-values for timepoints
            p_values_list <- list(TP0 = NA, TP1 = NA, TP2 = NA, TP3 = NA, TP4 = NA)
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            timepoints <- c("TP0", "TP1", "TP2", "TP3", "TP4")
            for (tp in timepoints) {
                good_values <- df_subpop_clean$mean_capZ[df_subpop_clean$Timepoint == tp & df_subpop_clean$Recovery == "Good"]
                bad_values <- df_subpop_clean$mean_capZ[df_subpop_clean$Timepoint == tp & df_subpop_clean$Recovery == "Bad"]
                
                # Perform Wilcoxon test for the comparison of Good vs Bad Recovery at each timepoint
                if (length(good_values) > 0 && length(bad_values) > 0) {
                    wilcox_result <- wilcox.test(good_values, bad_values)
                    p_values_list[[tp]] <- wilcox_result$p.value
                }
            }
            
            # Add results for Male vs Female for all timepoints as a single row
            results <- rbind(results, data.frame(WilcoxComparison = "Good vs Bad Recovery", Gene = gene, Subpopulation = sup, 
                                                 TP0_P_Value = p_values_list[["TP0"]],
                                                 TP1_P_Value = p_values_list[["TP1"]],
                                                 TP2_P_Value = p_values_list[["TP2"]],
                                                 TP3_P_Value = p_values_list[["TP3"]],
                                                 TP4_P_Value = p_values_list[["TP4"]]))
            
        }, error = function(e) {
            cat("Error processing", gene, ":", e$message, "\n")
        })
    }
}

Timecourse_Wilcox_Recovery_Z <- results
Timecourse_Wilcox_Recovery_Z$TP0_Significance <- sapply(Timecourse_Wilcox_Recovery_Z$TP0_P_Value, get_significance)
Timecourse_Wilcox_Recovery_Z$TP1_Significance <- sapply(Timecourse_Wilcox_Recovery_Z$TP1_P_Value, get_significance)
Timecourse_Wilcox_Recovery_Z$TP2_Significance <- sapply(Timecourse_Wilcox_Recovery_Z$TP2_P_Value, get_significance)
Timecourse_Wilcox_Recovery_Z$TP3_Significance <- sapply(Timecourse_Wilcox_Recovery_Z$TP3_P_Value, get_significance)
Timecourse_Wilcox_Recovery_Z$TP4_Significance <- sapply(Timecourse_Wilcox_Recovery_Z$TP4_P_Value, get_significance)
write.csv(Timecourse_Wilcox_Recovery_Z , file = "Timecourse_Wilcox_Recovery_Z_Results.csv", row.names = FALSE)

#### All results -----------------
TP_Wilcox_All_Z <- rbind(Timecourse_Wilcox_Recovery_Z, Timecourse_Wilcox_Category_Z,Timecourse_Wilcox_Sex_Z)
write.csv(TP_Wilcox_All_Z , file = "TP_comparisons_Wilcox_All_Z_Results.csv", row.names = FALSE)

TP_Wilcox_Allsig_Z <- TP_Wilcox_All_Z %>% filter(if_any(4:8, ~ . <= 0.05))
write.csv(TP_Wilcox_Allsig_Z , file = "TP_comparisons_Wilcox_Allsig_Z_Results.csv", row.names = FALSE)
