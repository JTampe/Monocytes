automate_ttest_extraction <- function(results_folder, plots_folder, results_name, plot_save, dataset, loop_vars, timepoint_col = "Timepoint") {
    # Create directories if they don't exist
    if (!dir.exists(results_folder)) {
        dir.create(results_folder)
    }
    
    if (!dir.exists(plots_folder)) {
        dir.create(plots_folder)
    }
    
    # Initialize an empty dataframe to store results
    results <- data.frame(Variable = character(), P_Value = numeric(),
                          TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                          TP3_P_Value = numeric(), TP4_P_Value = numeric())
    
    # Loop through the specified variables (columns) in the dataset
    for (var_name in loop_vars) {
        # Subset the dataset for the variable and timepoints
        df_test <- dataset[, c(var_name, timepoint_col)]
        
        # Filter out rows with missing values
        df_test <- df_test %>% filter(!is.na(get(var_name)), !is.na(get(timepoint_col)))
        
        # Ensure Timepoint is treated as a factor
        df_test[[timepoint_col]] <- factor(df_test[[timepoint_col]], levels = c("TP0", "TP1", "TP2", "TP3", "TP4"))
        
        # Check if we have sufficient data for all timepoints
        if (length(unique(df_test[[timepoint_col]])) < 2) {
            cat("Not enough timepoints with data for:", var_name, "\n")
            next
        }
        
        # Skip if not enough data
        if (nrow(df_test) < 3) {
            cat("Not enough data to run paired t-test for:", var_name, "\n")
            next
        }
        
        # Initialize timepoint p-values as NA
        tp1_p_value <- tp2_p_value <- tp3_p_value <- tp4_p_value <- NA
        
        # Perform paired t-test comparing each timepoint to TP0
        try({
            if ("TP1" %in% df_test[[timepoint_col]]) {
                tp1_p_value <- t.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP1"], 
                                      df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], 
                                      paired = TRUE, na.rm = TRUE)$p.value
            }
            if ("TP2" %in% df_test[[timepoint_col]]) {
                tp2_p_value <- t.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP2"], 
                                      df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], 
                                      paired = TRUE, na.rm = TRUE)$p.value
            }
            if ("TP3" %in% df_test[[timepoint_col]]) {
                tp3_p_value <- t.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP3"], 
                                      df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], 
                                      paired = TRUE, na.rm = TRUE)$p.value
            }
            if ("TP4" %in% df_test[[timepoint_col]]) {
                tp4_p_value <- t.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP4"], 
                                      df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], 
                                      paired = TRUE, na.rm = TRUE)$p.value
            }
        }, silent = TRUE)
        
        # Save the paired t-test p-values to the results dataframe
        results <- rbind(results, data.frame(Variable = var_name, 
                                             P_Value = min(tp1_p_value, tp2_p_value, tp3_p_value, tp4_p_value, na.rm = TRUE),
                                             TP1_P_Value = tp1_p_value, TP2_P_Value = tp2_p_value,
                                             TP3_P_Value = tp3_p_value, TP4_P_Value = tp4_p_value))
        
        # Generate and save the plots if required
        if (plot_save == "y") {
            median_TP0 <- median(df_test[df_test[[timepoint_col]] == "TP0", var_name], na.rm = TRUE)
            my_colors <- c("darkgreen", "orange", "red", "magenta", "purple")
            
            # Boxplot generation
            box_plot <- ggboxplot(df_test, x = timepoint_col, y = var_name, color = timepoint_col, 
                                  add = "jitter", legend = "none", 
                                  ylab = paste(var_name, "percentage"), width = 0.8, 
                                  add.params = list(size = 1, alpha = 0.5)) +
                geom_hline(yintercept = median_TP0, linetype = 2) +
                stat_compare_means(label = "p.signif", 
                                   method = "t.test", 
                                   ref.group = "TP0", 
                                   paired = TRUE, 
                                   hide.ns = TRUE, 
                                   label.y = max(df_test[[var_name]])) +
                scale_color_manual(values = my_colors)
            
            # Save boxplot
            boxplot_file_name <- file.path(plots_folder, paste(var_name, "_Boxplot_Paired_TTest.png", sep = ""))
            ggsave(filename = boxplot_file_name, plot = box_plot)
            
            # Barplot generation with SEM
            mean_values <- df_test %>%
                group_by(get(timepoint_col)) %>%
                summarise(mean_value = mean(get(var_name), na.rm = TRUE),
                          sd_value = sd(get(var_name), na.rm = TRUE),
                          n = n()) %>%
                mutate(SEM = sd_value / sqrt(n)) %>%
                rename(Timepoint = `get(timepoint_col)`)
            
            # Barplot with SEM and significance
            bar_plot <- ggbarplot(mean_values, x = "Timepoint", y = "mean_value", fill = "Timepoint",
                                  ylab = paste(var_name, "mean value"), 
                                  add = "mean_se", width = 0.8) +
                scale_fill_manual(values = my_colors) +
                geom_errorbar(aes(ymin = mean_value - SEM, ymax = mean_value + SEM), width = 0.2) +
                stat_compare_means(data = df_test, aes(x = get(timepoint_col), y = get(var_name)),
                                   method = "t.test", ref.group = "TP0", paired = TRUE, 
                                   hide.ns = TRUE, label = "p.signif", 
                                   label.y = max(mean_values$mean_value) + 0.1 * max(mean_values$mean_value))
            
            # Save barplot
            barplot_file_name <- file.path(plots_folder, paste(var_name, "_Barplot_Paired_TTest_SEM.png", sep = ""))
            ggsave(filename = barplot_file_name, plot = bar_plot)
            
        }
    }
    
    # Write the results to a CSV file
    results_file <- file.path(paste(results_name, "_Results.csv", sep = ""))
    write.csv(results, file = results_file, row.names = FALSE)
    
    # Return the results dataframe
    return(results)
}
