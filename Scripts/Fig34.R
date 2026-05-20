# FACS ANOVA
automate_anova_extraction <- function(results_folder, plots_folder, results_name, plot_save, dataset, loop_vars, timepoint_col = "Timepoint") {
    # Create directories if they don't exist
    if (!dir.exists(results_folder)) {dir.create(results_folder)}
    if (!dir.exists(plots_folder)) {dir.create(plots_folder)}
    
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
            cat("Not enough data to run Wilcoxon test for:", var_name, "\n")
            next
        }
        # Initialize timepoint p-values as NA
        tp1_p_value <- tp2_p_value <- tp3_p_value <- tp4_p_value <- NA
        
        # Perform Wilcoxon test comparing each timepoint to TP0
        try({
            if ("TP1" %in% df_test[[timepoint_col]]) {
                tp1_p_value <- wilcox.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP1"], 
                                           df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], na.rm = TRUE)$p.value
            }
            if ("TP2" %in% df_test[[timepoint_col]]) {
                tp2_p_value <- wilcox.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP2"], 
                                           df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], na.rm = TRUE)$p.value
            }
            if ("TP3" %in% df_test[[timepoint_col]]) {
                tp3_p_value <- wilcox.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP3"], 
                                           df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], na.rm = TRUE)$p.value
            }
            if ("TP4" %in% df_test[[timepoint_col]]) {
                tp4_p_value <- wilcox.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP4"], 
                                           df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], na.rm = TRUE)$p.value
            }
        }, silent = TRUE)
        
        # Save the Wilcoxon p-values to the results dataframe
        results <- rbind(results, data.frame(Variable = var_name, 
                                             P_Value = min(tp1_p_value, tp2_p_value, tp3_p_value, tp4_p_value, na.rm = TRUE),
                                             TP1_P_Value = tp1_p_value, TP2_P_Value = tp2_p_value,
                                             TP3_P_Value = tp3_p_value, TP4_P_Value = tp4_p_value))
        
        # Generate and save the plots if required
        if (plot_save == "y") {
            median_TP0 <- median(df_test[df_test[[timepoint_col]] == "TP0", var_name], na.rm = TRUE)
            my_colors <- c("darkgrey", "black", "black", "black", "black")
            df_test$Group <- ifelse(df_test[[timepoint_col]] == "TP0", "Controls", "Stroke Patients")
            
            my_fill_colors <- c(
                "TP0" = scales::alpha("darkgrey", 0.25),
                "TP1" = scales::alpha("black", 0.5),
                "TP2" = scales::alpha("black", 0.5),
                "TP3" = scales::alpha("black", 0.5),
                "TP4" = scales::alpha("black", 0.5)
            )
            
            box_plot <- ggplot(df_test,
                               aes(x = .data[[timepoint_col]],
                                   y = .data[[var_name]],
                                   color = Group)) +
                
                geom_boxplot(
                    aes(fill = .data[[timepoint_col]]),
                    width = 0.65,
                    outlier.shape = NA,
                    show.legend = TRUE
                ) +
                
                geom_jitter(
                    aes(fill = .data[[timepoint_col]]),
                    shape = 21,
                    size = 2,
                    alpha = 0.6,
                    width = 0.2,
                    show.legend = FALSE
                ) +
                
                geom_hline(yintercept = median_TP0, linetype = 2, color ="darkgrey") +
                
                scale_color_manual(values = c(
                    "Controls" = "darkgrey",
                    "Stroke Patients" = "black"
                )) +
                
                scale_fill_manual(values = my_fill_colors, guide = "none") +
                
                scale_x_discrete(labels = c(
                    "TP0" = "Control",
                    "TP1" = "24 hours",
                    "TP2" = "3-5 days",
                    "TP3" = "1 month",
                    "TP4" = "3 months"
                )) +
                
                #scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
                
                ylab(if (var_name == "Monocytes") {
                    "Percentage of PBMCs [%]"
                } else {
                    "Percentage of monocytes [%]"
                }) +
                
                xlab("Time after stroke") +
                
                theme_classic(base_size = 14) +
                theme(
                    axis.title.x = element_text(hjust = 1),
                    legend.position = "bottom",
                    legend.justification = c(0, 0.5),
                    legend.box.just = "left",
                    legend.margin = margin(t = -20),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12),
                    plot.margin = margin(t = 30, r = 5, b = 5, l = 5),
                    legend.title = element_blank()
                )
            
            tp_levels <- levels(df_test[[timepoint_col]])
            y_pos <- max(df_test[[var_name]], na.rm = TRUE) * 1.05
            
            p_vals <- c(tp1_p_value, tp2_p_value, tp3_p_value, tp4_p_value)
            names(p_vals) <- c("TP1","TP2","TP3","TP4")
            
            for (tp in names(p_vals)) {
                p <- p_vals[tp]
                
                if (!is.na(p) && p <= 0.05) {
                    symbol <- ifelse(p <= 0.001, "***",
                                     ifelse(p <= 0.01, "**", "*"))
                    
                    x_pos <- which(tp_levels == tp)
                    
                    box_plot <- box_plot +
                        annotate("text",
                                 x = x_pos,
                                 y = y_pos,
                                 label = symbol,
                                 size = 14,
                                 color = "black")
                }
            }
            
            
            # Save boxplot
            boxplot_file_name <- file.path(plots_folder, paste(var_name, "_Boxplot_Wilcox.png", sep = ""))
            ggsave(filename = boxplot_file_name, plot = box_plot, width = 7.5, height = 6.5)
            
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
                                  ylab = paste(var_name, "mean (+-SEM) percentage"), 
                                  add = "mean_se", width = 0.8) +
                scale_fill_manual(values = my_colors) +
                scale_x_discrete(labels = c(
                    "TP0" = "Control",
                    "TP1" = "24 hours",
                    "TP2" = "3-5 days",
                    "TP3" = "1 month",
                    "TP4" = "3 months"
                )) +
                geom_errorbar(aes(ymin = mean_value - SEM, ymax = mean_value + SEM), width = 0.2) +
                stat_compare_means(data = df_test, aes(x = get(timepoint_col), y = get(var_name)),
                                   method = "wilcox", ref.group = "TP0", hide.ns = TRUE,
                                   label = "p.signif", label.y = max(mean_values$mean_value) + 0.09 * max(mean_values$mean_value),
                                   size = 8)
            
            # Save barplot
            barplot_file_name <- file.path(plots_folder, paste(var_name, "_Barplot_Wilcox_SEM.png", sep = ""))
            ggsave(filename = barplot_file_name, plot = bar_plot)
        }
    }
    
    # Write the results to a CSV file
    results_file <- file.path(paste(results_name, "_Results.csv", sep = ""))
    write.csv(results, file = results_file, row.names = FALSE)
    
    # Return the results dataframe
    return(results)
}

automate_anova_extraction_Category <- function(results_folder, plots_folder, results_name, plot_save, dataset, loop_vars, timepoint_col = "Timepoint", color_var) {
    # Create directories if they don't exist
    if (!dir.exists(results_folder)) dir.create(results_folder)
    if (!dir.exists(plots_folder)) dir.create(plots_folder)
    
    # Initialize an empty dataframe to store results
    results <- data.frame(Category = character(), Variable = character(), P_Value = numeric(),
                          TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                          TP3_P_Value = numeric(), TP4_P_Value = numeric())
    
    # Loop through the specified variables (columns) in the dataset
    for (var_name in loop_vars) {
        # Subset the dataset for the variable and timepoints
        df_test <- dataset[, c(var_name, timepoint_col, color_var)]
        df_test <- df_test %>% filter(!is.na(get(var_name)), !is.na(get(timepoint_col)), !is.na(get(color_var)))
        
        # Ensure Timepoint is treated as a factor
        df_test[[timepoint_col]] <- factor(df_test[[timepoint_col]], levels = c("TP0", "TP1", "TP2", "TP3", "TP4"))
        
        # Check if we have sufficient data
        if (length(unique(df_test[[timepoint_col]])) < 2 || nrow(df_test) < 3) next
        
        for (color_vari in unique(df_test[[color_var]])) {
            df_color <- df_test %>% filter(get(color_var) == color_vari)
            tp_p_values <- sapply(c("TP1", "TP2", "TP3", "TP4"), function(tp) {
                if (tp %in% df_color[[timepoint_col]]) {
                    tryCatch(
                        wilcox.test(df_color[[var_name]][df_color[[timepoint_col]] == tp], 
                                    df_color[[var_name]][df_color[[timepoint_col]] == "TP0"])$p.value,
                        error = function(e) NA
                    )
                } else {
                    NA
                }
            })
            
            results <- rbind(results, data.frame(Category = color_vari, Variable = var_name, 
                                                 P_Value = min(tp_p_values, na.rm = TRUE),
                                                 TP1_P_Value = tp_p_values["TP1"], 
                                                 TP2_P_Value = tp_p_values["TP2"],
                                                 TP3_P_Value = tp_p_values["TP3"], 
                                                 TP4_P_Value = tp_p_values["TP4"]))
        }
        
        # Generate plots if required
        if (plot_save == "y") {
            median_TP0 <- median(df_test[df_test[[timepoint_col]] == "TP0", var_name], na.rm = TRUE)
            my_colors <- c("Female" = "#1D04C2","Male" = "#A67C00", 
                           "MINOR" = "#A67C00","MODERATE" = "#700606",
                           "Bad" = "#700606","Good" = "#4CAF50")
            
            # Filter my_colors to match the levels of color_var in df_test
            used_colors <- my_colors[unique(df_test[[color_var]])]
            
            # Boxplot with Correct Legend
            box_plot <- ggplot(df_test, aes(x = .data[[timepoint_col]],
                                            y = .data[[var_name]],
                                            color = .data[[color_var]])) +
                
                # Boxplots with different fill alpha per group
                geom_boxplot(
                    aes(fill = .data[[color_var]],
                        alpha = .data[[color_var]]),
                    position = position_dodge(width = 0.75),
                    width = 0.65,
                    outlier.shape = NA
                ) +
                
                # Jitter points (unchanged look)
                geom_jitter(
                    aes(fill = .data[[color_var]]),
                    shape = 21,
                    size = 2,
                    alpha = 0.6,
                    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
                ) +
                
                # Your colors
                scale_color_manual(values = my_colors) +
                scale_fill_manual(values = my_colors, guide = "none") +
                
                # Alpha ONLY for boxplots
                scale_alpha_manual(values = c("Female" = 0.5, "Male" = 0.15, "Bad" = 0.5, "Good" = 0.15),
                                   guide = "none") +
                
                scale_x_discrete(labels = c(
                    "TP0" = "Control",
                    "TP1" = "24 hours",
                    "TP2" = "3-5 days",
                    "TP3" = "1 month",
                    "TP4" = "3 months"
                )) +
                
                ylab( if (var_name == "Monocytes") {"Percentage of PBMCs [%]"} 
                    else {"Percentage of monocytes [%]"})+
                xlab("Time after stroke") +
            
                theme_classic(base_size = 14) +
                theme(axis.title.x = element_text(hjust = 1),
                      legend.position = "bottom",
                      legend.justification = c(0, 0.5),
                      legend.box.just = "left",
                      legend.margin = margin(t = -20),
                      axis.text.x = element_text(size = 12),  
                      axis.text.y = element_text(size = 12),
                      plot.margin = margin(t = 30, r = 5, b = 5, l = 5),
                      legend.title = element_blank()
                      )
            
            # Annotate significant comparisons for each category
            for (color in unique(df_test[[color_var]])) {
                df_color <- results[results$Category == color & results$Variable == var_name, ]
                y_base <- max(df_test[[var_name]], na.rm = TRUE) + 5
                for (tp in c("TP1", "TP2", "TP3", "TP4")) {
                    p_value <- df_color[[paste(tp, "P_Value", sep = "_")]]
                    if (!is.na(p_value) && p_value <= 0.05) {
                        asterisk <- ifelse(p_value <= 0.001, "***",
                                           ifelse(p_value <= 0.01, "**", "*"))
                        x_position <- which(levels(df_test[[timepoint_col]]) == tp)
                        box_plot <- box_plot +
                            annotate("text", x = x_position, y = y_base, 
                                     label = asterisk, size = 12, color = my_colors[color])
                    }
                }
            }
            ggsave(filename = file.path(plots_folder, paste(var_name, "_Boxplot_Wilcox.png", sep = "")), 
                   plot = box_plot, width = 7.5, height = 6.5)
        }
    }
    # Save results to a CSV file
    write.csv(results, file.path(results_folder, paste(results_name, "_results.csv", sep = "")), row.names = FALSE)
    
    # Return the results dataframe
    return(results)
}

# *FACS Statistics* ----------------
FACSdata_copy <- FACSdata 
#### 1WAY ANOVA* ----------------------------------------------------------------------------------------------------

ANOVA_FACSdata <- FACSdata %>% filter(!(Timepoint == "TP0" & SampleID %in% Unmatched_TP0_FACS))
write_xlsx(ANOVA_FACSdata, "ANOVA_FACSdata.xlsx")

# with wilcox - unpaired and normality not fully given
ANOVA_FACSall <- automate_anova_extraction(output_location, "Wilcox_FACS_Plots", 
                                           "Wilcox_FACS", "y", ANOVA_FACSdata, 
                                           colnames(ANOVA_FACSdata)[9:12], 
                                           "Timepoint")

ANOVA_FACSdata <- ANOVA_FACSdata[order(ANOVA_FACSdata$Sex == "Female", decreasing = TRUE), ]
ANOVA_FACSsex <- automate_anova_extraction_Category(output_location, "Wilcox_FACS_Plots_Sex",
                                                    "Wilcox_FACS_Sex", "y", ANOVA_FACSdata,
                                                    colnames(ANOVA_FACSdata)[9:12],
                                                    c(colnames(ANOVA_FACSdata)[2]),c(colnames(ANOVA_FACSdata)[22]))

ANOVA_FACSdata <- ANOVA_FACSdata[order(ANOVA_FACSdata$Recovery == "Bad", decreasing = TRUE), ]
ANOVA_FACSrecovery <- automate_anova_extraction_Category(output_location, "Wilcox_FACS_Plots_Recovery",
                                                         "Wilcox_FACS_Recovery", "y", ANOVA_FACSdata,
                                                         colnames(ANOVA_FACSdata)[9:12],
                                                         c(colnames(ANOVA_FACSdata)[2]),c(colnames(ANOVA_FACSdata)[25]))