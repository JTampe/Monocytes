## Set Up of Environment
# adjust it to all the other packages
if (requireNamespace("thematic")) 
  thematic::thematic_rmd(font = "auto")

today <- Sys.Date()

output_location <- paste(today,"_Stroke_Results", sep="")

setwd("/Users/ju5263ta/Github/Monocytes/Data/")
if (!dir.exists(output_location)) {
    dir.create(output_location, recursive = TRUE)
  }
getwd ()

rm(list=ls())

library(ggplot2)
library(tidyverse)
library(stats)
library(skimr)
library(sjPlot)
library(readxl)
library(thematic)
library(knitr)
library(lme4)
library(ggpubr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(stringr)
library(colorRamp2)
library(rstatix)
library(scales)
#library(emmeans)
library(pwr)
library(lmtest)
library(writexl)
library(openxlsx)
#library(nondetects)
library(skimr)
library(visdat) 
library(mice)
library(lmerTest)
library(beepr)
#library(Seurat)

thematic_on(bg = "white", fg = "black", accent = "blue")

### *Functions* ---------------------------------------------------------------------
my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
my_grey_scale <- c("grey", "black", "black", "black", "black") 

## Create Folder
createFolder <- function(folderName) {
  if (!dir.exists(folderName)) {
    dir.create(folderName, recursive = TRUE)
  }
}

## Split data
createSeparatedDataFiles <- function(data, categoryColumn, categoryFolder) {
  uniqueCategories <- unique(data[[categoryColumn]])
  
  for (category in uniqueCategories) {
    categoryData <- data %>% filter(data[[categoryColumn]] == category)
    categoryName <- as.character(category)
    filename <- file.path(categoryFolder, paste0(categoryName, ".txt"))
    write.table(categoryData, filename, sep = "\t", row.names = FALSE)
    assign(categoryName, categoryData, envir = .GlobalEnv)
  }
  
  return(uniqueCategories)
}
## Function to get the failes samples
getMeTheFailedSamples<- function(data, gene, cutoff) {
    # Extract the Sample IDs where the gene's value exceeds the cutoff
    unique(data$Sample[which(data$Gene == gene & (data$CT_Value > cutoff | is.na(data$CT_Value)))])
}
## Function to perform correlation and linear regression by Timepoint
calculate_linreg_by_timepoint <- function(data, cells_colname) {
    
    # Initialize an empty list to store results for each Timepoint
    results_list <- list()
    
    # Loop through each unique Timepoint and perform the correlation and regression
    for (timepoint in unique(data$Timepoint)) {
        
        # Subset data for the current Timepoint
        data_subset <- subset(data, Timepoint == timepoint)
        
        # Pearson correlation
        cor_test <- cor.test(data_subset[[cells_colname]], data_subset$Age)
        
        # Linear regression model
        lm_result <- lm(data_subset[[cells_colname]] ~ data_subset$Age)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        intercept <- summary(lm_result)$coefficients[1, 1]
        
        # Calculate Fold Induction
        Max <- max(data_subset$Age, na.rm = TRUE)
        Min <- min(data_subset$Age, na.rm = TRUE)
        FinalValue <- coef_x * Max + intercept
        InitialValue <- coef_x * Min + intercept
        Fold <- FinalValue / InitialValue
        
        # Store the results in a data frame for this Timepoint
        result <- data.frame(Cells = cells_colname,
                             Timepoint = timepoint,
                             N = length(data_subset[[cells_colname]]),
                             p.value = cor_test$p.value,
                             Coefficient = coef_x,
                             Down95 = ci[1],
                             Up95 = ci[2],
                             Intercept = intercept,
                             Fold = Fold,
                             Max_age = Max,
                             Min_age = Min)
        
        # Append to the results list
        results_list[[timepoint]] <- result
    }
    
    # Combine all results into a single data frame
    final_results <- do.call(rbind, results_list)
    
    return(final_results)
}
# Function to plot linear regression for a specific group
plot_linear_regression <- function(data, cells_colname, group_name, output_folder) {
    
    # Fit the initial linear model to identify outliers
    lm_initial <- lm(data[[cells_colname]] ~ data$Age)
    residuals_values <- residuals(lm_initial)
    outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
    
    # Filter out the outliers
    data_filtered <- data[-outlier_indices, ]
    
    # Check if there are enough data points to fit the model
    if (nrow(data_filtered) < 2) {
        warning(paste("Not enough data points after outlier removal for", group_name, "(", cells_colname, "). Skipping plot."))
        return(NULL)
    }
    
    # Calculate linear regression results
    rows_group <- calculate_linreg_by_timepoint(data_filtered, cells_colname)
    if (is.null(rows_group) || nrow(rows_group) == 0) {
        warning(paste("No results for", group_name, "(", cells_colname, "). Skipping plot."))
        return(NULL)
    }
    
    rows_group$Focus <- group_name  # Add the group label
    
    Age_correlation_FACS <<- rbind(Age_correlation_FACS, rows_group)  # Update global variable
    
    # Plot regression for the specific group with all timepoints
    p_group <- ggplot(data_filtered, aes(x = Age, y = .data[[cells_colname]], color = Timepoint)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", aes(color = Timepoint), se = FALSE, size = 1) +
        ggtitle(paste("Linear Regression -", group_name, "(", cells_colname, ")", sep = " ")) +
        xlab("Age") +
        ylab(cells_colname) +
        theme(text = element_text(size = 14)) +
        theme(
            panel.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "white", color = "black"),
            text = element_text(color = "black"),
            panel.grid.major = element_line(color = "gray", size = 0.2),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()
        ) +
        scale_color_manual(values = my_colors) +
        scale_fill_manual(values = my_colors)
    
    # Add confidence interval only if p-value is <= 0.05
    for (timepoint in unique(data_filtered$Timepoint)) {
        if (rows_group$p.value[rows_group$Timepoint == timepoint] <= 0.05) {
            p_group <- p_group + geom_smooth(data = subset(data_filtered, Timepoint == timepoint),
                                             method = "lm",
                                             aes(fill = Timepoint),
                                             se = TRUE,
                                             alpha = 0.25,
                                             color = NA)
        }
    }
    
    # Save the plot with all timepoints
    ggsave(filename = file.path(output_folder, paste0("LinReg_", group_name, "_", cells_colname, ".png")), plot = p_group)
    
    # Create an additional plot with only significant timepoints if there are any
    significant_timepoints <- rows_group$Timepoint[rows_group$p.value <= 0.05]
    
    if (length(significant_timepoints) > 0) {
        # Plot with all significant timepoints together
        p_significant <- ggplot(data_filtered %>% filter(Timepoint %in% significant_timepoints), 
                                aes(x = Age, y = .data[[cells_colname]], color = Timepoint)) +
            geom_point(size = 2) +
            geom_smooth(method = "lm", aes(color = Timepoint), se = TRUE, size = 1) +
            ggtitle(paste("Linear Regression - Significant Only -", group_name, "(", cells_colname, ")", sep = " ")) +
            xlab("Age") +
            ylab(cells_colname) +
            theme(text = element_text(size = 14)) +
            theme(
                panel.background = element_rect(fill = "white", color = "black"),
                plot.background = element_rect(fill = "white", color = "black"),
                text = element_text(color = "black"),
                panel.grid.major = element_line(color = "gray", size = 0.2),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_manual(values = my_colors) +
            scale_fill_manual(values = my_colors)
        
        # Save the plot with only significant timepoints together
        ggsave(filename = file.path(output_folder, paste0("LinReg_SignificantOnly_", group_name, "_", cells_colname, ".png")), plot = p_significant)
        
        # Plot individual plots for each significant timepoint
        for (sig_timepoint in significant_timepoints) {
            p_individual <- ggplot(data_filtered %>% filter(Timepoint == sig_timepoint), 
                                   aes(x = Age, y = .data[[cells_colname]], color = Timepoint)) +
                geom_point(size = 2) +
                geom_smooth(method = "lm", aes(color = Timepoint), se = TRUE, size = 1) +
                ggtitle(paste("Linear Regression -", group_name, "(", cells_colname, ")", " - Timepoint:", sig_timepoint)) +
                xlab("Age") +
                ylab(cells_colname) +
                theme(text = element_text(size = 14)) +
                theme(
                    panel.background = element_rect(fill = "white", color = "black"),
                    plot.background = element_rect(fill = "white", color = "black"),
                    text = element_text(color = "black"),
                    panel.grid.major = element_line(color = "gray", size = 0.2),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.x = element_blank()
                ) +
                scale_color_manual(values = my_colors[sig_timepoint]) +
                scale_fill_manual(values = my_colors[sig_timepoint])
            
            # Save the individual plot
            ggsave(filename = file.path(output_folder, paste0("LinReg_SignificantOnly_", group_name, "_", cells_colname, "_Timepoint_", sig_timepoint, ".png")), plot = p_individual)
        }
    }
}

plot_linear_regression_NHISS <- function(data, cells_colname, group_name, output_folder) {
    
    # Fit the initial linear model to identify outliers
    lm_initial <- lm(data[[cells_colname]] ~ data$NHISS)
    residuals_values <- residuals(lm_initial)
    outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
    
    # Filter out the outliers
    data_filtered <- data[-outlier_indices, ]
    
    # Check if there are enough data points to fit the model
    if (nrow(data_filtered) < 2) {
        warning(paste("Not enough data points after outlier removal for", group_name, "(", cells_colname, "). Skipping plot."))
        return(NULL)
    }
    
    # Calculate linear regression results on filtered data
    rows_group <- calculate_linreg_by_timepoint(data_filtered, cells_colname)
    if (is.null(rows_group) || nrow(rows_group) == 0) {
        warning(paste("No results for", group_name, "(", cells_colname, "). Skipping plot."))
        return(NULL)
    }
    
    rows_group$Focus <- group_name  # Add the group label
    
    NHISS_correlation_FACS <<- rbind(NHISS_correlation_FACS, rows_group)  # Update global variable
    
    # Plot regression for the specific group with all timepoints
    p_group <- ggplot(data_filtered, aes(x = NHISS, y = .data[[cells_colname]], color = Timepoint)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", aes(color = Timepoint), se = FALSE, size = 1) +
        ggtitle(paste("NHISS Linear Regression -", group_name, "(", cells_colname, ")", sep = " ")) +
        xlab("NHISS") +
        ylab(cells_colname) +
        theme(text = element_text(size = 14)) +
        theme(
            panel.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "white", color = "black"),
            text = element_text(color = "black"),
            panel.grid.major = element_line(color = "gray", size = 0.2),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()
        ) +
        scale_color_manual(values = my_colors) +
        scale_fill_manual(values = my_colors)
    
    # Add confidence interval only if p-value is <= 0.05
    for (timepoint in unique(data_filtered$Timepoint)) {
        if (rows_group$p.value[rows_group$Timepoint == timepoint] <= 0.05) {
            p_group <- p_group + geom_smooth(data = subset(data_filtered, Timepoint == timepoint),
                                             method = "lm",
                                             aes(fill = Timepoint),
                                             se = TRUE,
                                             alpha = 0.25,
                                             color = NA)
        }
    }
    
    # Save the plot with all timepoints
    ggsave(filename = file.path(output_folder, paste0("NHISS_LinReg_", group_name, "_", cells_colname, ".png")), plot = p_group)
    
    # Create an additional plot with only significant timepoints if there are any
    significant_timepoints <- rows_group$Timepoint[rows_group$p.value <= 0.05]
    
    if (length(significant_timepoints) > 0) {
        # Plot with all significant timepoints together
        p_significant <- ggplot(data_filtered %>% filter(Timepoint %in% significant_timepoints), 
                                aes(x = NHISS, y = .data[[cells_colname]], color = Timepoint)) +
            geom_point(size = 2) +
            geom_smooth(method = "lm", aes(color = Timepoint), se = TRUE, size = 1) +
            ggtitle(paste("NHISS Linear Regression - Significant Only -", group_name, "(", cells_colname, ")", sep = " ")) +
            xlab("NHISS") +
            ylab(cells_colname) +
            theme(text = element_text(size = 14)) +
            theme(
                panel.background = element_rect(fill = "white", color = "black"),
                plot.background = element_rect(fill = "white", color = "black"),
                text = element_text(color = "black"),
                panel.grid.major = element_line(color = "gray", size = 0.2),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_manual(values = my_colors) +
            scale_fill_manual(values = my_colors)
        
        # Save the plot with only significant timepoints together
        ggsave(filename = file.path(output_folder, paste0("NHISS_LinReg_SignificantOnly_", group_name, "_", cells_colname, ".png")), plot = p_significant)
        
        # Plot individual plots for each significant timepoint
        for (sig_timepoint in significant_timepoints) {
            p_individual <- ggplot(data_filtered %>% filter(Timepoint == sig_timepoint), 
                                   aes(x = NHISS, y = .data[[cells_colname]], color = Timepoint)) +
                geom_point(size = 2) +
                geom_smooth(method = "lm", aes(color = Timepoint), se = TRUE, size = 1) +
                ggtitle(paste("NHISS Linear Regression -", group_name, "(", cells_colname, ")", " - Timepoint:", sig_timepoint)) +
                xlab("NHISS") +
                ylab(cells_colname) +
                theme(text = element_text(size = 14)) +
                theme(
                    panel.background = element_rect(fill = "white", color = "black"),
                    plot.background = element_rect(fill = "white", color = "black"),
                    text = element_text(color = "black"),
                    panel.grid.major = element_line(color = "gray", size = 0.2),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.x = element_blank()
                ) +
                scale_color_manual(values = my_colors[sig_timepoint]) +
                scale_fill_manual(values = my_colors[sig_timepoint])
            
            # Save the individual plot
            ggsave(filename = file.path(output_folder, paste0("NHISS_LinReg_SignificantOnly_", group_name, "_", cells_colname, "_Timepoint_", sig_timepoint, ".png")), plot = p_individual)
        }
    }
}

# summary data.frame for dotplot
calculate_summary <- function(data, combo_col, gene_col) {
    # Calculate mean_capZ
    mean_summary <- data %>%
        group_by(!!sym(combo_col), !!sym(gene_col)) %>%
        summarize(
            mean_capZ = mean(mean_capZ, na.rm = TRUE),  
            count = n()  
        ) %>%
        filter(count > 1) 
    
    # Calculate sd_mean_capZ
    sd_summary <- data %>%
        group_by(!!sym(combo_col), !!sym(gene_col)) %>%
        summarize(
            sd_mean_capZ = sd(mean_capZ, na.rm = TRUE)  # Calculate standard deviation of mean_capZ
        ) 
    
    # Combine both summaries
    final_summary <- mean_summary %>%
        left_join(sd_summary, by = c(combo_col, gene_col)) %>%
        mutate(
            Dot_Size = scales::rescale(sd_mean_capZ, to = c(5, 1))  # Rescale SD to range 5 (low SD) to 1 (high SD)
        )
    
    return(final_summary)
}

check_sample_counts <- function(data, timepoint_col = "Timepoint", sampleid_col = "SampleID") {
    # List unique SampleID values for each timepoint
    unique_samples_tp0 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP0"])
    unique_samples_tp1 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP1"])
    unique_samples_tp2 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP2"])
    unique_samples_tp3 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP3"])
    unique_samples_tp4 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP4"])
    
    # Store lengths of unique samples for each timepoint
    lengths <- c(
        TP0 = length(unique_samples_tp0),
        TP1 = length(unique_samples_tp1),
        TP2 = length(unique_samples_tp2),
        TP3 = length(unique_samples_tp3),
        TP4 = length(unique_samples_tp4)
    )
    
    # Check if all timepoints have the same number of unique SampleIDs
    all_equal_lengths <- all(lengths == lengths[1])
    
    # Print result of the comparison
    if (all_equal_lengths) {
        cat("All timepoints have the same number of unique SampleIDs.\n")
    } else {
        cat("Not all timepoints have the same number of unique SampleIDs.\n")
    }
    
    # Optionally return the lengths for further inspection
    return(lengths)
}

# get significance stars
get_significance <- function(p_value) {
    ifelse(is.na(p_value), "NA",  # Handle NA values
           ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                         ifelse(p_value < 0.05, "*", "ns"))))
}

# extrTP p-values
get_p_value <- function(tp1, tp2, tukey_comparisons, tukey_result) {
    comparison <- paste(tp1, tp2, sep = "-")
    reverse_comparison <- paste(tp2, tp1, sep = "-")
    
    if (comparison %in% tukey_comparisons) {
        return(tukey_result$Timepoint[comparison, "p adj"])
    } else if (reverse_comparison %in% tukey_comparisons) {
        return(tukey_result$Timepoint[reverse_comparison, "p adj"])
    } else {
        return(NA)
    }
}

my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
# FACS ANOVA
automate_anova_extraction <- function(results_folder, plots_folder, results_name, plot_save, dataset, loop_vars, timepoint_col = "Timepoint") {
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
            my_colors <- c("darkgreen", "orange", "red", "magenta", "purple")
            
            # Boxplot generation
            box_plot <- ggboxplot(df_test, x = timepoint_col, y = var_name, color = timepoint_col, 
                                  add = "jitter", legend = "none", 
                                  ylab = paste(var_name, "percentage"), width = 0.8, 
                                  add.params = list(size = 1, alpha = 0.5)) +
                geom_hline(yintercept = median_TP0, linetype = 2) +
                stat_compare_means(label = "p.signif", 
                                   method = "wilcox", 
                                   ref.group = "TP0", 
                                   hide.ns = TRUE, 
                                   label.y = max(df_test[[var_name]])) +
                scale_color_manual(values = my_colors)
            
            # Save boxplot
            boxplot_file_name <- file.path(plots_folder, paste(var_name, "_Boxplot_Wilcox.png", sep = ""))
            ggsave(filename = boxplot_file_name, plot = box_plot)
            
            # Barplot generation with SEM
            mean_values <- df_test %>%
                group_by(get(timepoint_col)) %>%
                summarise(mean_value = mean(get(var_name), na.rm = TRUE),
                          sd_value = sd(get(var_name), na.rm = TRUE),
                          n = n()) %>%
                mutate(SEM = sd_value / sqrt(n)) %>%
                rename(Timepoint = `get(timepoint_col)`)
            
            #my_grey_scale <- c("grey", "black", "black", "black", "black") 
            
            # Barplot with SEM and significance
            bar_plot <- ggbarplot(mean_values, x = "Timepoint", y = "mean_value", fill = "Timepoint",
                                  ylab = paste(var_name, "mean value"), 
                                  add = "mean_se", width = 0.8) +
                scale_fill_manual(values = my_colors) +
                geom_errorbar(aes(ymin = mean_value - SEM, ymax = mean_value + SEM), width = 0.2) +
                stat_compare_means(data = df_test, aes(x = get(timepoint_col), y = get(var_name)),
                                   method = "wilcox", ref.group = "TP0", hide.ns = TRUE,
                                   label = "p.signif", label.y = max(mean_values$mean_value) + 0.1 * max(mean_values$mean_value))
            
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

demographics_N_Age_Sex <- function(df) {
    # Initialize an empty dataframe to store the results
    Demographics_FACS <- data.frame(Variable = character(), Controls = character(), Patients = character(), P_value = numeric())
    
    # Filter for controls and patients based on SampleID
    df_Controls <- df %>% filter(SampleID %in% c(unique(df$SampleID[grep("Ctr", df$SampleID)])))
    df_Patients <- df %>% filter(SampleID %in% c(unique(df$SampleID[grep("Patient", df$SampleID)])))
    
    # Calculate N (sample count)
    N_control <- length(unique(df_Controls$SampleID))
    N_patients <- length(unique(df_Patients$SampleID))
    row_N <- data.frame(Variable = "N", 
                        Controls = N_control, 
                        Patients = N_patients, 
                        P_value = NA)  
    
    # Calculate Age, mean (SD) for controls and patients
    mean_age_control <- round(mean(df_Controls$Age, na.rm = TRUE), 2)
    sd_age_control <- round(sd(df_Controls$Age, na.rm = TRUE), 2)
    mean_age_patients <- round(mean(df_Patients$Age, na.rm = TRUE), 2)
    sd_age_patients <- round(sd(df_Patients$Age, na.rm = TRUE), 2)
    
    # Perform t-test for Age comparison between controls and patients
    t_test_age <- t.test(df_Controls$Age, df_Patients$Age)
    p_value_age <- round(t_test_age$p.value, 3)
    
    row_Age <- data.frame(Variable = "Age", 
                          Controls = paste(mean_age_control, " ± ", sd_age_control, sep = ""), 
                          Patients = paste(mean_age_patients, " ± ", sd_age_patients, sep = ""),
                          P_value = p_value_age)  
    
    # Calculate percentage of females in controls and patients
    sex_table_controls <- table(df_Controls$Sex)
    sex_table_patients <- table(df_Patients$Sex)
    percent_female_controls <- round((sex_table_controls["F"] / sum(sex_table_controls)) * 100, 2)
    percent_female_patients <- round((sex_table_patients["F"] / sum(sex_table_patients)) * 100, 2)
    
    # Perform chi-square test for sex distribution comparison
    sex_table_combined <- table(df_Controls$Sex, df_Patients$Sex)
    chi_test_sex <- chisq.test(sex_table_combined)
    p_value_sex <- round(chi_test_sex$p.value, 3)
    
    row_Sex <- data.frame(Variable = "Female %", 
                          Controls = percent_female_controls, 
                          Patients = percent_female_patients, 
                          P_value = p_value_sex)
    
    # Combine rows into the Demographics table
    Demographics_FACS <- rbind(Demographics_FACS, row_N, row_Age, row_Sex)
    
    return(Demographics_FACS)
}

demographics_Compare_N_Age_Sex <- function(df, split_column) {
    # Initialize an empty dataframe to store the results
    Demographics_FACS <- data.frame(Variable = character(), Controls = character(), Patients = character(), P_value = numeric())
    
    # Split the data based on the specified column
    groups <- unique(df[[split_column]])
    if (length(groups) != 2) {
        stop("The specified column must contain exactly two groups (e.g., 'Control' and 'Patient').")
    }
    
    # Extract controls and patients based on the specified column
    df_Controls <- df %>% filter(!!sym(split_column) == groups[1])
    df_Patients <- df %>% filter(!!sym(split_column) == groups[2])
    
    # Calculate N (sample count)
    N_control <- length(unique(df_Controls$SampleID))
    N_patients <- length(unique(df_Patients$SampleID))
    row_N <- data.frame(Variable = "N", 
                        Controls = N_control, 
                        Patients = N_patients, 
                        P_value = NA)  
    
    # Calculate Age, mean (SD) for controls and patients
    mean_age_control <- round(mean(df_Controls$Age, na.rm = TRUE), 2)
    sd_age_control <- round(sd(df_Controls$Age, na.rm = TRUE), 2)
    mean_age_patients <- round(mean(df_Patients$Age, na.rm = TRUE), 2)
    sd_age_patients <- round(sd(df_Patients$Age, na.rm = TRUE), 2)
    
    # Perform t-test for Age comparison between controls and patients
    t_test_age <- t.test(df_Controls$Age, df_Patients$Age)
    p_value_age <- round(t_test_age$p.value, 3)
    
    row_Age <- data.frame(Variable = "Age", 
                          Controls = paste(mean_age_control, " ± ", sd_age_control, sep = ""), 
                          Patients = paste(mean_age_patients, " ± ", sd_age_patients, sep = ""),
                          P_value = p_value_age)  
    
    # Calculate percentage of females in controls and patients
    sex_table_controls <- table(df_Controls$Sex)
    sex_table_patients <- table(df_Patients$Sex)
    percent_female_controls <- round((sex_table_controls["F"] / sum(sex_table_controls)) * 100, 2)
    percent_female_patients <- round((sex_table_patients["F"] / sum(sex_table_patients)) * 100, 2)
    
    # Perform chi-square test for sex distribution comparison
    sex_table_combined <- table(df_Controls$Sex, df_Patients$Sex)
    chi_test_sex <- chisq.test(sex_table_combined)
    p_value_sex <- round(chi_test_sex$p.value, 3)
    
    row_Sex <- data.frame(Variable = "Female %", 
                          Controls = percent_female_controls, 
                          Patients = percent_female_patients, 
                          P_value = p_value_sex)
    
    # Combine rows into the Demographics table
    Demographics_FACS <- rbind(Demographics_FACS, row_N, row_Age, row_Sex)
    
    return(Demographics_FACS)
}

plot_my_pca = function(pca_result, data.meta, what, PCA_color=NULL, legend_pos=NULL){
    # Create a scatter plot of PC1 vs PC2
    if (is.null(PCA_color) ){
        PCA_color = rainbow( length(unique(data.meta[,what])))
    }
    
    plot(pca_result$x[, 1], pca_result$x[, 2], 
         xlab = "PC1", ylab = "PC2", 
         main = "PCA Scatter Plot (PC1 vs PC2)", 
         col = c(PCA_color)[as.numeric( factor( data.meta[,what]))], pch = 16)
    # Optional: Add legend for PCA_Sample colors
    if (!is.null(legend_pos) ){
    legend(legend_pos, legend = unique(data.meta[,what]), col = unique(PCA_color), pch = 16, title = "PCA Sample")
} }

## *Data* ---------------------------------------------------------------------

#### Loading of the Data -----------------------------------------
setwd("/Users/ju5263ta/Github/Monocytes/rawData_Stroke")
getwd ()
# make sure that before you check the files, that not an empty column has been added to the end of the file, remove and export as .csv otherwise.
files <- system( "ls *B2M.csv", intern=T)

# if you just have one data set:
i=1
raw_data <- read.csv(text=readLines(files[i])[-(1:11)], header = T, sep=',', dec = '.',
                     row.names = 1
                     )

# delete Call.1 column, and resent column names
if ( length(colnames(raw_data)) > 13) {
    raw_data <- raw_data[, !colnames(raw_data) %in% "Call.1"]
    colnames(raw_data) <- c("Name","Type","rConc","Name.1","Type.1","Reference","Value",
                            "Quality","Call","Threshold","Value.1","Quality.1","Call.1")
}
raw_data$Dataset<-files[i]

# for more the 1 data set:
i = 2
for (i in (2:length(files))) {
    y = raw_data
    x <- read.csv(text=readLines(files[i])[-(1:11)], header = T, sep=',', dec = '.',row.names = 1)
    
    if (length(colnames(x))>13) {
        x <- x[, !colnames(x) %in% "Call.1"]
        colnames(x) <- c("Name","Type","rConc","Name.1","Type.1","Reference","Value",
                         "Quality","Call","Threshold","Value.1","Quality.1","Call.1")
    }
    x$Dataset<-files[i]
    y <- rbind(y, x)
    i = i +1
    raw_data = y
} 

#### trimming of data ------------------------------------------------------
data <- cbind.data.frame(Sample=raw_data$Name,
                         Gene=raw_data$Name.1,
                         Value=raw_data$Value,
                         dCT_Value=raw_data$Value.1,
                         Call=raw_data$Call.1)

#### evaluate PASS FAIL------------------------------------------------
data$Value<-ifelse(data$Call=='Flag', yes=35, no=data$Value) #999 
data$Value<-ifelse(data$Call=='mFlag', yes=NA, no=data$Value) #999 
data$dCT_Value<-ifelse(data$Call=='Flag', yes=35, no=data$dCT_Value) #999
data$dCT_Value<-ifelse(data$Call=='mFlag', yes=NA, no=data$dCT_Value) #999 instead of NA

# sort through quality and give NA to everything that Failed (=999) or is above CT 35
data$Value<-ifelse(data$Value>=35, yes=35, no=data$Value) 
data$dCT_Value<-ifelse(data$dCT_Value>=35, yes=35, no=data$dCT_Value) 
#data$dCT_Value<-ifelse(data$dCT_Value>=35, yes=35, no=data$dCT_Value) 

# Trim your data
dataTRIM <-cbind.data.frame(Sample=data$Sample,
                            Gene1=data$Gene,
                            CT_Value=data$Value,
                            dCT_Value=data$dCT_Value,
                            Comment = NA)

# Information on the data ---------------------------------------------------------------------
# adjust the indicators for the genes list here (for the heatmap) - line 138:

dataTRIM$TechnicalReplicate <- raw_data$Name.1
info_genes = t(data.frame(lapply( dataTRIM$TechnicalReplicate, function(x){ 
    ret = unlist( stringr::str_split( x, "[_\\s\\.]")); 
    if ( length(ret) < 2){
        ret = c(ret, rep(0, 2-length(ret) ))
    } 
    ret  } ) ))
colnames(info_genes) = c("Gene", "TechnicalReplicate")

info_genes <- as.data.frame(info_genes)

dataTRIM$TechnicalReplicate <- NULL
dataTRIM = cbind(dataTRIM, info_genes)

# List of Genes
gene_names <- unique(dataTRIM$Gene)

# List of Genes
genes <- unique(dataTRIM$Gene)
genes

# Number of Genes
no_genes <- length(genes)
no_genes

# List of Samples
samples <- unique(dataTRIM$Sample)

# Split technical and (pseudo-)biological replicates....:

#without the B1/B2 segregation, we would have each read of the gene normalized against the mean of the two biological repeats of B2M. This is in theory correct but we lose the information on the mean/SD while doing so through the machine.
#without the T1/T2 separation, we have the plate as a whole normalized against all 4 reads... with it we can split the whole plate into two halves, each standing for its own technical read.

info_rows <- t(data.frame(
  sapply(dataTRIM$Sample, function(x) {
    ret <- unlist(stringr::str_split(x, "[_\\s\\.]"))
    if (length(ret) < 5) {
      ret <- c(ret, rep(NA, 5 - length(ret)))
    }
    ret[1:5]  # Select only the first five elements
  })
))
colnames(info_rows) = c("Subpopulation", "Cell_Count", "SampleID","Timepoint","BiologicalReplicate")

# NRT: no Polymerase, technical control 
# NTC: no Template, neg-Sample control   
# 10xLC:Linarity control
info_rows <- as.data.frame(info_rows)

# Label Subpopulations:
# remove the C in NC to avoid problems:
info_rows$Subpopulation <- gsub("NC","N",as.character(info_rows$Subpopulation))
info_rows$Subpopulation[grep("N",info_rows$Subpopulation)]  <- "nonclassical"
info_rows$Subpopulation[grep("M",info_rows$Subpopulation)]  <- "all"
info_rows$Subpopulation[grep("C",info_rows$Subpopulation)]  <- "classical"
info_rows$Subpopulation[grep("I",info_rows$Subpopulation)]  <- "intermediate"

dataTRIM = cbind(dataTRIM, info_rows)

## remove unnessessary info (for now)
#It's easier to make the imputation on less columnes/variables

# Remove Extra sorts
dataTRIM <- dataTRIM %>% filter(BiologicalReplicate != "Extra")
dataTRIM <- dataTRIM %>% filter(Subpopulation != "noRT")
dataTRIM <- dataTRIM %>% filter(Cell_Count != "200")
# remove all Timepoint 5
dataTRIM <- dataTRIM %>% filter(Timepoint != "TP5")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test01")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test02")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test03")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test04")

dataTRIM$Cell_Count <- NULL

#### Expand data --------------------------------------------------------
unique_SampleID <- unique(dataTRIM$SampleID[!grepl("^Ctr", dataTRIM$SampleID)])
unique_Gene <- unique(dataTRIM$Gene)
unique_Subpopulation <- unique(dataTRIM$Subpopulation)
unique_Timepoint <- unique(dataTRIM$Timepoint[!grepl("TP0", dataTRIM$Timepoint)])
unique_TechnicalReplicate <- unique(dataTRIM$TechnicalReplicate)
unique_BiologicalReplicate <- unique(dataTRIM$BiologicalReplicate)

expanded_df <- expand.grid(SampleID = unique_SampleID,
                           Gene = unique_Gene,
                           Subpopulation = unique_Subpopulation,
                           Timepoint = unique_Timepoint,
                           TechnicalReplicate = unique_TechnicalReplicate,
                           BiologicalReplicate = unique_BiologicalReplicate)

dataEXPAND <- merge(expanded_df, dataTRIM, by = c("SampleID", "Gene", "Subpopulation", "Timepoint", "TechnicalReplicate", "BiologicalReplicate"), all.x = TRUE)

dataEXPAND$SampleID <- as.character(dataEXPAND$SampleID)
dataEXPAND$Gene <- as.character(dataEXPAND$Gene)
dataEXPAND$Subpopulation <- as.character(dataEXPAND$Subpopulation)
dataEXPAND$Timepoint <- as.character(dataEXPAND$Timepoint)
dataEXPAND$TechnicalReplicate <- as.character(dataEXPAND$TechnicalReplicate)
dataEXPAND$BiologicalReplicate <- as.character(dataEXPAND$BiologicalReplicate)

dataCTR <- dataTRIM[grep("^Ctr", dataTRIM$SampleID), ]

dataEXPAND <- rbind(dataEXPAND, dataCTR)

dataEXPAND$Comment <- ifelse(is.na(dataEXPAND$Sample), yes = "No Fludigm run, imputed Values",  no = dataEXPAND$Comment)

#### Metadata & FACSdata---------------------------------------------------------------
setwd("/Users/ju5263ta/Github/Monocytes/rawData_Stroke")
metadataP <- read_excel("Metadata_PatientID.xlsx")
Metadata_GeneID <- read_excel("Metadata_GeneID.xlsx")
Plate_info <- read_excel("241014_PlateID.xlsx")
Metadata_NHISS <- read_excel("Metadata_NHISS.xlsx")

dataAll <- merge(dataEXPAND, metadataP, by = "SampleID", all.x = TRUE)
#### match the right controls
#Matched_TP0_Gene <- c("Ctr01", "Ctr02", "Ctr08", "Ctr09", "Ctr11", "Ctr14", "Ctr16", "Ctr17", "Ctr18", "Ctr25", "Ctr27", "Ctr29", "Ctr31","Ctr41", "Ctr43")
Matched_TP0_Gene <- c("Ctr20", "Ctr07", "Ctr10", "Ctr14", "Ctr25", "Ctr24", "Ctr18", "Ctr19", 
                      "Ctr08", "Ctr26", "Ctr17", "Ctr09", "Ctr27","Ctr16", "Ctr11")

Unmatched_TP0_FACS <- c("Ctr16", "Ctr21", "Ctr23", "Ctr30", "Ctr39", "Ctr44", "Ctr33")


FACSdata <- as.data.frame(read_excel("241008_StrokeControlMFI_combined.xls", # 241008_StrokeControlMFI_combined.xls for updated one
                                     col_types = c("text", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric")))

# split Sample ID into all info
FACS_info_rows <- t(data.frame(sapply(FACSdata$...1, function(x) {
    ret <- unlist(stringr::str_split(x, "[_\\s\\.]"))
    if (length(ret) < 4) {
        ret <- c(ret, rep(NA, 4 - length(ret)))
    }
    ret[1:4]  # Select only the first five elements
})
))

colnames(FACS_info_rows) = c("Type", "Plate", "Sorting Map","Sorted Sample")
FACS_info_rows <- as.data.frame(FACS_info_rows)
FACS_info_rows$PlateID <- paste(FACS_info_rows$Plate,FACS_info_rows$`Sorted Sample` )
FACSdata <- cbind(FACSdata, FACS_info_rows)

# Match PlateID with the right sample with the right data
FACSdata <- merge(FACSdata, Plate_info, by = "PlateID")
# remove the 5th timepoint
FACSdata <- FACSdata %>% filter(Timepoint != "TP5")

# add the metadata of Patients
FACSdata <- merge(FACSdata, metadataP, by = "SampleID", all.x = TRUE)

# add the metadata of NHISS score of Patients
FACSdata <- merge(FACSdata, Metadata_NHISS, by = c("SampleID", "Timepoint"), all.x = TRUE)

# Combine the data frames based on common values in the "id" column
# --> here i loose all the test because they are not in the metadata
dataAll <- merge(dataEXPAND, metadataP, by = "SampleID", all.x = TRUE)

#relabel Male & Female
dataAll$Sex[grep("M",dataAll$Sex)]  <- "Male"
dataAll$Sex[grep("F",dataAll$Sex)]  <- "Female"

FACSdata$Sex[grep("M",FACSdata$Sex)]  <- "Male"
FACSdata$Sex[grep("F",FACSdata$Sex)]  <- "Female"

dataAll <- dataAll %>% filter(!is.na(Age), !is.na(SampleID), !is.na(Category), !is.na(Sex))

#### failed Genes --------------------------------------------------------
fails <- dataAll %>%
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarize(
        N = n(),
        Count_Above_35_CT = sum(CT_Value >= 35, na.rm = TRUE) + sum(is.na(CT_Value)),  # Include NAs in the count
    )

fails$Per_Above_35_CT = fails$Count_Above_35_CT/fails$N * 100

failed_genes <- fails %>% filter(Per_Above_35_CT >= 90) %>% select(Gene)
wrong_threshold <- c("CLEC7A", "S100A8", "S100A9")
exclude_genes  <- unique(c("Xeno",failed_genes$Gene, wrong_threshold))
paste("The gene",unique(failed_genes$Gene), "was excluded, as more than 90% of the reads failed.")
paste("The gene",wrong_threshold, "was excluded, as the global threshold setting did not allow proper evaluation.")

# remove the failed genes:
dataWORKING <- dataAll %>%  filter(!Gene %in% exclude_genes)

#### failed Samples --------------------------------------------------------
fails_Sample <- dataWORKING %>%
    group_by(Sample) %>%
    summarize(
        N = n(),
        Count_Above_35_CT <- sum(CT_Value >= 35, na.rm = TRUE),
        Per_Above_35_CT = mean(CT_Value >= 35, na.rm = TRUE) * 100,
        Count_Above_35_dCT = sum(dCT_Value >= 35, na.rm = TRUE),
        Per_Above_35_dCT = mean(dCT_Value >= 35, na.rm = TRUE) * 100
    )

failed_samples <- fails_Sample %>% filter(Per_Above_35_CT >= 30) %>% select(Sample)
failed_samples <- failed_samples$Sample
paste("The Sample ",failed_samples, "was excluded, as more than 30% of the reads failed.")

# Set all Values from failed Samples (> 30%) NA:
dataWORKING$CT_Value <- ifelse(dataWORKING$Sample %in% failed_samples, yes = NA, no = dataWORKING$CT_Value)
dataWORKING$dCT_Value <- ifelse(dataWORKING$Sample %in% failed_samples, yes = NA, no = dataWORKING$dCT_Value)
dataWORKING$Comment <- ifelse(dataWORKING$Sample %in% failed_samples, yes = "More 30% of CT values were 35 and above, imputed Values", no = dataWORKING$Comment)
unique(dataWORKING$Comment)

#### Evaluate Housekeeping genes -----------
FailedHK <-table(c( 
    getMeTheFailedSamples(dataWORKING, "ACTB", 30),
    getMeTheFailedSamples(dataWORKING, "B2M", 30),
    getMeTheFailedSamples(dataWORKING, "GAPDH", 30)))

hist(as.numeric(FailedHK), 
     main = "Histogram of Failed Samples", 
     xlab = "Failed Samples Count", 
     breaks = 10)

failedSAMPLES <- names(FailedHK) [ which(FailedHK > 2) ]

# Set CT & dCT NA for the Sampels with failed HK gene
# Set data$Value to NA for the samples that are in the failedSAMPLES vector
dataWORKING$CT_Value <- ifelse(dataWORKING$Sample %in% failedSAMPLES, yes = NA, no = dataWORKING$CT_Value)
dataWORKING$dCT_Value <- ifelse(dataWORKING$Sample %in% failedSAMPLES, yes = NA, no = dataWORKING$dCT_Value)
dataWORKING$Comment <- ifelse(dataWORKING$Sample %in% failedSAMPLES, yes = "More than 2 houskeeing gene reads failed", no = dataWORKING$Comment)

# Filter the data for the desired genes
filtered_data <- dataWORKING %>%
    filter(Gene %in% c("B2M", "ACTB", "GAPDH"))

# Calculate the CV for each group
variance_table <- filtered_data %>%
    group_by(Gene, Timepoint, Subpopulation) %>%
    summarize(CV = sd(CT_Value, na.rm = TRUE) / mean(CT_Value, na.rm = TRUE)) %>%
    ungroup()

# Display the variance table
print(variance_table)

# Plotting the variances
ggplot(variance_table, aes(x = Timepoint, y = CV, fill = Subpopulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Gene, scales = "fixed") +
    labs(title = "Coefficient of Variation Across Timepoints and Subpopulations",
         x = "Timepoint", y = "Coefficient of Variation (CV)") +
    theme_minimal() +
    scale_fill_manual(values = c("darkgreen", "orange", "red", "magenta"))

# *Test: Pool TP1 & TP2 and TP3 & TP4* --------------------------------------------------------------------------------------------
# add TP2 to TP1
# dataWORKING$Timepoint[grep("TP2",dataWORKING$Timepoint)]  <- "TP1"
# add TP2 to TP1
# dataWORKING$Timepoint[grep("TP4",dataWORKING$Timepoint)]  <- "TP3"


# *Imputation* --------------------------------------------------------------------------------------------
#### Prepare Matrix for imputation ---------------------------------------------------------------------
# visualize the missing values:
vis_miss(dataWORKING)
#hist(dataWORKING$dCT_Value, main = "Histogram of Original Data")

initial2<-mice(dataWORKING, maxit=0, print=F)
# Check
initial2$method
# exclude form the imputation
initial2$method["Sample"]<-"none"
initial2$method["Gene1"]<-"none"
initial2$method["TechnicalReplicate"]<-"none"
initial2$method["BiologicalReplicate"]<-"none"
initial2$method["Comment"]<-"none"

#Check again
initial2$method

# convert to factors
dataWORKING$Sex <- as.factor(dataWORKING$Sex)
dataWORKING$Subpopulation <- as.factor(dataWORKING$Subpopulation)
dataWORKING$Gene <- as.factor(dataWORKING$Gene)
dataWORKING$Category <- as.factor(dataWORKING$Category) # check the effect of ordered
dataWORKING$Recovery <- as.factor(dataWORKING$Recovery) # check the effect of ordered
dataWORKING$Timepoint <- factor(dataWORKING$Timepoint, ordered = TRUE)

# Matrix of predictors:
initial2$predictorMatrix

#### Impute Gene expression --------------------------------------------------------------------------------------------
data_imputed<-mice(dataWORKING, m=20, maxit=10, seed=1234, meth=initial2$method, pred=initial2$predictorMatrix)

# extract the imputed data into a “normal” data frame, run the analyses on each imputation separately and do model diagnostics.
finaldata <- complete(data_imputed, "long")
#View(finaldata)
table(finaldata$.imp)

data_imputed_mean <- finaldata  %>%
    group_by(SampleID, Gene, Subpopulation,Timepoint, TechnicalReplicate, BiologicalReplicate, Sample, Gene1, Age, Sex, Category, Recovery) %>%    
    summarize(CT = mean(CT_Value, na.rm = TRUE),
             # CT_sd = sd(CT_Value, na.rm = TRUE), 
             # CT_N = sum(!is.na(CT_Value)), 
             # dCT_B2M = mean(dCT_Value, na.rm = TRUE), 
             # dCT_B2M_sd = sd(dCT_Value, na.rm = TRUE),
              .groups = "drop")

dataFINAL <- as.data.frame(data_imputed_mean)

# *Normalize the data vs the  Geometric mean of each Sample* ----------------

# Calculate the CT median or mean
dataFINAL$CTmedian <- NA
#dataFINAL$CTmean <- NA
#dataFINAL$CTgmean <- NA  # For geometric mean

dataFINAL <- dataFINAL %>%
    group_by(Sample) %>%
    mutate(
       CTmedian = median(CT[CT < 35], na.rm = TRUE),
       #CTmean = mean(CT[CT < 35], na.rm = TRUE),
       #CTgmean = exp(mean(log(CT[CT < 35]), na.rm = TRUE))  # Geometric mean
    )
# normalize all CT value to the Media CT value of each Sample
dataFINAL$dCT_median <- dataFINAL$CT - dataFINAL$CTmedian
#dataFINAL$dCT_mean <- dataFINAL$CT - dataFINAL$CTmean
#dataFINAL$dCT_gmean <- dataFINAL$CT - dataFINAL$CTgmean

# *Zscore: Normalize the data vs Gene expression ----------------
# as such: apply( data_summary , 1 , function(x) { (x["mean_Z"]- x["MeanForGene"])/ x["sd_for_gene"] } )
dataFINAL <- dataFINAL %>%
    group_by(Gene) %>%
    mutate(
        mean_dCT_Gene = mean(dCT_median, na.rm = TRUE),  # Mean of mean_capZ by Gene
        sd_dCT_Gene = sd(dCT_median, na.rm = TRUE),      # Standard deviation of mean_capZ by Gene
        z_score_dCT = (dCT_median - mean_dCT_Gene) / sd_dCT_Gene,  # Z-score calculation
        inverted_z_score = -z_score_dCT,
        capped_z_score = pmin(pmax(inverted_z_score, -3), 3) 
        ) %>%
    ungroup() 

#### Gene & NHISS Metadata  ----------------
## Gene Metadata & other edits:
Gene_Info <- data.frame(Gene = Metadata_GeneID$Gene, GeneID = Metadata_GeneID$GeneID, Gene_Group = Metadata_GeneID$Gene_Group)
dataFINAL <- merge(dataFINAL, Gene_Info, by = "Gene")
# add the metadata of NHISS score of Patients
dataFINAL <- merge(dataFINAL, Metadata_NHISS, by = c("SampleID", "Timepoint"), all.x = TRUE)

# List of Patients (2x 15))
patients <- unique(dataFINAL$SampleID)
no_patients <- length(patients)

#### Days Post-Stroke data ----------------
# DaysPS <- read_excel("Metadata_DaysPS.xls")
dataFINAL$DaysPS <- as.numeric(NA)
dataFINAL$DaysPS[grep("TP0",dataFINAL$Timepoint)]  <- 0
dataFINAL$DaysPS[grep("TP1",dataFINAL$Timepoint)]  <- 1
dataFINAL$DaysPS[grep("TP2",dataFINAL$Timepoint)]  <- 4
dataFINAL$DaysPS[grep("TP3",dataFINAL$Timepoint)]  <- 30
dataFINAL$DaysPS[grep("TP4",dataFINAL$Timepoint)]  <- 90
#dataFINAL$DaysPS[grep("TP5",dataFINAL$Timepoint)]  <- 360


## MEAN values for the ANOVA & so ----------------
dataFINALmean <- dataFINAL %>%
        group_by(SampleID, Sex, Age, Category, Timepoint, DaysPS,
                 Subpopulation, Gene, GeneID, Gene_Group, Recovery, NHISS) %>%
        summarise(mean_Z = mean(inverted_z_score), sd_Z = sd(inverted_z_score),
                  mean_capZ = mean(capped_z_score), sd_capZ = sd(capped_z_score)
                  )%>%
  ungroup()

#### File Location ----------------
today <- Sys.Date()
output_location <- paste(today,"_Stroke_Results", sep="")
setwd(paste("/Users/ju5263ta/Github/Monocytes/Data/",output_location, "/",sep=""))
getwd()

write_xlsx(fails, "Fail_rates.xlsx")
write_xlsx(dataFINAL, "dataFINAL.xlsx")
write_xlsx(dataFINALmean, "dataFINALmean.xlsx")

# *Dotplot of expressed Genes* -----------------------------------------------------------------------------
createFolder("Dotplots_Gene_summary")

# Create a new column for Timepoint and Subpopulation combination
dataFINALmean <- dataFINALmean %>%
    mutate(Sample_Combo = paste(Timepoint, Subpopulation, sep = "_"))

dataFINALmean <- dataFINALmean %>%
    mutate(Sample_Combo2 = paste(Subpopulation, Timepoint, sep = "_"))

# Create a combined column for Gene and GeneID for better sorting
dataFINALmean <- dataFINALmean %>%
    mutate(Gene_Combined = paste(GeneID, Gene, sep = "_"))  # Combine GeneID and Gene

# Using the function for data_summary and data_summary2
data_summary <- calculate_summary(dataFINALmean, "Sample_Combo", "Gene_Combined")
data_summary2 <- calculate_summary(dataFINALmean, "Sample_Combo2", "Gene_Combined")

data_summary <- data_summary %>%
    mutate(Gene_Combined = factor(Gene_Combined, levels = rev(unique(Gene_Combined))))

data_summary2 <- data_summary2 %>%
    mutate(Gene_Combined = factor(Gene_Combined, levels = rev(unique(Gene_Combined))))

# Create a linear gradient function for the original plots
linear_gradient <- function() {
    # Define the colors for the gradient
    colors <- c("blue", "white", "red")
    
    # Set the breakpoints for the gradient (-1.2 to 1.2)
    scale_color_gradientn(colors = colors,
                          limits = c(-1.2, 1.2),  # Set the limits from -1.2 to 1.2
                          guide = "colorbar",
                          na.value = "grey50")  # Color for NA values
}

ggplot(data_summary, aes(x = Sample_Combo, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_capZ)) +  # Use Dot_Size for scaling
    scale_size(range = c(1, 5)) +  # Size range reflects 1-5 scaling
    linear_gradient() +  # Use the custom linear gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of qPCR Data (high Expr. in red; low Expr. in blue)",
        color = "Gene Expr. \n(Zscored)",
        size = "high SD (small)\nlow SD (big)"  # Use "\n" for line break
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "white"),  # Set the plot background to white
        panel.background = element_rect(fill = "white"),  # Set the panel background to white
        panel.grid.major = element_line(color = "grey", size = 0.5),  # Optional: adjust grid lines
        panel.grid.minor = element_line(color = "lightgrey", size = 0.25)  # Optional: adjust minor grid lines
    )
ggsave(filename = "Dotplots_Gene_summary/Exp_by_TP.png")

ggplot(data_summary2, aes(x = Sample_Combo2, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_capZ)) +  # Use Dot_Size for scaling
    scale_size(range = c(1, 5)) +  # Size range reflects 1-5 scaling
    linear_gradient() +  # Use the custom linear gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of qPCR Data (high Expr. in red; low Expr. in blue)",
        color = "Gene Expr. \n(Zscored)",
        size = "high SD (small)\nlow SD (big)"  # Use "\n" for line break
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "white"),  # Set the plot background to white
        panel.background = element_rect(fill = "white"),  # Set the panel background to white
        panel.grid.major = element_line(color = "grey", size = 0.5),  # Optional: adjust grid lines
        panel.grid.minor = element_line(color = "lightgrey", size = 0.25)  # Optional: adjust minor grid lines
    )
ggsave(filename = "Dotplots_Gene_summary/Exp_by_Suptype.png")

# Create a linear gradient function for the original plots
linear_gradient <- function() {
    # Define the colors for the gradient
    colors <- c("white", "grey", "black")
    
    # Set the breakpoints for the gradient (-1.2 to 1.2)
    scale_color_gradientn(colors = colors,
                          limits = c(-1.2, 1.2),  # Set the limits from -1.2 to 1.2
                          guide = "colorbar",
                          na.value = "red")  # Color for NA values
}

ggplot(data_summary, aes(x = Sample_Combo, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_capZ)) +  # Use Dot_Size for scaling
    scale_size(range = c(1, 5)) +  # Size range reflects 1-5 scaling
    linear_gradient() +  # Use the custom linear gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of qPCR Data (high Expr. in red; low Expr. in blue)",
        color = "Gene Expr. \n(Zscored)",
        size = "high SD (small)\nlow SD (big)"  # Use "\n" for line break
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "white"),  # Set the plot background to white
        panel.background = element_rect(fill = "white"),  # Set the panel background to white
        panel.grid.major = element_line(color = "grey", size = 0.5),  # Optional: adjust grid lines
        panel.grid.minor = element_line(color = "lightgrey", size = 0.25)  # Optional: adjust minor grid lines
    )
ggsave(filename = "Dotplots_Gene_summary/Exp_by_TP_grey.png")

ggplot(data_summary2, aes(x = Sample_Combo2, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_capZ)) +  # Use Dot_Size for scaling
    scale_size(range = c(1, 5)) +  # Size range reflects 1-5 scaling
    linear_gradient() +  # Use the custom linear gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of qPCR Data (high Expr. in red; low Expr. in blue)",
        color = "Gene Expr. \n(Zscored)",
        size = "high SD (small)\nlow SD (big)"  # Use "\n" for line break
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "white"),  # Set the plot background to white
        panel.background = element_rect(fill = "white"),  # Set the panel background to white
        panel.grid.major = element_line(color = "grey", size = 0.5),  # Optional: adjust grid lines
        panel.grid.minor = element_line(color = "lightgrey", size = 0.25)  # Optional: adjust minor grid lines
    )
ggsave(filename = "Dotplots_Gene_summary/Exp_by_Suptype_grey.png")

#### PCA -------
dataFINALmean <- dataFINALmean %>%
    mutate(PCA_Sample = paste(Timepoint, Subpopulation, SampleID, sep = "_"))

# Reverse melt to create matrix for PCA
data_wide <- dcast(dataFINALmean, PCA_Sample ~ Gene, value.var = "mean_capZ")
rownames(data_wide) <- data_wide$PCA_Sample
data_wide <- data_wide[, -1]  # Remove the Sample column to leave only the matrix

data.meta <- t( data.frame(stringr::str_split( rownames( data_wide), "_")))
colnames(data.meta) <- c("Timepoint", "Subpopulation", "SampleID")
data.meta <- merge(data.meta , metadataP, by = "SampleID", all.x = TRUE)
#data.meta <- cbind( data.meta,  dataFINALmean[rownames(data.meta), "Age"])

pca_result <- prcomp(data_wide, scale. = TRUE)

plot_my_pca(pca_result, data.meta, "Timepoint", my_colors, legend_pos = "topleft")

subtype_col <- c("grey", "red", "orange", "green")
plot_my_pca(pca_result, data.meta, "Subpopulation", subtype_col, "topleft")
Category_col <- c("red", "orange", "green")
plot_my_pca(pca_result, data.meta, "Category", Category_col, "topleft")
plot_my_pca(pca_result, data.meta, "Recovery", c("red", "darkgreen","grey"), "topleft")
plot_my_pca(pca_result, data.meta, "Sex", c("#A67C00", "#1D04C2"), "topleft")

data.meta$Age <- factor(data.meta$Age, levels = sort(unique(data.meta$Age)), ordered = TRUE)
# Define Gradient_colour based on the levels of Age
Gradient_colour <- setNames(
    colorRampPalette(c("blue", "red"))(length(levels(data.meta$Age))), 
    levels(data.meta$Age)
)

# Run the function to plot PCA, mapping Age levels to colors
plot_my_pca(pca_result, data.meta, "Age", Gradient_colour[data.meta$Age])


# *FACS Statistics* ----------------
#### Linear regression -------------------------------------
output_folder <- "LinReg_FACS_Plots"
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}

# Define custom colors for the Timepoints
my_colors <- c("TP0" = "darkgreen", "TP1" = "orange", "TP2" = "red", "TP3" = "magenta", "TP4" = "purple")

# Initialize data frames to store results
Age_correlation_FACS <- data.frame(Cells = character(), Timepoint = character(), Sex = character(), N = numeric(), 
                                   p.value = numeric(), Coefficient = numeric(), Down95 = numeric(), 
                                   Up95 = numeric(), Intercept = numeric(), Fold = numeric(), 
                                   Max_age = numeric(), Min_age = numeric())
# remove Patient40 as I do not have age and Sex from them
FACSdata <- FACSdata[!is.na(FACSdata$Age), ]

# Main loop for each column
for (i in 5:17) {
    cells_colname <- colnames(FACSdata)[i]
    
    # All samples
    plot_linear_regression(FACSdata, cells_colname, "All", output_folder)
    
    # Male samples
    FACSdata_male <- FACSdata[FACSdata$Sex == "Male", ]
    if (nrow(FACSdata_male) > 0) {
        plot_linear_regression(FACSdata_male, cells_colname, "Male", output_folder)
    }
    
    # Female samples
    FACSdata_female <- FACSdata[FACSdata$Sex == "Female", ]
    if (nrow(FACSdata_female) > 0) {
        plot_linear_regression(FACSdata_female, cells_colname, "Female", output_folder)
    }
    
    # Minor category
    FACSdata_Minor <- FACSdata[FACSdata$Category == "MINOR", ]
    if (nrow(FACSdata_Minor) > 0) {
        plot_linear_regression(FACSdata_Minor, cells_colname, "Minor", output_folder)
    }
    
    # Moderate category
    FACSdata_Moderate <- FACSdata[FACSdata$Category == "MODERATE", ]
    if (nrow(FACSdata_Moderate) > 0) {
        plot_linear_regression(FACSdata_Moderate, cells_colname, "Moderate", output_folder)
    }
    
    # Good recovery
    FACSdata_Good <- FACSdata[FACSdata$Recovery == "Good", ]
    FACSdata_Good <- FACSdata_Good[!is.na(FACSdata_Good$Age), ]
    if (nrow(FACSdata_Good) > 0) {
        plot_linear_regression(FACSdata_Good, cells_colname, "Good", output_folder)
    }
    
    # Bad recovery
    FACSdata_Bad <- FACSdata[FACSdata$Recovery == "Bad", ]
    FACSdata_Bad <- FACSdata_Bad[!is.na(FACSdata_Bad$Age), ]
    if (nrow(FACSdata_Bad) > 0) {
        plot_linear_regression(FACSdata_Bad, cells_colname, "Bad", output_folder)
    }
}

# Save results to CSV
write.csv(Age_correlation_FACS, "Age_correlation_FACS.csv", row.names = FALSE)

sigAge_correlation_FACS <- Age_correlation_FACS %>% filter(Age_correlation_FACS$p.value <0.05)
write.csv(sigAge_correlation_FACS, "sigAge_correlation_FACS.csv", row.names = FALSE)

#### NHISS linear regression -----------------
NHISS_correlation_FACS <- data.frame(Cells = character(), Timepoint = character(), Sex = character(), N = numeric(), 
                                   p.value = numeric(), Coefficient = numeric(), Down95 = numeric(), 
                                   Up95 = numeric(), Intercept = numeric(), Fold = numeric(), 
                                   Max_NHISS = numeric(), Min_NHISS = numeric())
# remove Patient40 as I do not have NHISS and Sex from them
FACSdata <- FACSdata[!is.na(FACSdata$NHISS), ]

# Main loop for each column
for (i in 5:17) {
    cells_colname <- colnames(FACSdata)[i]
    
    # All samples
    plot_linear_regression_NHISS(FACSdata, cells_colname, "All", output_folder)
    
    # Male samples
    FACSdata_male <- FACSdata[FACSdata$Sex == "Male", ]
    if (nrow(FACSdata_male) > 0) {
        plot_linear_regression_NHISS(FACSdata_male, cells_colname, "Male", output_folder)
    }
    
    # Female samples
    FACSdata_female <- FACSdata[FACSdata$Sex == "Female", ]
    if (nrow(FACSdata_female) > 0) {
        plot_linear_regression_NHISS(FACSdata_female, cells_colname, "Female", output_folder)
    }
    
    # Minor category
    FACSdata_Minor <- FACSdata[FACSdata$Category == "MINOR", ]
    if (nrow(FACSdata_Minor) > 0) {
        plot_linear_regression_NHISS(FACSdata_Minor, cells_colname, "Minor", output_folder)
    }
    
    # Moderate category
    FACSdata_Moderate <- FACSdata[FACSdata$Category == "MODERATE", ]
    if (nrow(FACSdata_Moderate) > 0) {
        plot_linear_regression_NHISS(FACSdata_Moderate, cells_colname, "Moderate", output_folder)
    }
    
    # Good recovery
    FACSdata_Good <- FACSdata[FACSdata$Recovery == "Good", ]
    FACSdata_Good <- FACSdata_Good[!is.na(FACSdata_Good$Age), ]
    if (nrow(FACSdata_Good) > 0) {
        plot_linear_regression_NHISS(FACSdata_Good, cells_colname, "Good", output_folder)
    }
    
    # Bad recovery
    FACSdata_Bad <- FACSdata[FACSdata$Recovery == "Bad", ]
    FACSdata_Bad <- FACSdata_Bad[!is.na(FACSdata_Bad$Age), ]
    if (nrow(FACSdata_Bad) > 0) {
        plot_linear_regression_NHISS(FACSdata_Bad, cells_colname, "Bad", output_folder)
    }
}

# Save results to CSV
write.csv(NHISS_correlation_FACS, "NHISS_correlation_FACS.csv", row.names = FALSE)

sigNHISS_correlation_FACS <- NHISS_correlation_FACS %>% filter(NHISS_correlation_FACS$p.value <0.05)
write.csv(sigNHISS_correlation_FACS, "sigNHISS_correlation_FACS.csv", row.names = FALSE)

#### 1WAY ANOVA* ----------------------------------------------------------------------------------------------------
# with wilcox
ANOVA_FACSall <- automate_anova_extraction(output_location, "Wilcox_FACS_Plots", 
                                                  "Wilcox_FACS", "y", FACSdata, colnames(FACSdata)[5:17], 
                                                  "Timepoint")
# paired t.test
# FACSdata_matched <- FACSdata %>%
#     filter(!SampleID %in% Unmatched_TP0_FACS)
# 
# paired_Ttest_FACS <- automate_ttest_extraction(output_location, "paired_Ttest_FACS_Plots",
#                           "paired_Ttest_FACS", "y", FACSdata_matched, colnames(FACSdata_matched)[4:16],
#                           "Timepoint")

#### Mean & SEM FACS----------------------------------------------------------------------------------------------------
folder <- "MEAN_SEM_FACS_Sex"
createFolder(folder)
plot_save <- "y"

# Initialize an empty dataframe to store Wilcoxon test results
results_FACS_Wilcox <- data.frame(Cells = character(), Focus = character(), 
                                  TP0_P_Value = numeric(), TP1_P_Value = numeric(), 
                                  TP2_P_Value = numeric(), TP3_P_Value = numeric(), 
                                  TP4_P_Value = numeric(), stringsAsFactors = FALSE)

for (i in 5:17) {
    cells_colname <- colnames(FACSdata)[i]
    
    # Filter, clean, and rename the data column for the current column
    df_subpop_clean <- FACSdata %>%
        filter(!is.na(Sex), !is.na(Timepoint)) %>%
        select(Timepoint, Sex, !!sym(cells_colname)) %>%
        rename(data_column = !!sym(cells_colname))
    
    p_values <- c(TP0_P_Value = NA, TP1_P_Value = NA, TP2_P_Value = NA, 
                  TP3_P_Value = NA, TP4_P_Value = NA)
    
    # Calculate Wilcoxon test p-values for each timepoint
    for (tp in unique(df_subpop_clean$Timepoint)) {
        df_timepoint <- subset(df_subpop_clean, Timepoint == tp)
        
        if (length(unique(df_timepoint$Sex)) == 2) { # Ensure both sexes are present
            wilcox_test <- wilcox.test(data_column ~ Sex, data = df_timepoint)
            p_values[paste0(tp, "_P_Value")] <- wilcox_test$p.value
        }
    }
    
    # Append p-values to results dataframe
    results_FACS_Wilcox <- rbind(results_FACS_Wilcox, 
                                 data.frame(Cells = cells_colname, Focus = "Male vs. Female", 
                                            t(p_values), stringsAsFactors = FALSE))
    
    # Generate line plot with mean ± SEM for each sex across Timepoints
    p <- ggline(df_subpop_clean, 
                x = "Timepoint", 
                y = "data_column", 
                color = "Sex", 
                add = "mean_se",  # Add mean and standard error
                size = 1.2) +  
        labs(y = paste(cells_colname, "percentage"), 
             title = paste("Mean ± SEM of", cells_colname, "by Sex across Timepoints")) +
        scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2")) +
        theme_bw() +
        theme(
            panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.line = element_line(color = "black")
        ) +
        stat_compare_means(method = "wilcox.test",  # or "t.test"
                           aes(group = Sex),
                           label = "p.signif",  
                           size = 5) 
    
    # Save plot if specified
    if (plot_save == "y") {
        file_name_facet <- paste0(folder, "/", cells_colname, "_Mean_SEM_SexComparison.png")
        ggsave(filename = file_name_facet, plot = p)
    }
}


folder <- "MEAN_SEM_FACS_Category"
createFolder(folder)
plot_save <- "y"

for (i in 5:17) {
    cells_colname <- colnames(FACSdata)[i]
    
    # Filter, clean, and rename the data column for the current column
    df_subpop_clean <- FACSdata %>%
        filter(!is.na(Category), !is.na(Timepoint)) %>%
        select(Timepoint, Category, !!sym(cells_colname)) %>%
        rename(data_column = !!sym(cells_colname))
    
    p_values <- c(TP0_P_Value = NA, TP1_P_Value = NA, TP2_P_Value = NA, 
                  TP3_P_Value = NA, TP4_P_Value = NA)
    
    # Calculate Wilcoxon test p-values for each timepoint
    for (tp in unique(df_subpop_clean$Timepoint)) {
        df_timepoint <- subset(df_subpop_clean, Timepoint == tp)
        
        if (length(unique(df_timepoint$Category)) == 2) { # Ensure both sexes are present
            wilcox_test <- wilcox.test(data_column ~ Category, data = df_timepoint)
            p_values[paste0(tp, "_P_Value")] <- wilcox_test$p.value
        }
    }
    
    # Append p-values to results dataframe
    results_FACS_Wilcox <- rbind(results_FACS_Wilcox, 
                                 data.frame(Cells = cells_colname, Focus = "Minor vs. Moderate", 
                                            t(p_values), stringsAsFactors = FALSE))
    
    # Generate line plot with mean ± SEM for each sex across Timepoints
    p <- ggline(df_subpop_clean, 
                x = "Timepoint", 
                y = "data_column", 
                color = "Category", 
                add = "mean_se",  # Add mean and standard error
                size = 1.2) +  
        labs(y = paste(cells_colname, "percentage"), 
             title = paste("Mean ± SEM of", cells_colname, "by Category across Timepoints")) +
        scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606")) +
        theme_bw() +
        theme(
            panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.line = element_line(color = "black")
        ) +
        stat_compare_means(method = "wilcox.test",  # or "t.test"
                           aes(group = Category),
                           label = "p.signif",  
                           size = 5) 
    
    # Save plot if specified
    if (plot_save == "y") {
        file_name_facet <- paste0(folder, "/", cells_colname, "_Mean_SEM_CategoryComparison.png")
        ggsave(filename = file_name_facet, plot = p)
    }
}

folder <- "MEAN_SEM_FACS_Recovery"
createFolder(folder)
plot_save <- "y"

for (i in 5:17) {
    cells_colname <- colnames(FACSdata)[i]
    
    # Filter, clean, and rename the data column for the current column
    df_subpop_clean <- FACSdata %>%
        filter(!is.na(Recovery), !is.na(Timepoint)) %>%
        select(Timepoint, Recovery, !!sym(cells_colname)) %>%
        rename(data_column = !!sym(cells_colname))
    
    p_values <- c(TP0_P_Value = NA, TP1_P_Value = NA, TP2_P_Value = NA, 
                  TP3_P_Value = NA, TP4_P_Value = NA)
   
    for (tp in unique(df_subpop_clean$Timepoint)) {
        df_timepoint <- subset(df_subpop_clean, Timepoint == tp)
        
        if (length(unique(df_timepoint$Recovery)) == 2) { # Ensure both sexes are present
            wilcox_test <- wilcox.test(data_column ~ Recovery, data = df_timepoint)
            p_values[paste0(tp, "_P_Value")] <- wilcox_test$p.value
        }
    }
    
    results_FACS_Wilcox <- rbind(results_FACS_Wilcox, 
                                 data.frame(Cells = cells_colname, Focus = "Good vs. Bad Recovery", 
                                            t(p_values), stringsAsFactors = FALSE))
    p <- ggline(df_subpop_clean, 
                x = "Timepoint", 
                y = "data_column", 
                color = "Recovery", 
                add = "mean_se",  # Add mean and standard error
                size = 1.2) +  
        labs(y = paste(cells_colname, "percentage"), 
             title = paste("Mean ± SEM of", cells_colname, "by Recovery across Timepoints")) +
        scale_color_manual(values = c("Good" = "darkgreen", "okay"= "#A67C00", "Bad" = "#700606")) +
        theme_bw() +
        theme(
            panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.line = element_line(color = "black")
        ) +
        stat_compare_means(method = "wilcox.test",  # or "t.test"
                           aes(group = Recovery),
                           label = "p.signif",  
                           size = 5) 
    
    # Save plot if specified
    if (plot_save == "y") {
        file_name_facet <- paste0(folder, "/", cells_colname, "_Mean_SEM_RecoveryComparison.png")
        ggsave(filename = file_name_facet, plot = p)
    }
}

results_FACS_Wilcox$TP0_Significance <- sapply(results_FACS_Wilcox$TP0_P_Value, get_significance)
results_FACS_Wilcox$TP1_Significance <- sapply(results_FACS_Wilcox$TP1_P_Value, get_significance)
results_FACS_Wilcox$TP2_Significance <- sapply(results_FACS_Wilcox$TP2_P_Value, get_significance)
results_FACS_Wilcox$TP3_Significance <- sapply(results_FACS_Wilcox$TP3_P_Value, get_significance)
results_FACS_Wilcox$TP4_Significance <- sapply(results_FACS_Wilcox$TP4_P_Value, get_significance)
write.csv(results_FACS_Wilcox, file = "results_FACS_Wilcox.csv", row.names = FALSE)

# *Gene Expr. Statistics* ----------------------------------------------------------------------------------------------------
# As they are non-normal disributed, dependent with similar variance
# Adjust to yes if you want plots to be saved or n if not
plot_save <- "y"

# remove unmatched CTR for ANOVA too 
data_mean_matched <- dataFINALmean %>% filter(!(Timepoint == "TP0" & !SampleID %in% Matched_TP0_Gene))
# Check if all groups have the same size:
check_sample_counts(data_mean_matched)

TP_colors <- c("TP0" = "darkgreen", "TP1" = "orange", "TP2" = "red", "TP3" = "magenta", "TP4" = "purple")

#### Wilcox test for Z of Genes --------------------------------------------------
folder <- "Wilcox_Z_Plots_unpaired"
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
        plot_save <- "n"
        if (p_value_min <= 0.05) {
            plot_save <- "y"
        }
        
        if (plot_save == "y") {
            #my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
            
            # Create a boxplot with significance levels
            p <- ggboxplot(df, 
                           x = "Timepoint", 
                           y = "mean_capZ", 
                           color = "Timepoint", 
                           add = "jitter", 
                           legend = "none", 
                           ylab = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                           width = 0.8, 
                           add.params = list(size = 1, alpha = 0.5)) +  
                geom_hline(yintercept = median_TP0, linetype = 2) +
                scale_color_manual(values = TP_colors)
            
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

#### Wilcox test for Z of Genes --------------------------------------------------
folder <- "Wilcox_Z_Plots_SubtypesCombined_unpaired"
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


#*Wilcox for Comparison at individual TPs -------------------------------------

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
folder <- "Mean_SEM_Category_capZ_Plots"
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
folder <- "Mean_SEM_Recovery_capZ_Plots"
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
                ggsave(filename = file_name_facet)
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

# *Age Regression* --------------------------------------------------------
#### Pearson Gene vs Age --------------------------------------------------------
createFolder("LinReg_Results_capZ")
createFolder("LinReg_Plots_Age_capZ")

# Pearson's Correlation uses linear relationship to correlate the data
# Initialize dataframes to store correlation results
Age_correlation_capZ <- data.frame(Subpopulation = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

Age_correlation_Sex_capZ <- data.frame(Subpopulation = character(), Sex = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

Age_correlation_Category_capZ <- data.frame(Subpopulation = character(), Category = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

for (gene in unique(dataFINALmean$Gene)) {
    df <- dataFINALmean %>% filter(Gene == gene)
    df <- df %>% filter(!is.na(Timepoint), !is.na(mean_capZ), !is.na(Subpopulation), !is.na(Sex), !is.na(Category))
    
    for (subpop in unique(df$Subpopulation)) {
        df_subpop <- df %>% filter(Subpopulation == subpop)
        
        for (tp in unique(df_subpop$Timepoint)) {
            df_tp <- df_subpop %>% filter(Timepoint == tp)
            if (nrow(df_tp) <= 0) next
            
            # Outlier removal
            lm_initial <- lm(mean_capZ ~ Age, data = df_tp)
            residuals_values <- residuals(lm_initial)
            outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
            # Remove outliers only if any are detected
            if (length(outlier_indices) > 0) {
                df_tp_filtered <- df_tp[-outlier_indices, ]
            } else {
                df_tp_filtered <- df_tp  # No outliers, keep the original dataframe
            }
            
            if (nrow(df_tp_filtered) < 2) next  # Ensure there are enough data points after outlier removal
            
            tryCatch({
                # Pearson's correlation
                x <- cor.test(df_tp_filtered$mean_capZ, df_tp_filtered$Age)
                row <- data.frame(
                    Subpopulation = subpop,
                    Timepoint = tp,
                    Gene = gene,
                    GeneID = df_tp_filtered$GeneID[1],
                    Gene_Group = df_tp_filtered$Gene_Group[1],
                    N = nrow(df_tp_filtered),
                    p.value = x$p.value,
                    Estimate = x$estimate,
                    Est_CI_Lower = x$conf.int[1],
                    Est_CI_Upper = x$conf.int[2],
                    Coefficient = NA,
                    Coef_CI_Lower = NA,
                    Coef_CI_Upper = NA  
                )
                
                # Linear regression for coefficient and CI
                lm_result <- lm(mean_capZ ~ Age, data = df_tp_filtered)
                coef_x <- summary(lm_result)$coefficients[2, 1]
                ci <- confint(lm_result)[2, ]
                
                row$Coefficient <- coef_x
                row$Coef_CI_Lower <- ci[1]
                row$Coef_CI_Upper <- ci[2]
                
                if (row$p.value<=0.05){
                    titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," (All Samples)", sep="")
                    file_name <- paste("LinReg_Plots_Age_capZ/",titel, ".png", sep = "")
                    ggplot(data = df_tp_filtered, aes(x = Age, y = mean_capZ)) +
                        geom_smooth(method = "glm", color = "black") +
                        geom_point(aes(color = Timepoint), size = 2) +
                        ggtitle(titel) +
                        xlab("Age") +
                        ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                        scale_color_manual(values = my_colors) +
                        theme(text = element_text(size=14)) +
                        theme(
                            panel.background = element_rect(fill = "white", color = "black"),
                            plot.background = element_rect(fill = "white", color = "black"),
                            text = element_text(color = "black"),
                            panel.grid.major = element_line(color = "gray", size = 0.2),
                            panel.grid.minor = element_blank(),
                            panel.grid.major.x = element_blank()
                        )
                    ggsave(filename=file_name)
                }
                Age_correlation_capZ <- rbind(Age_correlation_capZ, row)
            }, error = function(e) {
                cat("Error in Pearson's Correlation for gene", gene, "& Subpopulation with Age", subpop, "- skipping this comparison.\n")
            })
            
            for (sex in unique(df_tp$Sex)) {
                df_sex <- df_tp %>% filter(Sex == sex)
                # Outlier removal
                lm_initial <- lm(mean_capZ ~ Age, data = df_sex)
                residuals_values <- residuals(lm_initial)
                outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
                # Remove outliers only if any are detected
                if (length(outlier_indices) > 0) {
                    df_sex_filtered <- df_sex[-outlier_indices, ]
                } else {
                    df_sex_filtered <- df_sex  # No outliers, keep the original dataframe
                }
                
                if (nrow(df_sex_filtered) < 2) next  # Ensure there are enough data points after outlier removal
                
                
                tryCatch({
                    # Pearson's correlation for Sex
                    x <- cor.test(df_sex_filtered$mean_capZ, df_sex_filtered$Age)
                    row <- data.frame(
                        Subpopulation = subpop,
                        Sex = sex,
                        Timepoint = tp,
                        Gene = gene,
                        GeneID = df_sex_filtered$GeneID[1],
                        Gene_Group = df_sex_filtered$Gene_Group[1],
                        N = nrow(df_sex_filtered),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA,
                        Coef_CI_Lower = NA,
                        Coef_CI_Upper = NA  
                    )
                    
                    # Linear regression for coefficient and CI
                    lm_result_sex <- lm(mean_capZ ~ Age, data = df_sex_filtered)
                    coef_x <- summary(lm_result_sex)$coefficients[2, 1]
                    ci <- confint(lm_result_sex)[2, ]
                    
                    row$Coefficient <- coef_x
                    row$Coef_CI_Lower <- ci[1]
                    row$Coef_CI_Upper <- ci[2]
                    
                    if (row$p.value<=0.05){
                        titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," ( ", sex, " Samples)", sep="")
                        file_name <- paste("LinReg_Plots_Age_capZ/",titel, ".png", sep = "")
                        ggplot(data = df_sex_filtered, aes(x = Age, y = mean_capZ)) +
                            geom_smooth(method = "glm", color = "black") +
                            geom_point(aes(color = Sex), size = 2) +
                            ggtitle(titel) +
                            xlab("Age") +
                            ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                            scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2")) +
                            theme(text = element_text(size=14)) +
                            theme(
                                panel.background = element_rect(fill = "white", color = "black"),
                                plot.background = element_rect(fill = "white", color = "black"),
                                text = element_text(color = "black"),
                                panel.grid.major = element_line(color = "gray", size = 0.2),
                                panel.grid.minor = element_blank(),
                                panel.grid.major.x = element_blank()
                            )
                        ggsave(filename=file_name)
                    }
                    Age_correlation_Sex_capZ <- rbind(Age_correlation_Sex_capZ, row)
                }, error = function(e) {
                    cat("Error for", gene, tp, sex, "&", subpop, "- skipping this comparison.\n")
                }) 
            }
            for (cat in unique(df_tp$Category)) {
                df_cat <- df_tp %>% filter(Category == cat)
                # Outlier removal
                lm_initial <- lm(mean_capZ ~ Age, data = df_cat)
                residuals_values <- residuals(lm_initial)
                outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
                # Remove outliers only if any are detected
                if (length(outlier_indices) > 0) {
                    df_cat_filtered <- df_cat[-outlier_indices, ]
                } else {
                    df_cat_filtered <- df_cat  # No outliers, keep the original dataframe
                }
                
                if (nrow(df_cat_filtered) < 2) next  # Ensure there are enough data points after outlier removal
                
                tryCatch({
                    # Pearson's correlation for Category
                    x <- cor.test(df_cat_filtered$mean_capZ, df_cat_filtered$Age)
                    row <- data.frame(
                        Subpopulation = subpop,
                        Category = cat,
                        Timepoint = tp,
                        Gene = gene,
                        GeneID = df_cat_filtered$GeneID[1],
                        Gene_Group = df_cat_filtered$Gene_Group[1],
                        N = nrow(df_cat_filtered),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA,
                        Coef_CI_Lower = NA,
                        Coef_CI_Upper = NA  
                    )
                    
                    # Linear regression for coefficient and CI
                    lm_result_cat <- lm(mean_capZ ~ Age, data = df_cat_filtered)
                    coef_x <- summary(lm_result_cat)$coefficients[2, 1]
                    ci <- confint(lm_result_cat)[2, ]
                    
                    row$Coefficient <- coef_x
                    row$Coef_CI_Lower <- ci[1]
                    row$Coef_CI_Upper <- ci[2]
                    
                    if (row$p.value<=0.05){
                        titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," ( ", cat, " Samples)", sep="")
                        file_name <- paste("LinReg_Plots_Age_capZ/",titel, ".png", sep = "")
                        ggplot(data = df_cat_filtered, aes(x = Age, y = mean_capZ)) +
                            geom_smooth(method = "glm", color = "black") +
                            geom_point(aes(color = Category), size = 2) +
                            ggtitle(titel) +
                            xlab("Age") +
                            ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                            scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606")) +
                            theme(text = element_text(size=14)) +
                            theme(
                                panel.background = element_rect(fill = "white", color = "black"),
                                plot.background = element_rect(fill = "white", color = "black"),
                                text = element_text(color = "black"),
                                panel.grid.major = element_line(color = "gray", size = 0.2),
                                panel.grid.minor = element_blank(),
                                panel.grid.major.x = element_blank()
                            )
                        ggsave(filename=file_name)
                    }
                    Age_correlation_Category_capZ <- rbind(Age_correlation_Category_capZ, row)
                }, error = function(e) {
                    cat("Error for", gene, tp, cat, "&", subpop, "- skipping this comparison.\n")
                }) 
            }
        }
    }
}

file_name <- "LinReg_Results_capZ/Pearson's Correlation_Age_capZ.csv"
write.csv(Age_correlation_capZ, file_name, row.names = FALSE)

file_name_s <- "LinReg_Results_capZ/Pearson's Correlation_Age_Sex_capZ.csv"
write.csv(Age_correlation_Sex_capZ, file_name_s, row.names = FALSE)

file_name_s <- "LinReg_Results_capZ/Pearson's Correlation_Age_Category_capZ.csv"
write.csv(Age_correlation_Category_capZ, file_name_s, row.names = FALSE)

#### Extract significant LinRegs --------------------------------------------------------

Age_correlation_capZ$Sex <- "All"
Age_correlation_capZ$Category <- "All"
Age_correlation_Sex_capZ$Category <- "All"
Age_correlation_Category_capZ$Sex<- "All"

All_Age_Pearson_capZ <- rbind(Age_correlation_capZ,Age_correlation_Sex_capZ, Age_correlation_Category_capZ)
Age_sigGenes_Pearson <- All_Age_Pearson_capZ %>% filter(All_Age_Pearson_capZ$p.value <0.05)
write.csv(Age_sigGenes_Pearson, "LinReg_Results_capZ/Age_sigGenes_Pearson.csv", row.names = FALSE)

# *NHISS Regression* --------------------------------------------------------
#### Pearson Gene vs NHISS --------------------------------------------------------
createFolder("LinReg_Plots_NHISS_capZ")

# Pearson's Correlation uses linear relationship to correlate the data
# Initialize dataframes to store correlation results
NHISS_correlation_capZ <- data.frame(Subpopulation = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

NHISS_correlation_Recovery_capZ <- data.frame(Subpopulation = character(), Recovery = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

NHISS_correlation_Category_capZ <- data.frame(Subpopulation = character(), Category = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

for (gene in unique(dataFINALmean$Gene)) {
    df <- dataFINALmean %>% filter(Gene == gene)
    df <- df %>% filter(!is.na(Timepoint), !is.na(mean_capZ), !is.na(Subpopulation), !is.na(Recovery), !is.na(Category))
    
    for (subpop in unique(df$Subpopulation)) {
        df_subpop <- df %>% filter(Subpopulation == subpop)
        
        for (tp in unique(df_subpop$Timepoint)) {
            df_tp <- df_subpop %>% filter(Timepoint == tp, !is.na(NHISS))
            if (nrow(df_tp) <= 0) next
            
            # Outlier removal
            lm_initial <- lm(mean_capZ ~ NHISS, data = df_tp)
            residuals_values <- residuals(lm_initial)
            outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
            # Remove outliers only if any are detected
            if (length(outlier_indices) > 0) {
                df_tp_filtered <- df_tp[-outlier_indices, ]
            } else {
                df_tp_filtered <- df_tp  # No outliers, keep the original dataframe
            }
            
            if (nrow(df_tp_filtered) < 2) next  # Ensure there are enough data points after outlier removal
            
            tryCatch({
                # Pearson's correlation
                x <- cor.test(df_tp_filtered$mean_capZ, df_tp_filtered$NHISS)
                row <- data.frame(
                    Subpopulation = subpop,
                    Timepoint = tp,
                    Gene = gene,
                    GeneID = df_tp_filtered$GeneID[1],
                    Gene_Group = df_tp_filtered$Gene_Group[1],
                    N = nrow(df_tp_filtered),
                    p.value = x$p.value,
                    Estimate = x$estimate,
                    Est_CI_Lower = x$conf.int[1],
                    Est_CI_Upper = x$conf.int[2],
                    Coefficient = NA,
                    Coef_CI_Lower = NA,
                    Coef_CI_Upper = NA  
                )
                
                # Linear regression for coefficient and CI
                lm_result <- lm(mean_capZ ~ NHISS, data = df_tp_filtered)
                coef_x <- summary(lm_result)$coefficients[2, 1]
                ci <- confint(lm_result)[2, ]
                
                row$Coefficient <- coef_x
                row$Coef_CI_Lower <- ci[1]
                row$Coef_CI_Upper <- ci[2]
                
                if (row$p.value<=0.05){
                    titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," (All Samples)", sep="")
                    file_name <- paste("LinReg_Plots_NHISS_capZ/",titel, ".png", sep = "")
                    ggplot(data = df_tp_filtered, aes(x = NHISS, y = mean_capZ)) +
                        geom_smooth(method = "glm", color = "black") +
                        geom_point(aes(color = Timepoint), size = 2) +
                        ggtitle(titel) +
                        xlab("NHISS") +
                        ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                        scale_color_manual(values = my_colors) +
                        theme(text = element_text(size=14)) +
                        theme(
                            panel.background = element_rect(fill = "white", color = "black"),
                            plot.background = element_rect(fill = "white", color = "black"),
                            text = element_text(color = "black"),
                            panel.grid.major = element_line(color = "gray", size = 0.2),
                            panel.grid.minor = element_blank(),
                            panel.grid.major.x = element_blank()
                        )
                    ggsave(filename=file_name)
                }
                NHISS_correlation_capZ <- rbind(NHISS_correlation_capZ, row)
            }, error = function(e) {
                cat("Error in Pearson's Correlation for gene", gene, "& Subpopulation with NHISS", subpop, "- skipping this comparison.\n")
            })
            
            for (recovery in unique(df_tp$Recovery)) {
                df_recovery <- df_tp %>% filter(Recovery == recovery, !is.na(NHISS))
                # Outlier removal
                lm_initial <- lm(mean_capZ ~ NHISS, data = df_recovery)
                residuals_values <- residuals(lm_initial)
                outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
                # Remove outliers only if any are detected
                if (length(outlier_indices) > 0) {
                    df_recovery_filtered <- df_recovery[-outlier_indices, ]
                } else {
                    df_recovery_filtered <- df_recovery  # No outliers, keep the original dataframe
                }
                
                if (nrow(df_recovery_filtered) < 2) next  # Ensure there are enough data points after outlier removal
                
                
                tryCatch({
                    # Pearson's correlation for recovery
                    x <- cor.test(df_recovery_filtered$mean_capZ, df_recovery_filtered$NHISS)
                    row <- data.frame(
                        Subpopulation = subpop,
                        Recovery = recovery,
                        Timepoint = tp,
                        Gene = gene,
                        GeneID = df_recovery_filtered$GeneID[1],
                        Gene_Group = df_recovery_filtered$Gene_Group[1],
                        N = nrow(df_recovery_filtered),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA,
                        Coef_CI_Lower = NA,
                        Coef_CI_Upper = NA  
                    )
                    
                    # Linear regression for coefficient and CI
                    lm_result_recovery <- lm(mean_capZ ~ NHISS, data = df_recovery_filtered)
                    coef_x <- summary(lm_result_recovery)$coefficients[2, 1]
                    ci <- confint(lm_result_recovery)[2, ]
                    
                    row$Coefficient <- coef_x
                    row$Coef_CI_Lower <- ci[1]
                    row$Coef_CI_Upper <- ci[2]
                    
                    if (row$p.value<=0.05){
                        titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," ( ", recovery, " Samples)", sep="")
                        file_name <- paste("LinReg_Plots_NHISS_capZ/",titel, ".png", sep = "")
                        ggplot(data = df_recovery_filtered, aes(x = NHISS, y = mean_capZ)) +
                            geom_smooth(method = "glm", color = "black") +
                            geom_point(aes(color = recovery), size = 2) +
                            ggtitle(titel) +
                            xlab("NHISS") +
                            ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                            scale_color_manual(values = c("Good" = "darkgreen", "Bad" = "#700606")) +
                            theme(text = element_text(size=14)) +
                            theme(
                                panel.background = element_rect(fill = "white", color = "black"),
                                plot.background = element_rect(fill = "white", color = "black"),
                                text = element_text(color = "black"),
                                panel.grid.major = element_line(color = "gray", size = 0.2),
                                panel.grid.minor = element_blank(),
                                panel.grid.major.x = element_blank()
                            )
                        ggsave(filename=file_name)
                    }
                    NHISS_correlation_Recovery_capZ <- rbind(NHISS_correlation_Recovery_capZ, row)
                }, error = function(e) {
                    cat("Error for", gene, tp, recovery, "&", subpop, "- skipping this comparison.\n")
                }) 
            }
            
            for (cat in unique(df_tp$Category)) {
                df_cat <- df_tp %>% filter(Category == cat, !is.na(NHISS))
                # Outlier removal
                lm_initial <- lm(mean_capZ ~ NHISS, data = df_cat)
                residuals_values <- residuals(lm_initial)
                outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
                # Remove outliers only if any are detected
                if (length(outlier_indices) > 0) {
                    df_cat_filtered <- df_cat[-outlier_indices, ]
                } else {
                    df_cat_filtered <- df_cat  # No outliers, keep the original dataframe
                }
                
                if (nrow(df_cat_filtered) < 2) next  # Ensure there are enough data points after outlier removal
                
                tryCatch({
                    # Pearson's correlation for Category
                    x <- cor.test(df_cat_filtered$mean_capZ, df_cat_filtered$NHISS)
                    row <- data.frame(
                        Subpopulation = subpop,
                        Category = cat,
                        Timepoint = tp,
                        Gene = gene,
                        GeneID = df_cat_filtered$GeneID[1],
                        Gene_Group = df_cat_filtered$Gene_Group[1],
                        N = nrow(df_cat_filtered),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA,
                        Coef_CI_Lower = NA,
                        Coef_CI_Upper = NA  
                    )
                    
                    # Linear regression for coefficient and CI
                    lm_result_cat <- lm(mean_capZ ~ NHISS, data = df_cat_filtered)
                    coef_x <- summary(lm_result_cat)$coefficients[2, 1]
                    ci <- confint(lm_result_cat)[2, ]
                    
                    row$Coefficient <- coef_x
                    row$Coef_CI_Lower <- ci[1]
                    row$Coef_CI_Upper <- ci[2]
                    
                    if (row$p.value<=0.05){
                        titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," ( ", cat, " Samples)", sep="")
                        file_name <- paste("LinReg_Plots_NHISS_capZ/",titel, ".png", sep = "")
                        ggplot(data = df_cat_filtered, aes(x = NHISS, y = mean_capZ)) +
                            geom_smooth(method = "glm", color = "black") +
                            geom_point(aes(color = Category), size = 2) +
                            ggtitle(titel) +
                            xlab("NHISS") +
                            ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                            scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606")) +
                            theme(text = element_text(size=14)) +
                            theme(
                                panel.background = element_rect(fill = "white", color = "black"),
                                plot.background = element_rect(fill = "white", color = "black"),
                                text = element_text(color = "black"),
                                panel.grid.major = element_line(color = "gray", size = 0.2),
                                panel.grid.minor = element_blank(),
                                panel.grid.major.x = element_blank()
                            )
                        ggsave(filename=file_name)
                    }
                    NHISS_correlation_Category_capZ <- rbind(NHISS_correlation_Category_capZ, row)
                }, error = function(e) {
                    cat("Error for", gene, tp, cat, "&", subpop, "- skipping this comparison.\n")
                }) 
            }
        }
    }
}

file_name <- "LinReg_Results_capZ/Pearson's Correlation_NHISS_capZ.csv"
write.csv(NHISS_correlation_capZ, file_name, row.names = FALSE)

file_name_s <- "LinReg_Results_capZ/Pearson's Correlation_NHISS_Recovery_capZ.csv"
write.csv(NHISS_correlation_Recovery_capZ, file_name_s, row.names = FALSE)

file_name_s <- "LinReg_Results_capZ/Pearson's Correlation_NHISS_Category_capZ.csv"
write.csv(NHISS_correlation_Category_capZ, file_name_s, row.names = FALSE)

#### Extract significant LinRegs --------------------------------------------------------

NHISS_correlation_capZ$Recovery <- "All"
NHISS_correlation_capZ$Category <- "All"
NHISS_correlation_Recovery_capZ$Category <- "All"
NHISS_correlation_Category_capZ$Recovery<- "All"

All_NHISS_Pearson_capZ <- rbind(NHISS_correlation_capZ,NHISS_correlation_Recovery_capZ, NHISS_correlation_Category_capZ)
NHISS_sigGenes_Pearson <- All_NHISS_Pearson_capZ %>% filter(All_NHISS_Pearson_capZ$p.value <0.05)
write.csv(NHISS_sigGenes_Pearson, "LinReg_Results_capZ/NHISS_sigGenes_Pearson.csv", row.names = FALSE)


# *Heatmaps* ----------------------------------
new_folder <- "Pearsson_Heatmaps"
createFolder(new_folder)

#### Aging (Estimate) ---------------
Age_correlation_capZ <- Age_correlation_capZ %>%
    mutate(Gene_Ordered = paste(GeneID, Gene, sep = "_"))
Age_correlation_capZ <- Age_correlation_capZ %>%
    mutate(Groups = paste(Subpopulation, Timepoint, sep = "_"))

Age_correlation_capZ$Significance <- ifelse(Age_correlation_capZ$p.value < 0.001, "***",
                                              ifelse(Age_correlation_capZ$p.value < 0.01, "**",
                                                     ifelse(Age_correlation_capZ$p.value < 0.05, "*", "")))
# Estimate: - 0.7 to 0.79 
name_plot <- "Heatmap of Age Correlation Estimates"
ggplot(Age_correlation_capZ, aes(Gene_Ordered, Groups, fill= Estimate)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-0.7, 0.804)) +
    geom_text(aes(label = Significance), color = "white", size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, x = "Gene", y = "Timepoint & Suptype", fill = "Estimate")
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
    plot = last_plot(), width = 10, height = 5, dpi = 300)

name_plot <- paste(name_plot," (transposed)")
ggplot(Age_correlation_capZ, aes(Groups, Gene_Ordered, fill= Estimate)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-0.7, 0.804)) +
    geom_text(aes(label = Significance), color = "white", size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, x = "Timepoint & Suptype", y = "Gene", fill = "Estimate")
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 5, height = 10, dpi = 300)

#### NHISS (Estimate) ---------------
NHISS_correlation_capZ <- NHISS_correlation_capZ %>%
    mutate(Gene_Ordered = paste(GeneID, Gene, sep = "_"))
NHISS_correlation_capZ <- NHISS_correlation_capZ %>%
    mutate(Groups = paste(Subpopulation, Timepoint, sep = "_"))

# Estimate: - 0.99 to 0.92
name_plot <- paste("Heatmap of NHISS Correlation Estimates")
# Add significance levels
NHISS_correlation_capZ$Significance <- ifelse(NHISS_correlation_capZ$p.value < 0.001, "***",
                                              ifelse(NHISS_correlation_capZ$p.value < 0.01, "**",
                                                     ifelse(NHISS_correlation_capZ$p.value < 0.05, "*", "")))
# Create the heatmap with significance levels
ggplot(NHISS_correlation_capZ, aes(x = Gene_Ordered, y = Groups, fill = Estimate)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-0.99, 0.92)) +
    geom_text(aes(label = Significance), color = "black", size = 3) + # Add significance levels with white text
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, x = "Gene", y = "Timepoint & Subpopulation", fill = "Estimate")
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 10, height = 5, dpi = 300)

name_plot <- paste(name_plot," (transposed)")
ggplot(NHISS_correlation_capZ, aes(Groups, Gene_Ordered, fill= Estimate)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-0.99, 0.92)) +
    geom_text(aes(label = Significance), color = "black", size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, x = "Timepoint & Subpopulation", y = "Timepoint", fill = "Estimate")
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 5, height = 10, dpi = 300)

#### Aging  Split by Subpopulation (Estimate) ---------------
for (sup in unique(Age_correlation_capZ$Subpopulation)) {
    df <- filter(Age_correlation_capZ, Subpopulation == sup)
    sigGenes <- filter(Age_correlation_capZ, p.value <= 0.05)
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of Aging Correlation Estimates - ", sup, " Monocytes") 
    ggplot(df, aes(Gene_Ordered, Timepoint, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Timepoint", fill = "Coefficient")+
        scale_y_discrete(labels = timepoint_labels)+ # Apply custom y-axis labels
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}

for (sup in unique(Age_correlation_capZ$Subpopulation)) {
    df <- filter(Age_correlation_capZ, Subpopulation == sup)
    # sigGenes <- filter(Age_correlation_capZ, p.value <= 0.05) #if you want to compare all of them
    sigGenes <- filter(df, p.value <= 0.05)
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of Aging Correlation Estimates - ", sup, " Monocytes (only sigGenes)") 
    ggplot(df, aes(Gene_Ordered, Timepoint, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Timepoint", fill = "Coefficient")+
        scale_y_discrete(labels = timepoint_labels)+# Apply custom y-axis labels
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}

#### NHISS Split by Subpopulation (Estimate) ---------------
for (sup in unique(NHISS_correlation_capZ$Subpopulation)) {
    df <- filter(NHISS_correlation_capZ, Subpopulation == sup)
    sigGenes <- filter(NHISS_correlation_capZ, p.value <= 0.05) 
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of NHISS Correlation Estimates - ", sup, " Monocytes") 
    ggplot(df, aes(Gene_Ordered, Timepoint, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Timepoint", fill = "Coefficient")+
        scale_y_discrete(labels = timepoint_labels)+# Apply custom y-axis labels
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}
for (sup in unique(NHISS_correlation_capZ$Subpopulation)) {
    df <- filter(NHISS_correlation_capZ, Subpopulation == sup)
    # sigGenes <- filter(NHISS_correlation_capZ, p.value <= 0.05) if you want to compare all of them
    sigGenes <- filter(df, p.value <= 0.05)
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of NHISS Correlation Estimates - ", sup, " Monocytes (only sigGenes)") 
    ggplot(df, aes(Gene_Ordered, Timepoint, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Timepoint", fill = "Coefficient")+
        scale_y_discrete(labels = timepoint_labels)+# Apply custom y-axis labels
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}

#### Aging  Split by Timepoint (Estimate) ---------------
for (tp in unique(Age_correlation_capZ$Timepoint)) {
    df <- filter(Age_correlation_capZ, Timepoint == tp)
    sigGenes <- filter(Age_correlation_capZ, p.value <= 0.05)
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of Aging Correlation Estimates - ", tp) 
    ggplot(df, aes(Gene_Ordered, Subpopulation, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Monocyte Suptype", fill = "Coefficient")+
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}

for (tp in unique(Age_correlation_capZ$Timepoint)) {
    df <- filter(Age_correlation_capZ, Timepoint == tp)
    # sigGenes <- filter(Age_correlation_capZ, p.value <= 0.05) #if you want to compare all of them
    sigGenes <- filter(df, p.value <= 0.05)
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of Aging Correlation Estimates - ", tp, " (only sigGenes)") 
    ggplot(df, aes(Gene_Ordered, Subpopulation, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Monocyte Suptype", fill = "Coefficient")+
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}

#### NHISS Split by Timepoint (Estimate) ---------------
for (tp in unique(NHISS_correlation_capZ$Timepoint)) {
    df <- filter(NHISS_correlation_capZ, Timepoint == tp)
    sigGenes <- filter(NHISS_correlation_capZ, p.value <= 0.05) 
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of NHISS Correlation Estimates - ", tp) 
    ggplot(df, aes(Gene_Ordered, Subpopulation, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Monocyte Suptype", fill = "Coefficient")+
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}
for (tp in unique(NHISS_correlation_capZ$Timepoint)) {
    df <- filter(NHISS_correlation_capZ, Timepoint == tp)
    # sigGenes <- filter(NHISS_correlation_capZ, p.value <= 0.05) if you want to compare all of them
    sigGenes <- filter(df, p.value <= 0.05)
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of NHISS Correlation Estimates - ", tp, " (only sigGenes)") 
    ggplot(df, aes(Gene_Ordered, Subpopulation, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Monocyte Suptype", fill = "Coefficient")+
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}

sigGenes <- filter(Age_correlation_capZ, p.value <= 0.05) 
unique(sigGenes$Gene)

sigGenes <- filter(NHISS_correlation_capZ, p.value <= 0.05) 
unique(sigGenes$Gene)

# *Demographics of Ctr & Patients* -----------------------------------------------------------------------------

# for FACS
metadataP_FACS <- metadataP %>% filter(!SampleID %in% Unmatched_TP0_FACS)
Demographics_FACS <- demographics_N_Age_Sex(metadataP_FACS)
write.csv(Demographics_FACS, file = "Demographics_FACS.csv", row.names = FALSE)
# for Gene expr.
metadataP_Gene <- metadataP %>% filter(SampleID %in% unique(data_mean_matched$SampleID))
Demographics_Gene <- demographics_N_Age_Sex(metadataP_Gene)
write.csv(Demographics_Gene, file = "Demographics_Gene.csv", row.names = FALSE)

# *Age & NHISS correlation* -----------------------------------------------------------------------------
folder <- "LinReg_AgeVsNHISS"
createFolder(folder)

Metadata_NHISS <- merge(Metadata_NHISS, metadataP, by = "SampleID")
Metadata_NHISS  <- Metadata_NHISS  %>% filter(!is.na(NHISS), !is.na(Age), !is.na(Timepoint))

AgeVsNHISS_correlation <- data.frame(Group = character(), Timepoint = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

for (tp in unique(Metadata_NHISS$Timepoint)) {
    Metadata_NHISS_tp <- Metadata_NHISS %>% filter(Timepoint == tp)
    
    if (nrow(Metadata_NHISS_tp) <= 0) next
    # Outlier removal
    lm_initial <- lm(NHISS ~ Age, data = Metadata_NHISS_tp)
    residuals_values <- residuals(lm_initial)
    outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
    # Remove outliers only if any are detected
    if (length(outlier_indices) > 0) {
        Metadata_NHISS_tp_filtered <- Metadata_NHISS_tp[-outlier_indices, ]
    } else {
        Metadata_NHISS_tp_filtered <- Metadata_NHISS_tp  # No outliers, keep the original dataframe
    }
    
    if (nrow(Metadata_NHISS_tp_filtered) < 2) next  # Ensure there are enough data points after outlier removal
    
    tryCatch({
        # Pearson's correlation
        x <- cor.test(Metadata_NHISS_tp_filtered$NHISS, Metadata_NHISS_tp_filtered$Age)
        row <- data.frame(
            Group = "All",
            Timepoint = tp,
            N = nrow(Metadata_NHISS_tp_filtered),
            p.value = x$p.value,
            Estimate = x$estimate,
            Est_CI_Lower = x$conf.int[1],
            Est_CI_Upper = x$conf.int[2],
            Coefficient = NA,
            Coef_CI_Lower = NA,
            Coef_CI_Upper = NA  
        )
        
        # Linear regression for coefficient and CI
        lm_result <- lm(NHISS ~ Age, data = Metadata_NHISS_tp_filtered)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        
        row$Coefficient <- coef_x
        row$Coef_CI_Lower <- ci[1]
        row$Coef_CI_Upper <- ci[2]
        
            titel <- paste("Age vs NHISS at ", tp," (All Samples)", sep="")
            file_name <- paste("LinReg_AgeVsNHISS/",titel, ".png", sep = "")
            ggplot(data = Metadata_NHISS_tp_filtered, aes(x = Age, y = NHISS)) +
                geom_smooth(method = "glm", color = "black") +
                geom_point(aes(color = Timepoint), size = 2) +
                ggtitle(titel) +
                xlab("Age") +
                ylab("NHISS") +
                scale_color_manual(values = my_colors) +
                theme(text = element_text(size=14)) +
                theme(
                    panel.background = element_rect(fill = "white", color = "black"),
                    plot.background = element_rect(fill = "white", color = "black"),
                    text = element_text(color = "black"),
                    panel.grid.major = element_line(color = "gray", size = 0.2),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.x = element_blank()
                )
            ggsave(filename=file_name)
            
        AgeVsNHISS_correlation <- rbind(AgeVsNHISS_correlation, row)
    }, error = function(e) {
        cat("Error in Pearson's Correlation for timepoint", tp, "- skipping this comparison.\n")
    })
    for (sex in unique(Metadata_NHISS_tp$Sex)) {
        Metadata_NHISS_sex <- Metadata_NHISS_tp %>% filter(Sex == sex)
        # Outlier removal
        lm_initial <- lm(NHISS ~ Age, data = Metadata_NHISS_sex)
        residuals_values <- residuals(lm_initial)
        outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
        # Remove outliers only if any are detected
        if (length(outlier_indices) > 0) {
            Metadata_NHISS_sex_filtered <- Metadata_NHISS_sex[-outlier_indices, ]
        } else {
            Metadata_NHISS_sex_filtered <- Metadata_NHISS_sex  # No outliers, keep the original dataframe
        }
        
        if (nrow(Metadata_NHISS_sex_filtered) < 2) next  # Ensure there are enough data points after outlier removal
        
        
        tryCatch({
            # Pearson's correlation for Sex
            x <- cor.test(Metadata_NHISS_sex_filtered$NHISS, Metadata_NHISS_sex_filtered$Age)
            row <- data.frame(
                Group = sex,
                Timepoint = tp,
                N = nrow(Metadata_NHISS_sex_filtered),
                p.value = x$p.value,
                Estimate = x$estimate,
                Est_CI_Lower = x$conf.int[1],
                Est_CI_Upper = x$conf.int[2],
                Coefficient = NA,
                Coef_CI_Lower = NA,
                Coef_CI_Upper = NA  
            )
            
            # Linear regression for coefficient and CI
            lm_result_sex <- lm(NHISS ~ Age, data = Metadata_NHISS_sex_filtered)
            coef_x <- summary(lm_result_sex)$coefficients[2, 1]
            ci <- confint(lm_result_sex)[2, ]
            
            row$Coefficient <- coef_x
            row$Coef_CI_Lower <- ci[1]
            row$Coef_CI_Upper <- ci[2]
            
            titel <- paste("Age vs NHISS at ", tp," (",sex, " Samples)", sep="")
            file_name <- paste("LinReg_AgeVsNHISS/",titel, ".png", sep = "")
            ggplot(data = Metadata_NHISS_sex_filtered, aes(x = Age, y = NHISS)) +
                    geom_smooth(method = "glm", color = "black") +
                    geom_point(aes(color = Sex), size = 2) +
                    ggtitle(titel) +
                    xlab("Age") +
                    ylab("NHISS") +
                    scale_color_manual(values = c("M" = "#A67C00", "F" = "#1D04C2")) +
                    theme(text = element_text(size=14)) +
                    theme(
                        panel.background = element_rect(fill = "white", color = "black"),
                        plot.background = element_rect(fill = "white", color = "black"),
                        text = element_text(color = "black"),
                        panel.grid.major = element_line(color = "gray", size = 0.2),
                        panel.grid.minor = element_blank(),
                        panel.grid.major.x = element_blank()
                    )
                ggsave(filename=file_name)
                AgeVsNHISS_correlation <- rbind(AgeVsNHISS_correlation, row)
        }, error = function(e) {
            cat("Error for", tp, "&", sex, "- skipping this comparison.\n")
        }) 
    }
    for (cat in unique(Metadata_NHISS_tp$Category)) {
        Metadata_NHISS_cat <- Metadata_NHISS_tp %>% filter(Category == cat)
        # Outlier removal
        lm_initial <- lm(NHISS ~ Age, data = Metadata_NHISS_cat)
        residuals_values <- residuals(lm_initial)
        outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
        # Remove outliers only if any are detected
        if (length(outlier_indices) > 0) {
            Metadata_NHISS_cat_filtered <- Metadata_NHISS_cat[-outlier_indices, ]
        } else {
            Metadata_NHISS_cat_filtered <- Metadata_NHISS_cat  # No outliers, keep the original dataframe
        }
        
        if (nrow(Metadata_NHISS_cat_filtered) < 2) next  # Ensure there are enough data points after outlier removal
        
        tryCatch({
            # Pearson's correlation for Category
            x <- cor.test(Metadata_NHISS_cat_filtered$NHISS, Metadata_NHISS_cat_filtered$Age)
            row <- data.frame(
                Group = cat,
                Timepoint = tp,
                N = nrow(Metadata_NHISS_cat_filtered),
                p.value = x$p.value,
                Estimate = x$estimate,
                Est_CI_Lower = x$conf.int[1],
                Est_CI_Upper = x$conf.int[2],
                Coefficient = NA,
                Coef_CI_Lower = NA,
                Coef_CI_Upper = NA  
            )
            
            # Linear regression for coefficient and CI
            lm_result_cat <- lm(NHISS ~ Age, data = Metadata_NHISS_cat_filtered)
            coef_x <- summary(lm_result_cat)$coefficients[2, 1]
            ci <- confint(lm_result_cat)[2, ]
            
            row$Coefficient <- coef_x
            row$Coef_CI_Lower <- ci[1]
            row$Coef_CI_Upper <- ci[2]
            
            titel <- paste("Age vs NHISS at ", tp, " (", cat, " Samples)", sep="")
            file_name <- paste("LinReg_AgeVsNHISS/", titel, ".png", sep="")
            ggplot(data = Metadata_NHISS_cat_filtered, aes(x = Age, y = NHISS)) +
                geom_smooth(method = "glm", color = "black") +
                geom_point(aes(color = Category), size = 2) +
                ggtitle(titel) +
                xlab("Age") +
                ylab("NHISS") +
                scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606")) +
                theme(text = element_text(size=14)) +
                theme(
                    panel.background = element_rect(fill = "white", color = "black"),
                    plot.background = element_rect(fill = "white", color = "black"),
                    text = element_text(color = "black"),
                    panel.grid.major = element_line(color = "gray", size = 0.2),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.x = element_blank()
                )
            ggsave(filename = file_name)
            AgeVsNHISS_correlation <- rbind(AgeVsNHISS_correlation, row)
        }, error = function(e) {
            cat("Error for", tp, "&", cat, "- skipping this comparison.\n")
        }) 
    }
}
    
file_name_s<- "LinReg_AgeVsNHISS/Pearson's Correlation_AgeVsNHISS.csv"
write.csv(AgeVsNHISS_correlation, file_name_s, row.names = FALSE)

# *SPAN Regression --------------------------
folder <- "LinReg_SPAN"
createFolder(folder)

# Create the SPAN variable as per Almekhlafi et. al 2014:
dataFINALmean <- dataFINALmean %>%
    mutate(SPAN = Age + NHISS)
# Initialize dataframes to store correlation results
SPAN_correlation_capZ <- data.frame(Subpopulation = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

SPAN_correlation_Recovery_capZ <- data.frame(Subpopulation = character(), Recovery = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

SPAN_correlation_Category_capZ <- data.frame(Subpopulation = character(), Category = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

for (gene in unique(dataFINALmean$Gene)) {
    df <- dataFINALmean %>% filter(Gene == gene)
    df <- df %>% filter(!is.na(Timepoint), !is.na(mean_capZ), !is.na(Subpopulation), !is.na(Recovery), !is.na(Category))
    
    for (subpop in unique(df$Subpopulation)) {
        df_subpop <- df %>% filter(Subpopulation == subpop)
        
        for (tp in unique(df_subpop$Timepoint)) {
            df_tp <- df_subpop %>% filter(Timepoint == tp, !is.na(SPAN))
            if (nrow(df_tp) <= 0) next
            
            # Outlier removal
            lm_initial <- lm(mean_capZ ~ SPAN, data = df_tp)
            residuals_values <- residuals(lm_initial)
            outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
            # Remove outliers only if any are detected
            if (length(outlier_indices) > 0) {
                df_tp_filtered <- df_tp[-outlier_indices, ]
            } else {
                df_tp_filtered <- df_tp  # No outliers, keep the original dataframe
            }
            
            if (nrow(df_tp_filtered) < 2) next  # Ensure there are enough data points after outlier removal
            
            tryCatch({
                # Pearson's correlation
                x <- cor.test(df_tp_filtered$mean_capZ, df_tp_filtered$SPAN)
                row <- data.frame(
                    Subpopulation = subpop,
                    Timepoint = tp,
                    Gene = gene,
                    GeneID = df_tp_filtered$GeneID[1],
                    Gene_Group = df_tp_filtered$Gene_Group[1],
                    N = nrow(df_tp_filtered),
                    p.value = x$p.value,
                    Estimate = x$estimate,
                    Est_CI_Lower = x$conf.int[1],
                    Est_CI_Upper = x$conf.int[2],
                    Coefficient = NA,
                    Coef_CI_Lower = NA,
                    Coef_CI_Upper = NA  
                )
                
                # Linear regression for coefficient and CI
                lm_result <- lm(mean_capZ ~ SPAN, data = df_tp_filtered)
                coef_x <- summary(lm_result)$coefficients[2, 1]
                ci <- confint(lm_result)[2, ]
                
                row$Coefficient <- coef_x
                row$Coef_CI_Lower <- ci[1]
                row$Coef_CI_Upper <- ci[2]
                
                if (row$p.value<=0.05){
                    titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," (All Samples)", sep="")
                    file_name <- paste("LinReg_Plots_SPAN_capZ/",titel, ".png", sep = "")
                    ggplot(data = df_tp_filtered, aes(x = SPAN, y = mean_capZ)) +
                        geom_smooth(method = "glm", color = "black") +
                        geom_point(aes(color = Timepoint), size = 2) +
                        ggtitle(titel) +
                        xlab("SPAN") +
                        ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                        scale_color_manual(values = my_colors) +
                        theme(text = element_text(size=14)) +
                        theme(
                            panel.background = element_rect(fill = "white", color = "black"),
                            plot.background = element_rect(fill = "white", color = "black"),
                            text = element_text(color = "black"),
                            panel.grid.major = element_line(color = "gray", size = 0.2),
                            panel.grid.minor = element_blank(),
                            panel.grid.major.x = element_blank()
                        )
                    ggsave(filename=file_name)
                }
                SPAN_correlation_capZ <- rbind(SPAN_correlation_capZ, row)
            }, error = function(e) {
                cat("Error in Pearson's Correlation for gene", gene, "& Subpopulation with SPAN", subpop, "- skipping this comparison.\n")
            })
            
            for (recovery in unique(df_tp$Recovery)) {
                df_recovery <- df_tp %>% filter(Recovery == recovery, !is.na(SPAN))
                # Outlier removal
                lm_initial <- lm(mean_capZ ~ SPAN, data = df_recovery)
                residuals_values <- residuals(lm_initial)
                outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
                # Remove outliers only if any are detected
                if (length(outlier_indices) > 0) {
                    df_recovery_filtered <- df_recovery[-outlier_indices, ]
                } else {
                    df_recovery_filtered <- df_recovery  # No outliers, keep the original dataframe
                }
                
                if (nrow(df_recovery_filtered) < 2) next  # Ensure there are enough data points after outlier removal
                
                
                tryCatch({
                    # Pearson's correlation for recovery
                    x <- cor.test(df_recovery_filtered$mean_capZ, df_recovery_filtered$SPAN)
                    row <- data.frame(
                        Subpopulation = subpop,
                        Recovery = recovery,
                        Timepoint = tp,
                        Gene = gene,
                        GeneID = df_recovery_filtered$GeneID[1],
                        Gene_Group = df_recovery_filtered$Gene_Group[1],
                        N = nrow(df_recovery_filtered),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA,
                        Coef_CI_Lower = NA,
                        Coef_CI_Upper = NA  
                    )
                    
                    # Linear regression for coefficient and CI
                    lm_result_recovery <- lm(mean_capZ ~ SPAN, data = df_recovery_filtered)
                    coef_x <- summary(lm_result_recovery)$coefficients[2, 1]
                    ci <- confint(lm_result_recovery)[2, ]
                    
                    row$Coefficient <- coef_x
                    row$Coef_CI_Lower <- ci[1]
                    row$Coef_CI_Upper <- ci[2]
                    
                    if (row$p.value<=0.05){
                        titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," ( ", recovery, " Samples)", sep="")
                        file_name <- paste("LinReg_Plots_SPAN_capZ/",titel, ".png", sep = "")
                        ggplot(data = df_recovery_filtered, aes(x = SPAN, y = mean_capZ)) +
                            geom_smooth(method = "glm", color = "black") +
                            geom_point(aes(color = recovery), size = 2) +
                            ggtitle(titel) +
                            xlab("SPAN") +
                            ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                            scale_color_manual(values = c("Good" = "darkgreen", "Bad" = "#700606")) +
                            theme(text = element_text(size=14)) +
                            theme(
                                panel.background = element_rect(fill = "white", color = "black"),
                                plot.background = element_rect(fill = "white", color = "black"),
                                text = element_text(color = "black"),
                                panel.grid.major = element_line(color = "gray", size = 0.2),
                                panel.grid.minor = element_blank(),
                                panel.grid.major.x = element_blank()
                            )
                        ggsave(filename=file_name)
                    }
                    SPAN_correlation_Recovery_capZ <- rbind(SPAN_correlation_Recovery_capZ, row)
                }, error = function(e) {
                    cat("Error for", gene, tp, recovery, "&", subpop, "- skipping this comparison.\n")
                }) 
            }
            
            for (cat in unique(df_tp$Category)) {
                df_cat <- df_tp %>% filter(Category == cat, !is.na(SPAN))
                # Outlier removal
                lm_initial <- lm(mean_capZ ~ SPAN, data = df_cat)
                residuals_values <- residuals(lm_initial)
                outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
                # Remove outliers only if any are detected
                if (length(outlier_indices) > 0) {
                    df_cat_filtered <- df_cat[-outlier_indices, ]
                } else {
                    df_cat_filtered <- df_cat  # No outliers, keep the original dataframe
                }
                
                if (nrow(df_cat_filtered) < 2) next  # Ensure there are enough data points after outlier removal
                
                tryCatch({
                    # Pearson's correlation for Category
                    x <- cor.test(df_cat_filtered$mean_capZ, df_cat_filtered$SPAN)
                    row <- data.frame(
                        Subpopulation = subpop,
                        Category = cat,
                        Timepoint = tp,
                        Gene = gene,
                        GeneID = df_cat_filtered$GeneID[1],
                        Gene_Group = df_cat_filtered$Gene_Group[1],
                        N = nrow(df_cat_filtered),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA,
                        Coef_CI_Lower = NA,
                        Coef_CI_Upper = NA  
                    )
                    
                    # Linear regression for coefficient and CI
                    lm_result_cat <- lm(mean_capZ ~ SPAN, data = df_cat_filtered)
                    coef_x <- summary(lm_result_cat)$coefficients[2, 1]
                    ci <- confint(lm_result_cat)[2, ]
                    
                    row$Coefficient <- coef_x
                    row$Coef_CI_Lower <- ci[1]
                    row$Coef_CI_Upper <- ci[2]
                    
                    if (row$p.value<=0.05){
                        titel <- paste(gene, " expression in ", subpop, " monocytes at ", tp," ( ", cat, " Samples)", sep="")
                        file_name <- paste("LinReg_Plots_SPAN_capZ/",titel, ".png", sep = "")
                        ggplot(data = df_cat_filtered, aes(x = SPAN, y = mean_capZ)) +
                            geom_smooth(method = "glm", color = "black") +
                            geom_point(aes(color = Category), size = 2) +
                            ggtitle(titel) +
                            xlab("SPAN") +
                            ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                            scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606")) +
                            theme(text = element_text(size=14)) +
                            theme(
                                panel.background = element_rect(fill = "white", color = "black"),
                                plot.background = element_rect(fill = "white", color = "black"),
                                text = element_text(color = "black"),
                                panel.grid.major = element_line(color = "gray", size = 0.2),
                                panel.grid.minor = element_blank(),
                                panel.grid.major.x = element_blank()
                            )
                        ggsave(filename=file_name)
                    }
                    SPAN_correlation_Category_capZ <- rbind(SPAN_correlation_Category_capZ, row)
                }, error = function(e) {
                    cat("Error for", gene, tp, cat, "&", subpop, "- skipping this comparison.\n")
                }) 
            }
        }
    }
}

file_name <- "LinReg_Results_capZ/Pearson's Correlation_SPAN_capZ.csv"
write.csv(SPAN_correlation_capZ, file_name, row.names = FALSE)

file_name_s <- "LinReg_Results_capZ/Pearson's Correlation_SPAN_Recovery_capZ.csv"
write.csv(SPAN_correlation_Recovery_capZ, file_name_s, row.names = FALSE)

file_name_s <- "LinReg_Results_capZ/Pearson's Correlation_SPAN_Category_capZ.csv"
write.csv(SPAN_correlation_Category_capZ, file_name_s, row.names = FALSE)

#### Extract significant LinRegs --------------------------------------------------------

SPAN_correlation_capZ$Recovery <- "All"
SPAN_correlation_capZ$Category <- "All"
SPAN_correlation_Recovery_capZ$Category <- "All"
SPAN_correlation_Category_capZ$Recovery<- "All"

All_SPAN_Pearson_capZ <- rbind(SPAN_correlation_capZ,SPAN_correlation_Recovery_capZ, SPAN_correlation_Category_capZ)
SPAN_sigGenes_Pearson <- All_SPAN_Pearson_capZ %>% filter(All_SPAN_Pearson_capZ$p.value <0.05)
write.csv(SPAN_sigGenes_Pearson, "LinReg_Results_capZ/SPAN_sigGenes_Pearson.csv", row.names = FALSE)

# *Sig Regression - Gene Overview* --------------------------
# with a mean for all Timepoints
Gene_Overview_Age <- Age_correlation_capZ %>%
    filter(p.value <= 0.05) %>% 
    group_by(Subpopulation, Gene)%>%
    summarise(Age_Correlation = mean(Estimate, na.rm = TRUE))

Gene_Overview_NHISS <- NHISS_correlation_capZ %>%
    filter(p.value <= 0.05) %>% 
    group_by(Subpopulation, Gene)%>%
    summarise(NHISS_Correlation = mean(Estimate, na.rm = TRUE))

Gene_Overview_SPAN <- SPAN_correlation_capZ %>%
    filter(p.value <= 0.05) %>% 
    group_by(Subpopulation, Gene)%>%
    summarise(SPAN_Correlation = mean(Estimate, na.rm = TRUE))

Gene_Overview_woTP <- merge(Gene_Overview_Age,
                       Gene_Overview_NHISS, 
                       by = c("Gene", "Subpopulation"), 
                       all = TRUE)

Gene_Overview_woTP <- merge(Gene_Overview_woTP,
                       Gene_Overview_SPAN, 
                       by = c("Gene", "Subpopulation"), 
                       all = TRUE)

# table including the details for the timepoints
Gene_Overview_NHISS <- NHISS_correlation_capZ %>%
    filter(p.value <= 0.05) %>% 
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarise(NHISS_Correlation = Estimate)

Gene_Overview_SPAN <- SPAN_correlation_capZ %>%
    filter(p.value <= 0.05) %>% 
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarise(SPAN_Correlation = Estimate)

Gene_Overview_Age <- Age_correlation_capZ %>%
    filter(p.value <= 0.05) %>% 
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarise(Age_Correlation = Estimate)

Gene_Overview <- merge(Gene_Overview_Age,
                       Gene_Overview_NHISS, 
                       by = c("Gene", "Subpopulation", "Timepoint"), 
                       all = TRUE)

Gene_Overview <- merge(Gene_Overview,
                       Gene_Overview_SPAN, 
                       by = c("Gene", "Subpopulation", "Timepoint"), 
                       all = TRUE)
    
# Replace the Estimated with the arrows!    
Gene_Overview_woTP <- Gene_Overview_woTP %>%
    mutate(across(c(Age_Correlation, NHISS_Correlation, SPAN_Correlation), 
                  ~ ifelse(is.na(.), "-", ifelse(. > 0, "↑", "↓"))))
Gene_Overview <- Gene_Overview %>%
    mutate(across(c(Age_Correlation, NHISS_Correlation, SPAN_Correlation), 
                  ~ ifelse(is.na(.), "-", ifelse(. > 0, "↑", "↓"))))

# Save the files:
file_name_s <- "LinReg_Results_capZ/Gene_Overview_woTP.csv"
write.csv(Gene_Overview_woTP, file_name_s, row.names = FALSE)
file_name_s <- "LinReg_Results_capZ/Gene_Overview.csv"
write.csv(Gene_Overview, file_name_s, row.names = FALSE)
