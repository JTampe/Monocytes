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
plot_linear_regression1 <- function(data, cells_colname, group_name, output_folder) {
    
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
        # Save the plot with all timepoints
        ggsave(filename = file.path(output_folder, paste0("LinReg_", group_name, "_", cells_colname, ".png")), plot = p_group)
        
        p_individual <- ggplot(data_filtered %>% filter(Timepoint == timepoint), 
                               aes(x = Age, y = .data[[cells_colname]], color = Timepoint)) +
            geom_point(size = 2) +
            geom_smooth(method = "lm", aes(color = Timepoint), se = TRUE, size = 1) +
            ggtitle(paste("Linear Regression -", group_name, "(", cells_colname, ")", " - Timepoint:", timepoint)) +
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
            scale_color_manual(values = my_colors[timepoint]) +
            scale_fill_manual(values = my_colors[timepoint])
        
        # Save the individual plot
        ggsave(filename = file.path(output_folder, paste0("LinReg_", group_name, "_", cells_colname, "_Timepoint_", timepoint, ".png")), plot = p_individual,
               width = 7, height = 6)
        
    }
    
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
        # Save the plot with all timepoints
        ggsave(filename = file.path(output_folder, paste0("LinReg_", group_name, "_", cells_colname, ".png")), plot = p_group)
        
        p_individual <- ggplot(data_filtered %>% filter(Timepoint == timepoint), 
                               aes(x = Age, y = .data[[cells_colname]], color = Timepoint)) +
            geom_point(size = 2) +
            geom_smooth(method = "lm", aes(color = Timepoint), se = TRUE, size = 1) +
            ggtitle(paste("Linear Regression -", group_name, "(", cells_colname, ")", " - Timepoint:", timepoint)) +
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
            scale_color_manual(values = my_colors[timepoint]) +
            scale_fill_manual(values = my_colors[timepoint])
        
        # Save the individual plot
        ggsave(filename = file.path(output_folder, paste0("LinReg_", group_name, "_", cells_colname, "_Timepoint_", timepoint, ".png")), plot = p_individual,
               width = 8, height = 5)
        
    }
    
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


plot_linear_regression_NHISS <- function(data, cells_colname, group_name, output_folder, result, x_axis_var) {
    
    # Fit the initial linear model to identify outliers
    lm_initial <- lm(data[[cells_colname]] ~ data[[x_axis_var]])
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
    
    result <- rbind(result, rows_group)  # Update global variable
    
    # Plot regression for the specific group with all timepoints
    p_group <- ggplot(data_filtered, aes_string(x = x_axis_var, y = cells_colname, color = "Timepoint")) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", aes(color = Timepoint), se = FALSE, size = 1) +
        ggtitle(paste("Linear Regression -", group_name, "(", cells_colname, ")", sep = " ")) +
        xlab(x_axis_var) +
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
    ggsave(filename = file.path(output_folder, paste0("LinReg_", x_axis_var, "_", group_name, "_", cells_colname, ".png")), plot = p_group)
    
    # Create an additional plot with only significant timepoints if there are any
    significant_timepoints <- rows_group$Timepoint[rows_group$p.value <= 0.05]
    
    if (length(significant_timepoints) > 0) {
        # Plot with all significant timepoints together
        p_significant <- ggplot(data_filtered %>% filter(Timepoint %in% significant_timepoints), 
                                aes_string(x = x_axis_var, y = cells_colname, color = "Timepoint")) +
            geom_point(size = 2) +
            geom_smooth(method = "lm", aes(color = Timepoint), se = TRUE, size = 1) +
            ggtitle(paste("Linear Regression - Significant Only -", group_name, "(", cells_colname, ")", sep = " ")) +
            xlab(x_axis_var) +
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
        ggsave(filename = file.path(output_folder, paste0("LinReg_SignificantOnly_", x_axis_var, "_", group_name, "_", cells_colname, ".png")), plot = p_significant)
    }
    return(result)
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
            my_colors <- c("darkgrey", "black", "black", "black", "black")
            
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
                                   label.y = max(df_test[[var_name]]),
                                   size = 8) +
                scale_x_discrete(labels = c(
                    "TP0" = "Control",
                    "TP1" = "24 hours",
                    "TP2" = "3-5 days",
                    "TP3" = "1 month",
                    "TP4" = "3 months"
                )) +
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
                           "Bad" = "#700606","Good" = "#A67C00")
            
            # Filter my_colors to match the levels of color_var in df_test
            used_colors <- my_colors[unique(df_test[[color_var]])]
            
            # Boxplot with Correct Legend
            box_plot <- ggboxplot(df_test, x = timepoint_col, y = var_name, color = color_var, 
                                  add = "jitter", legend = "right", ylab = paste(var_name, "percentage")) +
                # geom_hline(yintercept = median_TP0, linetype = 2) +
                scale_x_discrete(labels = c("TP0" = "Control", "TP1" = "24 hours", 
                                            "TP2" = "3-5 days", "TP3" = "1 month", 
                                            "TP4" = "3 months")) +
                scale_color_manual(values = used_colors, 
                                   name = "Group",
                                   labels = names(used_colors))
            
            
            # Annotate significant comparisons for each category
            for (color in unique(df_test[[color_var]])) {
                df_color <- results[results$Category == color & results$Variable == var_name, ]
                y_base <- max(df_test[[var_name]], na.rm = TRUE) + 0.1
                y_step <- 0.15
                for (tp in c("TP1", "TP2", "TP3", "TP4")) {
                    p_value <- df_color[[paste(tp, "P_Value", sep = "_")]]
                    if (!is.na(p_value) && p_value <= 0.05) {
                        asterisk <- ifelse(p_value <= 0.001, "***",
                                           ifelse(p_value <= 0.01, "**", "*"))
                        x_position <- which(levels(df_test[[timepoint_col]]) == tp)
                        box_plot <- box_plot +
                            annotate("text", x = x_position, y = y_base, 
                                     label = asterisk, size = 12, color = my_colors[color]) #+
                        #annotate("segment", x = x_position - 0.15, xend = x_position + 0.15, 
                        #         y = y_base - 0.02, yend = y_base - 0.02, color = my_colors[color])
                        y_base <- y_base + y_step
                    }
                }
            }
            
            ggsave(filename = file.path(plots_folder, paste(var_name, "_Boxplot_Wilcox.png", sep = "")), 
                   plot = box_plot, width = 7.435, height = 5)
            
            # # Barplot with SEM
            # mean_values <- df_test %>%
            #     group_by(get(timepoint_col), get(color_var)) %>%
            #     summarise(mean_value = mean(get(var_name), na.rm = TRUE),
            #               sd_value = sd(get(var_name), na.rm = TRUE),
            #               n = n(), .groups = "drop") %>%
            #     mutate(SEM = sd_value / sqrt(n)) %>%
            #     rename(Timepoint = `get(timepoint_col)`, Group = `get(color_var)`)
            # 
            # dodge_position <- position_dodge(width = 0.8)
            # bar_plot <- ggbarplot(mean_values, x = "Timepoint", y = "mean_value", fill = "Group", 
            #                       ylab = paste(var_name, "mean (+-SEM) percentage"), 
            #                       add = "mean_se", width = 0.6, position = dodge_position) +
            #     scale_fill_manual(values = my_colors) +
            #     scale_x_discrete(labels = c("TP0" = "Control", "TP1" = "24 hours", 
            #                                 "TP2" = "3-5 days", "TP3" = "1 month", 
            #                                 "TP4" = "3 months")) +
            #     geom_errorbar(aes(ymin = mean_value - SEM, ymax = mean_value + SEM), 
            #                   width = 0.2, position = dodge_position) +
            #     theme_minimal()
            # 
            # ggsave(filename = file.path(plots_folder, paste(var_name, "_Barplot_Wilcox_SEM.png", sep = "")), plot = bar_plot)
        }
    }
    
    # Save results to a CSV file
    write.csv(results, file.path(results_folder, paste(results_name, "_results.csv", sep = "")), row.names = FALSE)
    
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

# Helper function to remove outliers
remove_outliers <- function(data, variable) {
    lm_initial <- lm(mean_capZ ~ get(variable), data = data)
    residuals_values <- residuals(lm_initial)
    outlier_indices <- which(abs(residuals_values) > (3 * sd(residuals_values)))
    if (length(outlier_indices) > 0) {
        data <- data[-outlier_indices, ]
    }
    return(data)
}

# Helper function to append correlation and regression results
append_correlation_results <- function(data, variable, gene, subpop, tp, subgroup, value, results_df, plot_folder) {
    tryCatch({
        if (nrow(data) < 2) {
            cat("Not enough data points for linear regression (less than 2) for", gene, subpop, "Monocytes at", tp, "for", subgroup, value, "- skipping this subset.\n")
            return(results_df)
        }
        
        # Pearson correlation
        cor_result <- cor.test(data$mean_capZ, data[[variable]])
        
        # Linear regression
        lm_result <- lm(mean_capZ ~ get(variable), data = data)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        IC <- summary(lm_result)$coefficients[1, 1]
        
        # Append to results dataframe
        row <- data.frame(
            Gene = gene,
            GeneID = data$GeneID[1],
            Gene_Group = data$Gene_Group[1],
            Subpopulation = subpop,
            Timepoint = tp,
            N = nrow(data),
            p.value = cor_result$p.value,
            Estimate = cor_result$estimate,
            Est_CI_Lower = cor_result$conf.int[1],
            Est_CI_Upper = cor_result$conf.int[2],
            Coefficient = coef_x,
            Coef_CI_Lower = ci[1],
            Coef_CI_Upper = ci[2],
            Intercept = IC,
            Subgroup = subgroup,
            Subgroup_Value = value
        )
        
        results_df <- rbind(results_df, row)
        
        if (row$p.value <= 0.05) {
            #titel <- paste(gene, " expression in ", subpop, " Monocytes at ", tp, " ( ", value, " Samples)", sep="")
            titel <- paste(value, " Samples -", gene," expr. of ", subpop, " Monocytes at ", tp, sep="")
            file_name <- paste0(plot_folder, "/",
                                value, "_Samples_",
                                gene, "_",
                                subpop, "Monocytes_", tp,
                                ".png", sep="")
            # Create the plot and save it
            plot <- ggplot(data = data, aes_string(x = variable, y = "mean_capZ")) +
                geom_smooth(method = "glm", color = "black") +
                geom_point(aes(color = Recovery), size = 2) +
                ggtitle(titel) +
                xlab(variable) +
                ylab(paste(gene, "expression (Z-scored from CT - CT Sample Median)")) +
                scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606",
                                              "Male" = "#A67C00", "Female" = "#1D04C2",
                                              "Good" = "darkgreen", "Bad" = "#700606")) +
                theme(text = element_text(size = 14)) +
                theme(
                    panel.background = element_rect(fill = "white", color = "black"),
                    plot.background = element_rect(fill = "white", color = "black"),
                    text = element_text(color = "black"),
                    panel.grid.major = element_line(color = "gray", size = 0.2),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.x = element_blank()
                )
            
            # Save plot with specified dimensions
            ggsave(filename = file_name, plot = plot, width = 8, height = 6, dpi = 300)
            cat("Plot saved successfully to:", file_name, "\n")
        } else {}
        
        return(results_df)
    }, error = function(e) {
        cat("Error for", gene, subpop, "Monocytes at", tp, "for", subgroup, value, "- skipping this subset.\n")
        return(results_df)
    })
}

perform_linear_regression_correlation <- function(data, variable, folder_prefix) {
    # Create the output folders if they do not exist
    plot_folder <- paste0("LinReg_Plots_", folder_prefix, "_capZ")
    result_folder <- paste0("LinReg_Results_capZ")
    if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)
    if (!dir.exists(result_folder)) dir.create(result_folder, recursive = TRUE)
    
    # Initialize dataframes for correlation results
    correlation_results <- data.frame(Subpopulation = character(), 
                                      Timepoint = character(),
                                      Gene = character(), 
                                      GeneID = numeric(), 
                                      Gene_Group = character(), 
                                      N = numeric(), 
                                      p.value = numeric(), 
                                      Estimate = numeric(), 
                                      Est_CI_Lower = numeric(), 
                                      Est_CI_Upper = numeric(), 
                                      Coefficient = numeric(),
                                      Coef_CI_Lower = numeric(),
                                      Coef_CI_Upper = numeric(),
                                      Intercept = numeric(),
                                      Subgroup = character(),
                                      Subgroup_Value = character())
    
    # Loop over genes and combinations of subpopulations, timepoints, etc.
    for (gene in unique(data$Gene)) {
        df <- data %>% filter(Gene == gene)
        df <- df %>% filter(!is.na(Timepoint), !is.na(mean_capZ), !is.na(get(variable))) # Clean data
        
        for (subpop in unique(df$Subpopulation)) {
            df_subpop <- df %>% filter(Subpopulation == subpop)
            for (tp in unique(df_subpop$Timepoint)) {
                df_tp <- df_subpop %>% filter(Timepoint == tp)
                if (nrow(df_tp) <= 0) next
                
                # Outlier removal for the overall group
                df_tp_filtered <- remove_outliers(df_tp, variable)
                if (nrow(df_tp_filtered) == 0) next
                
                # Perform analysis for the whole group
                correlation_results <- append_correlation_results(df_tp_filtered, variable, 
                                                                  gene, subpop, tp, "All", "All", correlation_results, plot_folder)
                
                # # Subgroup analysis for Sex, Category, and Recovery
                # for (subgroup in c("Sex", "Category", "Recovery")) {
                #     for (value in unique(df_tp[[subgroup]])) {
                #         df_subgroup <- df_tp %>% filter(.data[[subgroup]] == value)
                #         df_subgroup_filtered <- remove_outliers(df_subgroup, variable)
                #         if (nrow(df_subgroup_filtered) == 0) next
                #         
                #         correlation_results <- append_correlation_results(df_subgroup_filtered, variable, 
                #                                                           gene, subpop, tp, subgroup, value, correlation_results, plot_folder)
                #    }
                #}
            }
        }
    }
    # Write final results to CSV
    write.csv(correlation_results, file = paste0(result_folder, "/LinReg_", variable, "_capZ.csv"), row.names = FALSE)
    return(correlation_results)
}
