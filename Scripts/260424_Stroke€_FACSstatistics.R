# *FACS Statistics* ----------------
FACSdata_copy <- FACSdata 
#### 1WAY ANOVA* ----------------------------------------------------------------------------------------------------

ANOVA_FACSdata <- FACSdata %>% filter(!(Timepoint == "TP0" & SampleID %in% Unmatched_TP0_FACS))
write_xlsx(ANOVA_FACSdata, "ANOVA_FACSdata.xlsx")

# with wilcox
ANOVA_FACSall <- automate_anova_extraction(output_location, "Wilcox_FACS_Plots", 
                                           "Wilcox_FACS", "y", ANOVA_FACSdata, 
                                           colnames(ANOVA_FACSdata)[9:12], 
                                           "Timepoint")

ANOVA_FACSdata <- ANOVA_FACSdata[order(ANOVA_FACSdata$Sex == "Female", decreasing = TRUE), ]
# ANOVA_FACSsex <- automate_anova_extraction_Category(output_location, "Wilcox_FACS_Plots_Sex", 
#                                             "Wilcox_FACS_Sex", "y", ANOVA_FACSdata, 
#                                             colnames(ANOVA_FACSdata)[9:12], 
#                                             c(colnames(ANOVA_FACSdata)[2]),c(colnames(ANOVA_FACSdata)[22]))
# 
# ANOVA_FACSdata <- ANOVA_FACSdata[order(ANOVA_FACSdata$Category == "MODERATE", decreasing = FALSE), ]
# ANOVA_FACScategory <- automate_anova_extraction_Category(output_location, "Wilcox_FACS_Plots_Category", 
#                                                     "Wilcox_FACS_Category", "y", ANOVA_FACSdata, 
#                                                     colnames(ANOVA_FACSdata)[9:12], 
#                                                     c(colnames(ANOVA_FACSdata)[2]),c(colnames(ANOVA_FACSdata)[24]))
# 
# ANOVA_FACSdata <- ANOVA_FACSdata[order(ANOVA_FACSdata$Recovery == "Bad", decreasing = TRUE), ]
# ANOVA_FACSrecovery <- automate_anova_extraction_Category(output_location, "Wilcox_FACS_Plots_Recovery", 
#                                                     "Wilcox_FACS_Recovery", "y", ANOVA_FACSdata, 
#                                                     colnames(ANOVA_FACSdata)[9:12], 
#                                                     c(colnames(ANOVA_FACSdata)[2]),c(colnames(ANOVA_FACSdata)[25]))
# write_xlsx(ANOVA_FACSsex, "ANOVA_FACSsex.xlsx")
# write_xlsx(ANOVA_FACScategory, "ANOVA_FACScategory.xlsx")
# write_xlsx(ANOVA_FACSrecovery, "ANOVA_FACSrecovery.xlsx")
# 
# ANOVA_FACSall$Category <- "All"
# ANOVA_FACScombined <- rbind (ANOVA_FACSall,ANOVA_FACSsex)
# ANOVA_FACScombined <- rbind (ANOVA_FACScombined,ANOVA_FACSrecovery)
# ANOVA_FACScombined <- rbind (ANOVA_FACScombined,ANOVA_FACScategory)
# write_xlsx(ANOVA_FACScombined, "ANOVA_FACScombined.xlsx")

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
plot_save <- "n"

# Initialize an empty dataframe to store Wilcoxon test results
results_FACS_Wilcox <- data.frame(Cells = character(), Focus = character(), 
                                  TP0_P_Value = numeric(), TP1_P_Value = numeric(), 
                                  TP2_P_Value = numeric(), TP3_P_Value = numeric(), 
                                  TP4_P_Value = numeric(), stringsAsFactors = FALSE)

for (i in 9:12) {
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
        scale_x_discrete(labels = c(
            "TP0" = "Control",
            "TP1" = "24 hours",
            "TP2" = "3-5 days",
            "TP3" = "1 month",
            "TP4" = "3 months"
        )) +
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
plot_save <- "n"

for (i in 9:12) {
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
        scale_x_discrete(labels = c(
            "TP0" = "Control",
            "TP1" = "24 hours",
            "TP2" = "3-5 days",
            "TP3" = "1 month",
            "TP4" = "3 months"
        )) +
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
plot_save <- "n"

for (i in 9:12) {
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
    
    max_y_value <- mean(df_subpop_clean$data_column, na.rm = TRUE)+5
    
    p <- ggline(df_subpop_clean, 
                x = "Timepoint", 
                y = "data_column", 
                color = "Recovery", 
                add = "mean_se",  # Add mean and standard error
                size = 1.2) +  
        labs(y = paste(cells_colname, "percentage"), 
             title = paste("Mean ± SEM of", cells_colname, "by Recovery across Timepoints")) +
        scale_color_manual(values = c("Good" = "darkgreen",  "Bad" = "#1D04C2")) +
        theme_bw() +
        #scale_y_continuous(limits = c(0, max_y_value)) +  # Dynamically set y-axis limits
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

ggline(df_subpop_clean, 
       x = "Timepoint", 
       y = "data_column", 
       color = "Recovery", 
       add = "mean_se",  # Add mean and standard error
       size = 1.2,
       ylim = c(5, 25)) +  
    labs(y = paste(cells_colname, "percentage"), 
         title = paste("Mean ± SEM of", cells_colname, "by Recovery across Timepoints")) +
    scale_color_manual(values = c("Good" = "darkgreen",  "Bad" = "#1D04C2")) +
    theme_bw() +
    theme(
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black")
    ) +
    scale_x_discrete(labels = c(
        "TP0" = "Control",
        "TP1" = "24 hours",
        "TP2" = "3-5 days",
        "TP3" = "1 month",
        "TP4" = "3 months"
    )) +
    stat_compare_means(method = "wilcox.test",  # or "t.test"
                       aes(group = Recovery),
                       label = "p.signif",  
                       size = 5,
                       label.y = 25)

ggsave(filename = "NonClassical Monocytes_Mean_SEM_RecoveryComparison.png")

results_FACS_Wilcox$TP0_Significance <- sapply(results_FACS_Wilcox$TP0_P_Value, get_significance)
results_FACS_Wilcox$TP1_Significance <- sapply(results_FACS_Wilcox$TP1_P_Value, get_significance)
results_FACS_Wilcox$TP2_Significance <- sapply(results_FACS_Wilcox$TP2_P_Value, get_significance)
results_FACS_Wilcox$TP3_Significance <- sapply(results_FACS_Wilcox$TP3_P_Value, get_significance)
results_FACS_Wilcox$TP4_Significance <- sapply(results_FACS_Wilcox$TP4_P_Value, get_significance)
write.csv(results_FACS_Wilcox, file = "results_FACS_Wilcox.csv", row.names = FALSE)


#### Linear regression -------------------------------------
output_folder <- "LinReg_FACS_Plots_Age"
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
#FACSdata <- FACSdata[!is.na(FACSdata$Age), ]

# Main loop for each column
for (i in 9:12) {
    cells_colname <- colnames(FACSdata)[i]
    
    # All samples
    plot_linear_regression1(FACSdata, cells_colname, "All", output_folder)
    
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
output_folder <- "LinReg_FACS_Plots_NHISS"
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}

NHISS_correlation_FACS <- data.frame(Cells = character(), Timepoint = character(), Sex = character(), N = numeric(), 
                                     p.value = numeric(), Coefficient = numeric(), Down95 = numeric(), 
                                     Up95 = numeric(), Intercept = numeric(), Fold = numeric(), 
                                     Max_NHISS = numeric(), Min_NHISS = numeric())
# remove Patient40 as I do not have NHISS and Sex from them
#FACSdata <- FACSdata[!is.na(FACSdata$NHISS), ]

# Main loop for each column
for (i in 9:12) {
    cells_colname <- colnames(FACSdata)[i]
    
    # All samples
    NHISS_correlation_FACS <- plot_linear_regression_NHISS(FACSdata, cells_colname, "All", output_folder, NHISS_correlation_FACS, "NHISS")
    
    # # Male samples
    # FACSdata_male <- FACSdata[FACSdata$Sex == "Male", ]
    # if (nrow(FACSdata_male) > 0) {
    #     plot_linear_regression_NHISS(FACSdata_male, cells_colname, "Male", output_folder, NHISS_correlation_FACS)
    # }
    # 
    # # Female samples
    # FACSdata_female <- FACSdata[FACSdata$Sex == "Female", ]
    # if (nrow(FACSdata_female) > 0) {
    #     plot_linear_regression_NHISS(FACSdata_female, cells_colname, "Female", output_folder, NHISS_correlation_FACS)
    # }
    # 
    # # Minor category
    # FACSdata_Minor <- FACSdata[FACSdata$Category == "MINOR", ]
    # if (nrow(FACSdata_Minor) > 0) {
    #     plot_linear_regression_NHISS(FACSdata_Minor, cells_colname, "Minor", output_folder, NHISS_correlation_FACS)
    # }
    # 
    # # Moderate category
    # FACSdata_Moderate <- FACSdata[FACSdata$Category == "MODERATE", ]
    # if (nrow(FACSdata_Moderate) > 0) {
    #     plot_linear_regression_NHISS(FACSdata_Moderate, cells_colname, "Moderate", output_folder, NHISS_correlation_FACS)
    # }
    # 
    # # Good recovery
    # FACSdata_Good <- FACSdata[FACSdata$Recovery == "Good", ]
    # FACSdata_Good <- FACSdata_Good[!is.na(FACSdata_Good$Age), ]
    # if (nrow(FACSdata_Good) > 0) {
    #     plot_linear_regression_NHISS(FACSdata_Good, cells_colname, "Good", output_folder, NHISS_correlation_FACS)
    # }
    # 
    # # Bad recovery
    # FACSdata_Bad <- FACSdata[FACSdata$Recovery == "Bad", ]
    # FACSdata_Bad <- FACSdata_Bad[!is.na(FACSdata_Bad$Age), ]
    # if (nrow(FACSdata_Bad) > 0) {
    #     plot_linear_regression_NHISS(FACSdata_Bad, cells_colname, "Bad", output_folder, NHISS_correlation_FACS)
    #}
}

# Save results to CSV
write.csv(NHISS_correlation_FACS, "NHISS_correlation_FACS.csv", row.names = FALSE)

sigNHISS_correlation_FACS <- NHISS_correlation_FACS %>% filter(NHISS_correlation_FACS$p.value <0.05)
write.csv(sigNHISS_correlation_FACS, "sigNHISS_correlation_FACS.csv", row.names = FALSE)

#### NHISS_End linear regression -----------------
output_folder <- "LinReg_FACS_Plots_NHISS_End"
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}
# # Correlation with "final" NHISS
# # Initialize NHISS_End with NA
# FACSdata$NHISS_End <- NA
# 
# # Assign NHISS_End values for all timepoints based on TP4
# TP4_indices <- which(FACSdata$Timepoint == "TP4")
# TP4_values <- FACSdata$NHISS[TP4_indices]
# 
# # Define timepoints to adjust
# timepoints_to_adjust <- c("TP3", "TP2", "TP1")
# 
# # Loop through each timepoint and assign NHISS_End
# for (tp in timepoints_to_adjust) {
#     current_indices <- which(FACSdata$Timepoint == tp)
#     min_length <- min(length(current_indices), length(TP4_indices))
#     
#     # Assign TP4 values to the current timepoint
#     FACSdata$NHISS_End[current_indices[1:min_length]] <- TP4_values[1:min_length]
#     
#     # Fill remaining rows with NA if TP4 has fewer values
#     if (length(current_indices) > min_length) {
#         FACSdata$NHISS_End[current_indices[(min_length + 1):length(current_indices)]] <- NA
#     }
# }
# 
# # Ensure TP4 values are copied correctly
# FACSdata$NHISS_End[TP4_indices] <- TP4_values
FACSdata_copy <- FACSdata

NHISS_End_correlation_FACS <- data.frame(Cells = character(), Timepoint = character(), Sex = character(), N = numeric(), 
                                         p.value = numeric(), Coefficient = numeric(), Down95 = numeric(), 
                                         Up95 = numeric(), Intercept = numeric(), Fold = numeric(), 
                                         Max_NHISS_End = numeric(), Min_NHISS_End = numeric())
# Main loop for each column
for (i in 9:12) {
    cells_colname <- colnames(FACSdata)[i]
    # Remove rows with NA in NHISS_End
    FACSdata <- FACSdata[!is.na(FACSdata$NHISS_End), ]
    
    # All samples
    NHISS_End_correlation_FACS <- plot_linear_regression_NHISS(
        FACSdata,
        cells_colname,
        "All",
        output_folder,
        NHISS_End_correlation_FACS,
        "NHISS_End" # Pass NHISS_End as a string
    )
}

# Save results to CSV
write.csv(NHISS_End_correlation_FACS, "NHISS_End_correlation_FACS.csv", row.names = FALSE)

sigNHISS_End_correlation_FACS <- NHISS_End_correlation_FACS %>% filter(NHISS_End_correlation_FACS$p.value <0.05)
write.csv(sigNHISS_End_correlation_FACS, "sigNHISS_End_correlation_FACS.csv", row.names = FALSE)

#### NHISS_Diff linear regression -----------------
output_folder <- "LinReg_FACS_Plots_NHISS_Diff"
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}

FACSdata <- FACSdata_copy

# Initialize NHISS_Diff with NA
FACSdata$NHISS_Diff <- NA

# Calculate the difference for each SampleID
FACSdata <- FACSdata %>%
    group_by(SampleID) %>%
    mutate(
        NHISS_Diff = if (all(c("TP1", "TP4") %in% Timepoint)) {
            NHISS[Timepoint == "TP1"] - NHISS[Timepoint == "TP4"]
        } else {
            NA
        }
    ) %>%
    ungroup()

FACSdata_copy <- FACSdata

NHISS_Diff_correlation_FACS <- data.frame(Cells = character(), Timepoint = character(), Sex = character(), N = numeric(), 
                                          p.value = numeric(), Coefficient = numeric(), Down95 = numeric(), 
                                          Up95 = numeric(), Intercept = numeric(), Fold = numeric(), 
                                          Max_NHISS_Diff = numeric(), Min_NHISS_Diff = numeric())

# Main loop for each column
for (i in 9:12) {
    cells_colname <- colnames(FACSdata)[i]
    # Remove rows with NA in NHISS_End
    FACSdata <- FACSdata[!is.na(FACSdata$NHISS_Diff), ]
    
    # All samples
    NHISS_Diff_correlation_FACS <- plot_linear_regression_NHISS(
        FACSdata,
        cells_colname,
        "All",
        output_folder,
        NHISS_End_correlation_FACS,
        "NHISS_Diff" # Pass NHISS_End as a string
    )
}

# Save results to CSV
write.csv(NHISS_Diff_correlation_FACS, "NHISS_Diff_correlation_FACS.csv", row.names = FALSE)

sigNHISS_Diff_correlation_FACS <- NHISS_Diff_correlation_FACS %>% filter(NHISS_Diff_correlation_FACS$p.value <0.05)
write.csv(sigNHISS_Diff_correlation_FACS, "sigNHISS_Diff_correlation_FACS.csv", row.names = FALSE)

#### NHISS_Ratio linear regression -----------------
output_folder <- "LinReg_FACS_Plots_NHISS_Ratio"
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}

FACSdata <- FACSdata_copy
# Initialize NHISS_Ratio with NA
FACSdata$NHISS_Ratio <- NA

# Calculate the difference for each SampleID
FACSdata <- FACSdata %>%
    group_by(SampleID) %>%
    mutate(
        NHISS_Ratio = if (all(c("TP1", "TP4") %in% Timepoint)) {
            (NHISS[Timepoint == "TP1"] - NHISS[Timepoint == "TP4"])/NHISS[Timepoint == "TP1"]
        } else {
            NA
        }
    ) %>%
    ungroup()

FACSdata_copy <- FACSdata

NHISS_Ratio_correlation_FACS <- data.frame(Cells = character(), Timepoint = character(), Sex = character(), N = numeric(), 
                                           p.value = numeric(), Coefficient = numeric(), Down95 = numeric(), 
                                           Up95 = numeric(), Intercept = numeric(), Fold = numeric(), 
                                           Max_NHISS_Ratio = numeric(), Min_NHISS_Ratio = numeric())

# Main loop for each column
for (i in 9:12) {
    cells_colname <- colnames(FACSdata)[i]
    # Remove rows with NA in NHISS_End
    FACSdata <- FACSdata[!is.na(FACSdata$NHISS_Ratio), ]
    
    # All samples
    NHISS_Ratio_correlation_FACS <- plot_linear_regression_NHISS(
        FACSdata,
        cells_colname,
        "All",
        output_folder,
        NHISS_End_correlation_FACS,
        "NHISS_Ratio" # Pass NHISS_End as a string
    )
}


# Save results to CSV
write.csv(NHISS_Ratio_correlation_FACS, "NHISS_Ratio_correlation_FACS.csv", row.names = FALSE)

sigNHISS_Ratio_correlation_FACS <- NHISS_Ratio_correlation_FACS %>% filter(NHISS_Ratio_correlation_FACS$p.value <0.05)
write.csv(sigNHISS_Ratio_correlation_FACS, "sigNHISS_Ratio_correlation_FACS.csv", row.names = FALSE)

# restore original dataframe, as its randomly deleting things....
FACSdata <- FACSdata_copy