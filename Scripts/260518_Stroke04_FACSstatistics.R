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
# write_xlsx(ANOVA_FACSsex, "ANOVA_FACSsex.xlsx")
# write_xlsx(ANOVA_FACSrecovery, "ANOVA_FACSrecovery.xlsx")

ANOVA_FACSall$Category <- "All"
ANOVA_FACScombined <- rbind (ANOVA_FACSall,ANOVA_FACSsex)
ANOVA_FACScombined <- rbind (ANOVA_FACScombined,ANOVA_FACSrecovery)
write_xlsx(ANOVA_FACScombined, "ANOVA_FACScombined.xlsx")

#### Mean & SEM FACS----------------------------------------------------------------------------------------------------
# folder <- "MEAN_SEM_FACS_Sex"
# createFolder(folder)
# plot_save <- "y"
# 
# FACSdata <- FACSdata %>%
#     mutate(
#         Sex = if_else(Timepoint == "TP0" & Sex == "Female",
#                       "Female_Ctr",
#                       Sex)
#     )
# FACSdata <- FACSdata %>%
#     mutate(
#         Sex = if_else(Timepoint == "TP0" & Sex == "Male",
#                       "Male_Ctr",
#                       Sex)
#     )
# 
# # Initialize an empty dataframe to store Wilcoxon test results
# results_FACS_Wilcox <- data.frame(Cells = character(), Focus = character(), 
#                                   TP0_P_Value = numeric(), TP1_P_Value = numeric(), 
#                                   TP2_P_Value = numeric(), TP3_P_Value = numeric(), 
#                                   TP4_P_Value = numeric(), stringsAsFactors = FALSE)
# 
# for (i in 9:12) {
#     cells_colname <- colnames(FACSdata)[i]
#     
#     # Filter, clean, and rename the data column for the current column
#     df_subpop_clean <- FACSdata %>%
#         filter(!is.na(Sex), !is.na(Timepoint)) %>%
#         select(Timepoint, Sex, !!sym(cells_colname)) %>%
#         rename(data_column = !!sym(cells_colname))
#     
#     p_values <- c(TP0_P_Value = NA, TP1_P_Value = NA, TP2_P_Value = NA, 
#                   TP3_P_Value = NA, TP4_P_Value = NA)
#     
#     # Calculate Wilcoxon test p-values for each timepoint
#     for (tp in unique(df_subpop_clean$Timepoint)) {
#         df_timepoint <- subset(df_subpop_clean, Timepoint == tp)
#         
#         if (length(unique(df_timepoint$Sex)) == 2) { # Ensure both sexes are present
#             wilcox_test <- wilcox.test(data_column ~ Sex, data = df_timepoint)
#             p_values[paste0(tp, "_P_Value")] <- wilcox_test$p.value
#         }
#     }
#     
#     # Append p-values to results dataframe
#     results_FACS_Wilcox <- rbind(results_FACS_Wilcox, 
#                                  data.frame(Cells = cells_colname, Focus = "Male vs. Female", 
#                                             t(p_values), stringsAsFactors = FALSE))
#     
#     # Generate line plot with mean ± SEM for each sex across Timepoints
#     p <- ggline(df_subpop_clean, 
#                 x = "Timepoint", 
#                 y = "data_column", 
#                 color = "Sex", 
#                 add = "mean_se",  # Add mean and standard error
#                 size = 1.2) +  
#         labs(y = paste(cells_colname, "percentage"), 
#              title = paste("Mean ± SEM of", cells_colname, "by Sex across Timepoints")) +
#         scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2",
#                                       "Male_Ctr" = "#A67C00", "Female_Ctr" = "#1D04C2")) +
#         theme_bw() +
#         theme(
#             panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
#             panel.grid.minor = element_blank(), 
#             panel.border = element_blank(),
#             axis.line = element_line(color = "black")
#         ) +
#         scale_x_discrete(labels = c(
#             "TP0" = "Control",
#             "TP1" = "24 hours",
#             "TP2" = "3-5 days",
#             "TP3" = "1 month",
#             "TP4" = "3 months"
#         )) +
#         stat_compare_means(method = "wilcox.test",  # or "t.test"
#                            aes(group = Sex),
#                            label = "p.signif",  
#                            size = 5) 
#     
#     # Save plot if specified
#     if (plot_save == "y") {
#         file_name_facet <- paste0(folder, "/", cells_colname, "_Mean_SEM_SexComparison.png")
#         ggsave(filename = file_name_facet, plot = p)
#     }
# }
# 
# folder <- "MEAN_SEM_FACS_Category"
# createFolder(folder)
# plot_save <- "y"
# 
# FACSdata <- FACSdata %>%
#     mutate(
#         Category = if_else(Timepoint == "TP0",
#                       "Control",
#                       Category)
#     )
# 
# for (i in 9:12) {
#     cells_colname <- colnames(FACSdata)[i]
#     
#     # Filter, clean, and rename the data column for the current column
#     df_subpop_clean <- FACSdata %>%
#         filter(!is.na(Category), !is.na(Timepoint)) %>%
#         select(Timepoint, Category, !!sym(cells_colname)) %>%
#         rename(data_column = !!sym(cells_colname))
#     
#     p_values <- c(TP0_P_Value = NA, TP1_P_Value = NA, TP2_P_Value = NA, 
#                   TP3_P_Value = NA, TP4_P_Value = NA)
#     
#     # Calculate Wilcoxon test p-values for each timepoint
#     for (tp in unique(df_subpop_clean$Timepoint)) {
#         df_timepoint <- subset(df_subpop_clean, Timepoint == tp)
#         
#         if (length(unique(df_timepoint$Category)) == 2) { # Ensure both sexes are present
#             wilcox_test <- wilcox.test(data_column ~ Category, data = df_timepoint)
#             p_values[paste0(tp, "_P_Value")] <- wilcox_test$p.value
#         }
#     }
#     
#     # Append p-values to results dataframe
#     results_FACS_Wilcox <- rbind(results_FACS_Wilcox, 
#                                  data.frame(Cells = cells_colname, Focus = "Minor vs. Moderate", 
#                                             t(p_values), stringsAsFactors = FALSE))
#     
#     # Generate line plot with mean ± SEM for each sex across Timepoints
#     p <- ggline(df_subpop_clean, 
#                 x = "Timepoint", 
#                 y = "data_column", 
#                 color = "Category", 
#                 add = "mean_se",  # Add mean and standard error
#                 size = 1.2) +  
#         labs(y = paste(cells_colname, "percentage"), 
#              title = paste("Mean ± SEM of", cells_colname, "by Category across Timepoints")) +
#         scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606")) +
#         scale_x_discrete(labels = c(
#             "TP0" = "Control",
#             "TP1" = "24 hours",
#             "TP2" = "3-5 days",
#             "TP3" = "1 month",
#             "TP4" = "3 months"
#         )) +
#         theme_bw() +
#         theme(
#             panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
#             panel.grid.minor = element_blank(), 
#             panel.border = element_blank(),
#             axis.line = element_line(color = "black")
#         ) +
#         stat_compare_means(method = "wilcox.test",  # or "t.test"
#                            aes(group = Category),
#                            label = "p.signif",  
#                            size = 5) 
#     
#     # Save plot if specified
#     if (plot_save == "y") {
#         file_name_facet <- paste0(folder, "/", cells_colname, "_Mean_SEM_CategoryComparison.png")
#         ggsave(filename = file_name_facet, plot = p)
#     }
# }
# 
# folder <- "MEAN_SEM_FACS_Recovery"
# createFolder(folder)
# plot_save <- "y"
# 
# FACSdata <- FACSdata %>%
#     mutate(
#         Recovery = if_else(Timepoint == "TP0",
#                            "Control",
#                            Recovery)
#     )
# 
# 
# for (i in 9:12) {
#     cells_colname <- colnames(FACSdata)[i]
#     
#     # Filter, clean, and rename the data column for the current column
#     df_subpop_clean <- FACSdata %>%
#         filter(!is.na(Recovery), !is.na(Timepoint)) %>%
#         select(Timepoint, Recovery, !!sym(cells_colname)) %>%
#         rename(data_column = !!sym(cells_colname))
#     
#     p_values <- c(TP0_P_Value = NA, TP1_P_Value = NA, TP2_P_Value = NA, 
#                   TP3_P_Value = NA, TP4_P_Value = NA)
#     
#     for (tp in unique(df_subpop_clean$Timepoint)) {
#         df_timepoint <- subset(df_subpop_clean, Timepoint == tp)
#         
#         if (length(unique(df_timepoint$Recovery)) == 2) { # Ensure both sexes are present
#             wilcox_test <- wilcox.test(data_column ~ Recovery, data = df_timepoint)
#             p_values[paste0(tp, "_P_Value")] <- wilcox_test$p.value
#         }
#     }
#     
#     results_FACS_Wilcox <- rbind(results_FACS_Wilcox, 
#                                  data.frame(Cells = cells_colname, Focus = "Good vs. Bad Recovery", 
#                                             t(p_values), stringsAsFactors = FALSE))
#     
#     max_y_value <- mean(df_subpop_clean$data_column, na.rm = TRUE)+5
#     
#     p <- ggline(df_subpop_clean, 
#                 x = "Timepoint", 
#                 y = "data_column", 
#                 color = "Recovery", 
#                 add = "mean_se",  # Add mean and standard error
#                 size = 1.2) +  
#         labs(y = paste(cells_colname, "percentage"), 
#              title = paste("Mean ± SEM of", cells_colname, "by Recovery across Timepoints")) +
#         scale_color_manual(values = c("Good" = "darkgreen",  "Bad" = "#1D04C2",
#                                       "Control" = "grey")) +
#         theme_bw() +
#         #scale_y_continuous(limits = c(0, max_y_value)) +  # Dynamically set y-axis limits
#         theme(
#             panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
#             panel.grid.minor = element_blank(), 
#             panel.border = element_blank(),
#             axis.line = element_line(color = "black")
#         ) +
#         stat_compare_means(method = "wilcox.test",  # or "t.test"
#                            aes(group = Recovery),
#                            label = "p.signif",  
#                            size = 5) 
#     
#     # Save plot if specified
#     if (plot_save == "y") {
#         file_name_facet <- paste0(folder, "/", cells_colname, "_Mean_SEM_RecoveryComparison.png")
#         ggsave(filename = file_name_facet, plot = p)
#     }
# }
# 
# ggline(df_subpop_clean, 
#        x = "Timepoint", 
#        y = "data_column", 
#        color = "Recovery", 
#        add = "mean_se",  # Add mean and standard error
#        size = 1.2,
#        ylim = c(5, 25)) +  
#     labs(y = paste(cells_colname, "percentage"), 
#          title = paste("Mean ± SEM of", cells_colname, "by Recovery across Timepoints")) +
#     scale_color_manual(values = c("Good" = "darkgreen",  "Bad" = "#1D04C2")) +
#     theme_bw() +
#     theme(
#         panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank(),
#         axis.line = element_line(color = "black")
#     ) +
#     scale_x_discrete(labels = c(
#         "TP0" = "Control",
#         "TP1" = "24 hours",
#         "TP2" = "3-5 days",
#         "TP3" = "1 month",
#         "TP4" = "3 months"
#     )) +
#     stat_compare_means(method = "wilcox.test",  # or "t.test"
#                        aes(group = Recovery),
#                        label = "p.signif",  
#                        size = 5,
#                        label.y = 25)
# 
# ggsave(filename = "NonClassical Monocytes_Mean_SEM_RecoveryComparison.png")
# 
# results_FACS_Wilcox$TP0_Significance <- sapply(results_FACS_Wilcox$TP0_P_Value, get_significance)
# results_FACS_Wilcox$TP1_Significance <- sapply(results_FACS_Wilcox$TP1_P_Value, get_significance)
# results_FACS_Wilcox$TP2_Significance <- sapply(results_FACS_Wilcox$TP2_P_Value, get_significance)
# results_FACS_Wilcox$TP3_Significance <- sapply(results_FACS_Wilcox$TP3_P_Value, get_significance)
# results_FACS_Wilcox$TP4_Significance <- sapply(results_FACS_Wilcox$TP4_P_Value, get_significance)
# write.csv(results_FACS_Wilcox, file = "results_FACS_Wilcox.csv", row.names = FALSE)
