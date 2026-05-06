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

Age_correlation_capZ$Sex <- "All"
Age_correlation_capZ$Category <- "All"
Age_correlation_Sex_capZ$Category <- "All"
Age_correlation_Category_capZ$Sex<- "All"

LinReg_Age_capZ <- rbind(Age_correlation_capZ, Age_correlation_Sex_capZ, Age_correlation_Category_capZ)
LinReg_Age_sigGenes <- LinReg_Age_capZ %>% filter(LinReg_Age_capZ$p.value <0.05)
write.csv(LinReg_Age_sigGenes, "LinReg_Results_capZ/LinReg_Age_sigGenes.csv", row.names = FALSE)
