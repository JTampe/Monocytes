# Define the function with improved error handling and dynamic tuning grid
run_rf_cv_analysis <- function(data, features, target, prefix = "RF_CV_Model", plot_folder = "RF_Plots", imp_folder = "RF_Importance_Plots") {
    library(randomForest)
    library(ggplot2)
    library(dplyr)
    library(caret)
    
    tryCatch({
        if (!dir.exists(plot_folder)) dir.create(plot_folder)
        if (!dir.exists(imp_folder)) dir.create(imp_folder)
        
        data <- data %>%
            select(all_of(c(features, target))) %>%
            drop_na()
        
        if (nrow(data) < 10) {
            stop("Not enough data to train the model.")
        }
        
        set.seed(123)
        train_index <- createDataPartition(data[[target]], p = 0.8, list = FALSE)
        train_data <- data[train_index, ]
        test_data <- data[-train_index, ]
        
        trctrl <- trainControl(method = "cv", number = 5)
        
        # max_mtry <- max(1, floor(sqrt(length(features))))
        # min_mtry <- min(2, length(features))
        # tunegrid <- expand.grid(.mtry = min_mtry:max_mtry)
        max_mtry <- max(1, floor(sqrt(length(features))))
        min_mtry <- max(1, min(2, length(features)))  # Ensure min_mtry is at least 1
        tunegrid <- expand.grid(.mtry = min_mtry:max_mtry)
        
        set.seed(123)
        rf_model <- train(
            x = train_data[, features],
            y = train_data[[target]],
            method = "rf",
            metric = "RMSE",
            trControl = trctrl,
            tuneGrid = tunegrid,
            ntree = 500,
            importance = TRUE
        )
        
        predictions <- predict(rf_model, newdata = test_data[, features])
        
        # plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
        # p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
        #     geom_point(color = "steelblue", alpha = 0.7) +
        #     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
        #     labs(title = paste("RF CV Predictions vs Actual -", prefix),
        #          x = paste("Actual", target), y = paste("Predicted", target)) +
        #     theme_minimal()
        # ggsave(filename = file.path(plot_folder, paste0("RF_CV_Pred_vs_Actual_", prefix, ".png")), plot = p)
        # 
        var_imp <- varImp(rf_model)$importance
        imp_df <- data.frame(Feature = rownames(var_imp), Importance = var_imp[,1]) %>%
            arrange(desc(Importance))
        
        # write.csv(imp_df, file.path(imp_folder, paste0("RF_CV_Feature_Importances_", prefix, ".csv")), row.names = FALSE)
        # 
        # p_imp <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
        #     geom_bar(stat = "identity", fill = "steelblue") +
        #     coord_flip() +
        #     labs(title = paste("Variable Importance (CV) -", prefix), x = "Feature", y = "Importance") +
        #     theme_minimal()
        # ggsave(filename = file.path(imp_folder, paste0("VarImpPlot_CV_", prefix, "_ggplot.png")), plot = p_imp)
        
        return(list(rf_model = rf_model, imp_df = imp_df))
    }, error = function(e) {
        return(list(error = e$message))
    })
}

run_rf_feature_combinations <- function(data, features, target) {
    library(dplyr)
    library(purrr)
    library(Metrics)
    
    RF_Results <- data.frame(
        Prefix = character(),
        Features = character(),
        MAE = numeric(),
        MSE = numeric(),
        Comments = character(),
        stringsAsFactors = FALSE
    )
    
    for (k in 1:length(features)) {
        feature_combos <- combn(features, k, simplify = FALSE)
        
        for (feat_set in feature_combos) {
            prefix <- paste0("RF_k", k, "_", paste(feat_set, collapse = "_"))
            
            result <- run_rf_cv_analysis(
                data = data,
                features = feat_set,
                target = target,
                prefix = prefix
            )
            
            if (!is.null(result$error)) {
                RF_Results <- rbind(RF_Results, data.frame(
                    Prefix = prefix,
                    Features = paste(feat_set, collapse = ", "),
                    MAE = NA,
                    MSE = NA,
                    Comments = paste("Error:", result$error),
                    stringsAsFactors = FALSE
                ))
            } else {
                test_data <- data %>%
                    select(all_of(c(feat_set, target))) %>%
                    drop_na()
                
                set.seed(123)
                train_index <- createDataPartition(test_data[[target]], p = 0.8, list = FALSE)
                
                if (length(train_index) < nrow(test_data)) {
                    test_data <- test_data[-train_index, ]
                    
                    if (nrow(test_data) > 0) {
                        predictions <- predict(result$rf_model, newdata = test_data[, feat_set])
                        actuals <- test_data[[target]]
                        
                        mae_val <- mae(actuals, predictions)
                        mse_val <- mse(actuals, predictions)
                        
                        RF_Results <- rbind(RF_Results, data.frame(
                            Prefix = prefix,
                            Features = paste(feat_set, collapse = ", "),
                            MAE = mae_val,
                            MSE = mse_val,
                            Comments = "",
                            stringsAsFactors = FALSE
                        ))
                    } else {
                        RF_Results <- rbind(RF_Results, data.frame(
                            Prefix = prefix,
                            Features = paste(feat_set, collapse = ", "),
                            MAE = NA,
                            MSE = NA,
                            Comments = "Empty test set after partitioning",
                            stringsAsFactors = FALSE
                        ))
                    }
                } else {
                    RF_Results <- rbind(RF_Results, data.frame(
                        Prefix = prefix,
                        Features = paste(feat_set, collapse = ", "),
                        MAE = NA,
                        MSE = NA,
                        Comments = "All rows assigned to training set",
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }
    
    return(RF_Results)
}

# Prepare data
rf_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")

#write.csv(rf_data, "rf_data_TP2.csv", row.names = FALSE)

# Define feature set
clinical_features <- c("Barthel", "HADS_Depression","MoCA", "Sex","NHISS", "mRS", "HADS_Anxiety"
                       #,"SPAN", "Sex","NHISS", "mRS", "HADS_Anxiety"
                       )
gene_features <- c("CD33", "CD169", "TNFA", "TGFB1", "CD91" #, "STAT6"
                   )
all_features <- c(clinical_features, gene_features)

combined_RF_results2 <- run_rf_feature_combinations(
    data = rf_data,
    features = all_features,
    target = "NHISS_Diff"
)
View(combined_RF_results2)

# only known risk factors without testing:
clinical_RF_results <- run_rf_feature_combinations(
    data = rf_data,
    features = c("Age", "Sex"),
    target = "NHISS_Diff"
)
# run all clinical features
clinical_features <- c("Barthel", "HADS_Depression","MoCA", "NHISS", "mRS", "HADS_Anxiety")
clinical_RF_results <- run_rf_feature_combinations(
    data = rf_data,
    features = clinical_features,
    target = "NHISS_Diff"
)
View(clinical_RF_results)
write_xlsx(clinical_RF_results,"RF_result_clinical.xlsx")
### only gene results
gene_features <- c("CD33", "CD169", "TNFA", "TGFB1", "CD91", "STAT6")
gene_RF_results <- run_rf_feature_combinations(
    data = rf_data,
    features = gene_features,
    target = "NHISS_Diff"
)
View(gene_RF_results)
write_xlsx(gene_RF_results,"RF_result_gene.xlsx")
# clinical data plus CD33, CD91
clinical_RF_resultsCD33CD91 <- run_rf_feature_combinations(
    data = rf_data,
    features = c("Barthel", "HADS_Depression","MoCA","NHISS", "mRS", "HADS_Anxiety","CD33", "CD91"),
    target = "NHISS_Diff"
)
View(clinical_RF_resultsCD33CD91)
write_xlsx(clinical_RF_resultsCD33CD91,"RF_result_clinicalCD33CD91.xlsx")
# clinical data plus CD33, CD91, TGFB1
clinical_RF_resultsCD33CD91TGFB1 <- run_rf_feature_combinations(
    data = rf_data,
    features = c("Barthel", "HADS_Depression","MoCA","NHISS", "mRS", "HADS_Anxiety","CD33", "CD91","TGFB1"),
    target = "NHISS_Diff"
)
View(clinical_RF_resultsCD33CD91TGFB1)
write_xlsx(clinical_RF_resultsCD33CD91TGFB1,"RF_result_clinicalCD33CD91TGFB1.xlsx")

# clinical data plus CD33, CD91, TGFB1, CD169
clinical_RF_resultsCD33CD91TGFB1CD169TNFA <- run_rf_feature_combinations(
    data = rf_data,
    features = c("Barthel", "HADS_Depression","MoCA", "NHISS", "mRS", "HADS_Anxiety","CD33", "CD91", "TGFB1", "CD169"),
    target = "NHISS_Diff"
)
View(clinical_RF_resultsCD33CD91TGFB1CD169)
write_xlsx(clinical_RF_resultsCD33CD91TGFB1CD169,"RF_result_clinicalCD33CD91TGFB1CD169.xlsx")

# clinical data plus CD33, CD91, TGFB1, CD169, TNFA, STAT6
clinical_RF_resultsCD33CD91TGFB1CD169TNFASTAT6 <- run_rf_feature_combinations(
    data = rf_data,
    features = c("Barthel", "HADS_Depression","MoCA", "NHISS", "mRS", "HADS_Anxiety","CD33", "CD91", "TGFB1", "CD169","TNFA","STAT6"),
    target = "NHISS_Diff"
)
# View(clinical_RF_resultsCD33CD91TGFB1CD169TNFASTAT6)
# write_xlsx(clinical_RF_resultsCD33CD91TGFB1CD169TNFASTAT6,"RF_result_clinicalCD33CD91TGFB1CD169TNFASTAT6.xlsx")

# write_xlsx(clinical_RF_resultsCD33CD91TGFB1CD169TNFASTAT6, "clinical_RF_resultsCD33CD91TGFB1CD169TNFASTAT6.xlsx")
write_xlsx(clinical_RF_resultsCD33CD91TGFB1CD169TNFASTAT6, "RF_result_clincalCD33CD91TGFB1CD169TNFASTAT6.xlsx")

GenesAgeSex_RF_results1 <- run_rf_feature_combinations(
    data = rf_data,
    features = c("Age", "CD38","CD33",  "CD169", "TNFA", "TGFB1", "CD91"),
    target = "NHISS_Diff"
)
View(GenesAgeSex_RF_results1)
