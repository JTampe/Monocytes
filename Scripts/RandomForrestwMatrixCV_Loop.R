# Define the function
run_rf_cv_analysis <- function(data, features, target, prefix = "RF_CV_Model", plot_folder = "RF_Plots", imp_folder = "RF_Importance_Plots") {
    library(randomForest)
    library(ggplot2)
    library(dplyr)
    library(caret)
    
    if (!dir.exists(plot_folder)) dir.create(plot_folder)
    if (!dir.exists(imp_folder)) dir.create(imp_folder)
    
    data <- data %>%
        select(all_of(c(features, target))) %>%
        drop_na()
    
    if (nrow(data) < 10) {
        warning("Not enough data to train the model.")
        return(NULL)
    }
    
    set.seed(123)
    train_index <- createDataPartition(data[[target]], p = 0.8, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    trctrl <- trainControl(method = "cv", number = 5)
    
    if (length(features) >= 2) {
        tunegrid <- expand.grid(.mtry = seq(2, floor(sqrt(length(features))), by = 1))
    } else {
        tunegrid <- expand.grid(.mtry = 1)
    }
    
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
    
    plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
    p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
        geom_point(color = "steelblue", alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
        labs(title = paste("RF CV Predictions vs Actual -", prefix),
             x = paste("Actual", target), y = paste("Predicted", target)) +
        theme_minimal()
    ggsave(filename = file.path(plot_folder, paste0("RF_CV_Pred_vs_Actual_", prefix, ".png")), plot = p)
    
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
}

# Prepare data
rf_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")

# Load required libraries
library(dplyr)
library(purrr)
library(Metrics)

# Define feature set for the clinical features
#all_features <- c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression")

# Initialize results
RF_Results <- data.frame(
    Prefix = character(),
    Features = character(),
    MAE = numeric(),
    MSE = numeric(),
    Comments = character(),
    stringsAsFactors = FALSE
)

# Loop through combinations
for (k in 1:length(all_features)) {
    feature_combos <- combn(all_features, k, simplify = FALSE)
    
    for (features in feature_combos) {
        prefix <- paste0("RF_k", k, "_", paste(features, collapse = "_"))
        
        tryCatch({
            result <- run_rf_cv_analysis(
                data = rf_data,
                features = features,
                target = "NHISS_Diff",
                prefix = prefix
            )
            
            if (!is.null(result)) {
                test_data <- rf_data %>%
                    select(all_of(c(features, "NHISS_Diff"))) %>%
                    drop_na()
                
                set.seed(123)
                train_index <- createDataPartition(test_data[["NHISS_Diff"]], p = 0.8, list = FALSE)
                
                if (length(train_index) < nrow(test_data)) {
                    test_data <- test_data[-train_index, ]
                    
                    if (nrow(test_data) > 0) {
                        predictions <- predict(result$rf_model, newdata = test_data[, features])
                        actuals <- test_data[["NHISS_Diff"]]
                        
                        mae_val <- mae(actuals, predictions)
                        mse_val <- mse(actuals, predictions)
                        
                        RF_Results <- rbind(RF_Results, data.frame(
                            Prefix = prefix,
                            Features = paste(features, collapse = ", "),
                            MAE = mae_val,
                            MSE = mse_val,
                            Comments = "",
                            stringsAsFactors = FALSE
                        ))
                    } else {
                        RF_Results <- rbind(RF_Results, data.frame(
                            Prefix = prefix,
                            Features = paste(features, collapse = ", "),
                            MAE = NA,
                            MSE = NA,
                            Comments = "Empty test set after partitioning",
                            stringsAsFactors = FALSE
                        ))
                    }
                } else {
                    RF_Results <- rbind(RF_Results, data.frame(
                        Prefix = prefix,
                        Features = paste(features, collapse = ", "),
                        MAE = NA,
                        MSE = NA,
                        Comments = "All rows assigned to training set",
                        stringsAsFactors = FALSE
                    ))
                }
            } else {
                RF_Results <- rbind(RF_Results, data.frame(
                    Prefix = prefix,
                    Features = paste(features, collapse = ", "),
                    MAE = NA,
                    MSE = NA,
                    Comments = "Model training returned NULL",
                    stringsAsFactors = FALSE
                ))
            }
        }, error = function(e) {
            RF_Results <- rbind(RF_Results, data.frame(
                Prefix = prefix,
                Features = paste(features, collapse = ", "),
                MAE = NA,
                MSE = NA,
                Comments = paste("Error:", e$message),
                stringsAsFactors = FALSE
            ))
        })
    }
}

# Define feature set for the clinical features
all_features <- c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", #"SPAN",
                  "CD33", "CD169", "TNFA", "TGFB1", "CD91", "STAT6")

