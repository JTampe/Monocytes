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
        
        # Check for zero or near-zero variance
        if (any(sapply(features, function(f) length(unique(data[[f]])) < 2))) {
            stop("Feature(s) have zero or near-zero variance.")
        }
        
        set.seed(123)
        train_index <- createDataPartition(data[[target]], p = 0.8, list = FALSE)
        train_data <- data[train_index, ]
        test_data <- data[-train_index, ]
        
        trctrl <- trainControl(method = "cv", number = 5)
        
        # Handle single-feature case
        if (length(features) == 1) {
            tunegrid <- expand.grid(.mtry = 1)
        } else {
            max_mtry <- max(1, floor(sqrt(length(features))))
            min_mtry <- max(1, min(2, length(features)))
            tunegrid <- expand.grid(.mtry = min_mtry:max_mtry)
        }
        
        set.seed(123)
        rf_model <- train(
            x = train_data[, features, drop = FALSE],
            y = train_data[[target]],
            method = "rf",
            metric = "RMSE",
            trControl = trctrl,
            tuneGrid = tunegrid,
            ntree = 500,
            importance = TRUE
        )
        
        if (is.null(rf_model$finalModel)) {
            stop("Model training failed.")
        }
        
        predictions <- predict(rf_model, newdata = test_data[, features, drop = FALSE])
        
        # Calculate MAE
        mae <- mean(abs(predictions - test_data[[target]]))
        cat(paste("MAE for", prefix, ":", round(mae, 4), "\n"))
        
        
        plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
        p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
            geom_point(color = "steelblue", alpha = 0.7) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
            labs(title = paste("RF CV Predictions vs Actual -", prefix),
                 x = paste("Actual", target), y = paste("Predicted", target)) +
            theme_minimal()
        ggsave(filename = file.path(plot_folder, paste0("RF_CV_Pred_vs_Actual_", prefix, ".png")), plot = p)
 
        var_imp <- varImp(rf_model)$importance
        imp_df <- data.frame(Feature = rownames(var_imp), Importance = var_imp[, 1]) %>%
            arrange(desc(Importance))
        
        write.csv(imp_df, file.path(imp_folder, paste0("RF_CV_Feature_Importances_", prefix, ".csv")), row.names = FALSE)

        p_imp <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() +
            labs(title = paste("Variable Importance (CV) -", prefix), x = "Feature", y = "Importance") +
            theme_minimal()
        ggsave(filename = file.path(imp_folder, paste0("VarImpPlot_CV_", prefix, "_ggplot.png")), plot = p_imp)
        
        return(list(rf_model = rf_model, imp_df = imp_df, mae = mae))
    }, error = function(e) {
        return(list(error = e$message))
    })
}

rf_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")

# Age & CD33 RF model --------------------------
AgeCD33_RF_model <- run_rf_cv_analysis(data = matrix_meanCapZ %>% filter(Timepoint == "TP2"),
                                       features = c("Age","CD33"), 
                                       target = "NHISS_Diff", 
                                       prefix = "RF_AgeCD33_Model", 
                                       plot_folder = "RF_AgeCD33_Plots", 
                                       imp_folder = "RF_AgeCD33_Importance_Plots")

# get MAE just from the model:
AgeCD33_RF_model$rf_model[4]$results$MAE[1]

# Create the data frame
prediction_data <- expand.grid(
    Age = c(52, 60, 70, 80, 84),
    CD33 = c(-0.5, 0, 0.5, 1)
)

# column for the predicted values
prediction_data$NHISS_Diff_predicted <- "NA"

# example to calculate the NHISS_Diff for the first row of example values
# predict(AgeCD33_RF_model$rf_model, newdata = prediction_data[1,])

# Loop to predict for all rows
for (i in 1:nrow(prediction_data)) {
    prediction_data$NHISS_Diff_predicted[i] <- predict(AgeCD33_RF_model$rf_model, newdata = prediction_data[i, , drop = FALSE])
}

View(prediction_data)

# Age & Sex RF model --------------------------

AgeSexExtended_cleaned <- Metadata_NeuroTest %>%
    filter(Timepoint == "TP2") %>%
    select(Age, Sex, NHISS_Ratio, NHISS_Diff) %>%
    drop_na()

# Diff
AgeSexExtended_RF_model_Diff <- run_rf_cv_analysis(data = Metadata_NeuroTest %>% filter(Timepoint == "TP2"),
                                       features = c("Age","Sex", "Barthel", "HADS_Depression","MoCA","NHISS", "mRS", "HADS_Anxiety"), 
                                       target = "NHISS_Diff", 
                                       prefix = "RF_AgeSexExtended_Model", 
                                       plot_folder = "RF_AgeSexExtended_Plots", 
                                       imp_folder = "RF_AgeSexExtended_Importance_Plots")
#AgeSexExtended_RF_model_Diff$rf_model

# Ratio
AgeSexExtended_RF_model_Ratio <- run_rf_cv_analysis(data = Metadata_NeuroTest %>% filter(Timepoint == "TP2"),
                                                   features =  c("Age","Sex", "Barthel", "HADS_Depression","MoCA","NHISS", "mRS", "HADS_Anxiety"), 
                                                   target = "NHISS_Ratio", 
                                                   prefix = "RF_AgeSexExtended_Model", 
                                                   plot_folder = "RF_AgeSexExtended_Plots", 
                                                   imp_folder = "RF_AgeSexExtended_Importance_Plots")
AgeSexExtended_RF_model_Ratio$rf_model
