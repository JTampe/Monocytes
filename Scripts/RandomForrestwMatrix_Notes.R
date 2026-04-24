# function for random forrest
run_rf_analysis <- function(data, features, target, prefix = "RF_Model", plot_folder = "RF_Plots", imp_folder = "RF_Importance_Plots", ntree = 500) {
    library(randomForest)
    library(ggplot2)
    library(dplyr)
    library(caret)
    
    # Create output folders if they don't exist
    if (!dir.exists(plot_folder)) dir.create(plot_folder)
    if (!dir.exists(imp_folder)) dir.create(imp_folder)
    
    # Prepare data
    data <- data %>%
        select(all_of(c(features, target))) %>%
        drop_na()
    
    if (nrow(data) < 10) {
        warning("Not enough data to train the model.")
        return(NULL)
    }
    
    # Split data
    set.seed(123)
    train_index <- createDataPartition(data[[target]], p = 0.8, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    # Train model
    rf_model <- randomForest(as.formula(paste(target, "~ .")), data = train_data, importance = TRUE, ntree = ntree)
    
    # Predict
    predictions <- predict(rf_model, newdata = test_data)
    
    # Save prediction plot
    plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
    p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
        geom_point(color = "steelblue", alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
        labs(title = paste("RF Predictions vs Actual -", prefix),
             x = paste("Actual", target), y = paste("Predicted", target)) +
        theme_minimal()
    ggsave(filename = file.path(plot_folder, paste0("RF_Pred_vs_Actual_", prefix, ".png")), plot = p)
    
    # Extract importance
    imp <- importance(rf_model)[, 1]
    imp_df <- data.frame(Feature = names(imp), Importance = as.vector(imp)) %>%
        arrange(desc(Importance))
    
    # Save importance table
    write.csv(imp_df, file.path(imp_folder, paste0("RF_Feature_Importances_", prefix, ".csv")), row.names = FALSE)
    
    # Save ggplot version of importance
    p_imp <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = paste("Variable Importance -", prefix), x = "Feature", y = "Mean Decrease in Accuracy") +
        theme_minimal()
    ggsave(filename = file.path(imp_folder, paste0("VarImpPlot_", prefix, "_ggplot.png")), plot = p_imp)
    
    return(list(rf_model = rf_model, imp_df = imp_df))
}

### Random Forrest Analysis over all --------------------------------------------------

# extracte genes of the data set
genes_tp2 <- unique(dataFINALmean$Gene)
dataFINALmean  <- merge(dataFINALmean, Age_Sex, by = "SampleID")

# make a separate data frame for each value:
Data_meanCapZ <- dataFINALmean  %>% select(all_of(c("SampleID", "Timepoint", "Subpopulation", "Gene", "mean_capZ"))) 
Data_sdCapZ <- dataFINALmean  %>% select(all_of(c("SampleID", "Timepoint", "Subpopulation", "Gene", "sd_capZ"))) 
#Data_meanZ <- dataFINALmean  %>% select(all_of(c("SampleID", "Timepoint", "Subpopulation", "Gene", "mean_Z"))) 
#Data_sdZ <- dataFINALmean  %>% select(all_of(c("SampleID", "Timepoint", "Subpopulation", "Gene", "sd_Z"))) 

# make each gene a column and add them all listed to the features
matrix_meanCapZ <-  pivot_wider(Data_meanCapZ , names_from = Gene, values_from = mean_capZ)
matrix_sdCapZ <-  pivot_wider(Data_sdCapZ , names_from = Gene, values_from = sd_capZ)
#matrix_meanZ <-  pivot_wider(Data_meanZ , names_from = Gene, values_from = mean_Z)
#matrix_sdZ <-  pivot_wider(Data_sdZ , names_from = Gene, values_from = sd_Z)

# merge with the Neurological tests
matrix_meanCapZ <- merge(matrix_meanCapZ, Metadata_NeuroTest, by = c("SampleID", "Timepoint"), all.x = TRUE)
matrix_sdCapZ <- merge(matrix_sdCapZ, Metadata_NeuroTest, by = c("SampleID", "Timepoint"), all.x = TRUE)

rf_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")

result <- run_rf_analysis(
    data = rf_data,
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", 
                 as.vector(genes_tp2)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff"
)

# Access the model and importance table
rf_model <- result$rf_model
imp_df <- result$imp_df

important_features <- imp_df %>%
    filter(Importance > 2) %>%
    pull(Feature)

result_important <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2"),
    features = important_features,
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_important2"
)
