library(dplyr)
library(randomForest)
library(caret)
library(ggplot2)
library(tidyr)

# extracte genes of the data set
genes_tp2 <- unique(dataFINALmean$Gene)

# make a separate data frame for each value:
Data_meanCapZ <- dataFINALmean  %>% select(all_of(c("SampleID", "Timepoint", "Subpopulation", "Gene", "mean_capZ"))) 
Data_sdCapZ <- dataFINALmean  %>% select(all_of(c("SampleID", "Timepoint", "Subpopulation", "Gene", "sd_capZ"))) 
#Data_meanZ <- dataFINALmean  %>% select(all_of(c("SampleID", "Timepoint", "Subpopulation", "Gene", "mean_Z"))) 
#Data_sdZ <- dataFINALmean  %>% select(all_of(c("SampleID", "Timepoint", "Subpopulation", "Gene", "sd_Z"))) 

    #filter 
#Data_meanZ <- as.data.frame(SampleID=dataFINALmean$SampleID, Timepoint=dataFINALmean$Timepoint, Subpopulation=dataFINALmean$Subpopulation, Gene=dataFINALmean$Gene, 
#                    Value=dataFINALmean$mean_capZ)

# make each gene a column and add them all listed to the features
matrix_meanCapZ <-  pivot_wider(Data_meanCapZ , names_from = Gene, values_from = mean_capZ)
matrix_sdCapZ <-  pivot_wider(Data_sdCapZ , names_from = Gene, values_from = sd_capZ)
#matrix_meanZ <-  pivot_wider(Data_meanZ , names_from = Gene, values_from = mean_Z)
#matrix_sdZ <-  pivot_wider(Data_sdZ , names_from = Gene, values_from = sd_Z)

# merge with the Neurological tests
matrix_meanCapZ <- merge(matrix_meanCapZ, Metadata_NeuroTest, by = c("SampleID", "Timepoint"), all.x = TRUE)
matrix_sdCapZ <- merge(matrix_sdCapZ, Metadata_NeuroTest, by = c("SampleID", "Timepoint"), all.x = TRUE)

# add the metadata of of Patients
metadataP <- metadataP %>% select(all_of(c("SampleID", "Sex", "Age"))) 
matrix_meanCapZ <- merge(matrix_meanCapZ , metadataP, by = "SampleID", all.x = TRUE)
matrix_sdCapZ <- merge(matrix_sdCapZ , metadataP, by = "SampleID", all.x = TRUE)

# Create output folders
plot_folder <- "RF_Plots"
imp_folder <- "RF_Importance_Plots"
if (!dir.exists(plot_folder)) dir.create(plot_folder)
if (!dir.exists(imp_folder)) dir.create(imp_folder)


# Filter for TP2
model_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")

# Define features and target
features <- c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression",
              as.vector(genes_tp2))
#target <- "NHISS_Ratio"
target <- "NHISS_Diff"
#target <- "NHISS_End"

    # Convert categorical variables to factors
    model_data$Sex <- as.factor(model_data$Sex)
    #model_data$Category <- as.factor(model_data$Category)
    model_data$Subpopulation <- as.factor(model_data$Subpopulation)
    
    
    model_data <- model_data %>%
        select(all_of(c(features, target))) %>%
        drop_na()
    
    # Skip if not enough data --> make a warning if there are too little datasets
    # nrow(model_data) < 10
    
    # Split data
    set.seed(123)
    train_index <- createDataPartition(model_data[[target]], p = 0.8, list = FALSE)
    train_data <- model_data[train_index, ]
    test_data <- model_data[-train_index, ]
    
    # Train model
    rf_model <- randomForest(NHISS_Diff ~ ., data = train_data, importance = TRUE, ntree = 500)
    
    
    # Predict
    predictions <- predict(rf_model, newdata = test_data)
    
    # Save prediction plot
    plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
    p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
        geom_point(color = "steelblue", alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
        labs(title = paste("RF Predictions vs Actual - NHISS_Diff"),
             x = "Actual NHISS_Diff", y = "Predicted NHISS_Diff") +
        theme_minimal()
    ggsave(filename = file.path(plot_folder, paste0("RF_Pred_vs_Actual_NHISS_Diff.png")), plot = p)
    
    # Save variable importance plot
    png(filename = file.path(imp_folder, paste0("VarImpPlot_NHISS_Diff.png")), width = 800, height = 600)
    varImpPlot(rf_model, main = paste("Variable Importance - NHISS_Diff"))
    dev.off()
    
    # Extract importance
    imp <- importance(rf_model)[, 1]  # MeanDecreaseAccuracy
    imp_df <- data.frame(t(imp))

# Save to CSV
write.csv(imp_df, "RF_Feature_Importances_Pivot_NHISS_Diff.csv", row.names = FALSE)

imp_df <- data.frame(Feature = names(imp), Importance = as.vector(imp))

imp_df <- imp_df %>% arrange(desc(Importance))


imp_df <- data.frame(Feature = rownames(importance(rf_model)), Importance = importance(rf_model)[, 1])
p_imp <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Variable Importance", x = "Feature", y = "Mean Decrease in Accuracy") +
    theme_minimal()
ggsave(filename = file.path(imp_folder, "VarImpPlot_NHISS_Diff_ggplot.png"), plot = p_imp)

png(filename = file.path(imp_folder, paste0("VarImpPlot_NHISS_Diff_compared.png")), width = 800, height = 600)
varImpPlot(rf_model, main = paste("Variable Importance compared - NHISS_Diff"))
dev.off()

### for NHISS_End instead of Diff -----------------
target <- "NHISS_End"

model_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")

model_data <- model_data %>%
    select(all_of(c(features, target))) %>%
    drop_na()

# Skip if not enough data --> make a warning if there are too little datasets
# nrow(model_data) < 10

# Split data
set.seed(123)
train_index <- createDataPartition(model_data[[target]], p = 0.8, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Train model
rf_model <- randomForest(NHISS_End ~ ., data = train_data, importance = TRUE, ntree = 500)


# Predict
predictions <- predict(rf_model, newdata = test_data)

# Save prediction plot
plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
    geom_point(color = "steelblue", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
    labs(title = paste("RF Predictions vs Actual - NHISS_End"),
         x = "Actual NHISS_End", y = "Predicted NHISS_End") +
    theme_minimal()
ggsave(filename = file.path(plot_folder, paste0("RF_Pred_vs_Actual_NHISS_End.png")), plot = p)

# Save variable importance plot
png(filename = file.path(imp_folder, paste0("VarImpPlotNHISS_End.png")), width = 800, height = 600)
varImpPlot(rf_model, main = paste("Variable Importance"))
dev.off()

# Extract importance
imp <- importance(rf_model)[, 1]  # MeanDecreaseAccuracy
imp_df <- data.frame(t(imp))

# Save to CSV
write.csv(imp_df, "RF_Feature_Importances_Pivot_NHISS_End.csv", row.names = FALSE)

imp_df <- data.frame(Feature = names(imp), Importance = as.vector(imp))

imp_df <- imp_df %>% arrange(desc(Importance))


imp_df <- data.frame(Feature = rownames(importance(rf_model)), Importance = importance(rf_model)[, 1])
p_imp <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Variable Importance", x = "Feature", y = "Mean Decrease in Accuracy") +
    theme_minimal()
ggsave(filename = file.path(imp_folder, "VarImpPlot_NHISS_End_ggplot.png"), plot = p_imp)

png(filename = file.path(imp_folder, paste0("VarImpPlot_NHISS_End_compared.png")), width = 800, height = 600)
varImpPlot(rf_model, main = paste("Variable Importance compared - NHISS_End"))
dev.off()

# # Filter for TP2
# tp2_data <- dataFINALmean %>% filter(Timepoint == "TP2")
# 
# # Define features and target
# features <- c("Subpopulation", "Age", "Sex", "Category", "mean_capZ",
#               "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression")
# #target <- "NHISS_Ratio"
# target <- "NHISS_Diff"
# 
# # Initialize list to store importance
# importance_list <- list()
# 
# # Loop through each gene
# for (gene in unique(tp2_data$Gene)) {
#     gene_data <- tp2_data %>% filter(Gene == gene)
#     
#     # Prepare data
#     model_data <- gene_data %>%
#         select(all_of(c(features, target))) %>%
#         na.omit()
#     
#     # Convert categorical variables to factors
#     model_data$Sex <- as.factor(model_data$Sex)
#     model_data$Category <- as.factor(model_data$Category)
#     model_data$Subpopulation <- as.factor(model_data$Subpopulation)
#     
#     # Skip if not enough data
#     if (nrow(model_data) < 10) next
#     
#     # Split data
#     set.seed(123)
#     train_index <- createDataPartition(model_data[[target]], p = 0.8, list = FALSE)
#     train_data <- model_data[train_index, ]
#     test_data <- model_data[-train_index, ]
#     
#     # Train model
#     rf_model <- randomForest(NHISS_Diff ~ ., data = train_data, importance = TRUE, ntree = 500)
#     
#     
#     # Predict
#     predictions <- predict(rf_model, newdata = test_data)
#     
#     # Save prediction plot
#     plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
#     p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
#         geom_point(color = "steelblue", alpha = 0.7) +
#         geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
#         labs(title = paste("RF Predictions vs Actual for", gene),
#              x = "Actual NHISS_Ratio", y = "Predicted NHISS_Ratio") +
#         theme_minimal()
#     ggsave(filename = file.path(plot_folder, paste0("RF_Pred_vs_Actual_", gene, ".png")), plot = p)
#     
#     # Save variable importance plot
#     png(filename = file.path(imp_folder, paste0("VarImpPlot_", gene, ".png")), width = 800, height = 600)
#     varImpPlot(rf_model, main = paste("Variable Importance for", gene))
#     dev.off()
#     
#     # Extract importance
#     imp <- importance(rf_model)[, 1]  # MeanDecreaseAccuracy
#     imp_df <- data.frame(Gene = gene, t(imp))
#     importance_list[[gene]] <- imp_df
# }
# 
# # Combine all into one data frame
# importance_matrix <- bind_rows(importance_list)
# 
# # Save to CSV
# write.csv(importance_matrix, "RF_Feature_Importances_Pivot.csv", row.names = FALSE)
