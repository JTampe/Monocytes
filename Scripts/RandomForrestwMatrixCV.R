run_rf_cv <- function(data, features, target, prefix = "RF_CV_Model", plot_folder = "RF_Plots", imp_folder = "RF_Importance_Plots") {
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
    
    # Set up cross-validation and tuning grid
    trctrl <- trainControl(method = "cv", number = 5)
    tunegrid <- expand.grid(.mtry = seq(2, floor(sqrt(length(features))), by = 1))
    
    # Train model with cross-validation and tuning
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
    
    # Predict on test set
    predictions <- predict(rf_model, newdata = test_data[, features])
    
    # Save prediction plot
    plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
    p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
        geom_point(color = "steelblue", alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
        labs(title = paste("RF CV Predictions vs Actual -", prefix),
             x = paste("Actual", target), y = paste("Predicted", target)) +
        theme_minimal()
    ggsave(filename = file.path(plot_folder, paste0("RF_CV_Pred_vs_Actual_", prefix, ".png")), plot = p)
    
    # Extract importance (from final model)
    var_imp <- varImp(rf_model)$importance
    imp_df <- data.frame(Feature = rownames(var_imp), Importance = var_imp[,1]) %>%
        arrange(desc(Importance))
    
    # Save importance table
    write.csv(imp_df, file.path(imp_folder, paste0("RF_CV_Feature_Importances_", prefix, ".csv")), row.names = FALSE)
    
    # Save ggplot version of importance
    p_imp <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = paste("Variable Importance (CV) -", prefix), x = "Feature", y = "Importance") +
        theme_minimal()
    ggsave(filename = file.path(imp_folder, paste0("VarImpPlot_CV_", prefix, "_ggplot.png")), plot = p_imp)
    
    return(list(rf_model = rf_model, imp_df = imp_df))
}

rf_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")
# z-score age
#rf_data$Age <- scale(rf_data$Age)

working_genes <- setdiff(gene_names, exclude_genes)

all_features <- c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", #"SPAN",
                  working_genes)
                 # "CD33", "CD169", "TNFA", "TGFB1", "CD91", "STAT6")
features_to_remove <- c(
    "Barthel", "B2M", "SLC24A4", "CD32", "CD31", "TREM1", "TGM2",
    "NHISS_Ratio", "IL1B", "TGFB1", "mRS", "CD40", "IL10", "CD86",
    "ACTB", "ANXA1", "NHISS", "NHISS_End", "TSPO", "CX3CR1", "CD11C"
)

features_to_include <- c("BDNF", "CCR2", "CD11B", 
    "CD169", "CD33", "CD36", "CD38", "CD64", "CD91",
    "GAPDH", "IFITM2", "MARCO", "SIGLEC10", "STAT6", "TGFB1", "TLR8", "TNFA",
    "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", "Age"
    #, "NHISS_End", "NHISS_Ratio", "NHISS_Diff", "Sex", "SampleID", "Timepoint", "Subpopulation"
)

gene_features <- c("CD33", "CD169", "TNFA", "TGFB1", "CD91" #, "STAT6"
)

all_features <- setdiff(gene_names, exclude_genes)

result_cv <- run_rf_cv(
    data = rf_data,
    features = features_to_include,
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_1selection"
)

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_clinical"
)

result_cv <- run_rf_cv(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2"),
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_clinical"
)


gene_panel <- c("CD33", "CD169", "TNFA", "TGFB1", "CD91", "STAT6")

result_cv <- run_rf_cv(
    data = rf_data,
    features = c(as.vector(gene_panel)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_topGenes"
)

# wo Subpopulation
result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_X_clinical"
)

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", as.vector(gene_panel)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_X_topGenes"
)

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", "CD33"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_CD33"
)

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", "TGFB1"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_TGFB1"
)


result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", "TNFA"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_TNFA"
)

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", "CD91"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_CD91"
)

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", "CD169"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_CD169"
)

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", "CD169"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_CD169"
)

rf_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", "CD169"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_CD169"
)

result_cv <- run_rf_cv(
    data = rf_data,
    features = c("Age", "NHISS", "CD91"),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_CV_CD169"
)
