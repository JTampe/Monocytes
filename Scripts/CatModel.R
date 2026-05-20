# Install if not already installed
# install.packages("devtools")
# devtools::install_github('catboost/catboost', subdir = 'catboost/R-package')

library(catboost)
library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)

# Define target
target <- "NHISS_Diff"

# Filter and clean
model_data <- matrix_meanCapZ %>%
    filter(Timepoint == "TP2") %>%
    select(all_of(c(features, target))) %>%
    drop_na()

# Convert categorical variables
model_data$Sex <- as.factor(model_data$Sex)
model_data$Subpopulation <- as.factor(model_data$Subpopulation)

# Split data
set.seed(123)
train_index <- createDataPartition(model_data[[target]], p = 0.8, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Identify categorical features
cat_features <- which(sapply(train_data, is.factor))

# Convert to Pool format
train_pool <- catboost.load_pool(data = train_data[, -which(names(train_data) == target)],
                                 label = train_data[[target]],
                                 cat_features = cat_features)

test_pool <- catboost.load_pool(data = test_data[, -which(names(test_data) == target)],
                                label = test_data[[target]],
                                cat_features = cat_features)

# Train model
cat_model <- catboost.train(train_pool,
                            params = list(loss_function = 'RMSE',
                                          iterations = 500,
                                          learning_rate = 0.05,
                                          depth = 6,
                                          eval_metric = 'RMSE',
                                          random_seed = 123,
                                          verbose = 50))

# Predict
predictions <- catboost.predict(cat_model, test_pool)

# Plot predictions
plot_data <- data.frame(Actual = test_data[[target]], Predicted = predictions)
p <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
    geom_point(color = "steelblue", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
    labs(title = "CatBoost Predictions vs Actual - NHISS_Diff",
         x = "Actual", y = "Predicted") +
    theme_minimal()
ggsave("CatBoost_Pred_vs_Actual_NHISS_Diff.png", plot = p)


# Get importance
importance <- catboost.get_feature_importance(cat_model, pool = train_pool)
importance_df <- data.frame(Feature = colnames(train_data[, -which(names(train_data) == target)]),
                            Importance = importance)

# Plot
p_imp <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "CatBoost Feature Importance", x = "Feature", y = "Importance") +
    theme_minimal()
ggsave("CatBoost_Feature_Importance.png", plot = p_imp)

