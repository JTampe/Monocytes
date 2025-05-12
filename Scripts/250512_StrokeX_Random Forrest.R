# Load libraries
library(dplyr)
library(randomForest)
library(caret)
library(ggplot2)

# Load and filter FACS data
facs <- read.csv("FACSdata_clean.csv", stringsAsFactors = TRUE)
tp2_data <- facs %>% filter(Timepoint == "TP2")

# Define features and target
features <- c("Monocytes", "Classical_Monocytes", "Intermediate_Monocytes", "NonClassical_Monocytes",
              "Age", "Sex", "Category", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression")
target <- "NHISS_Ratio"

# Prepare data and drop rows with any missing values
model_data <- tp2_data %>%
    select(all_of(c(features, target))) %>%
    na.omit()

# Convert categorical variables to factors
model_data$Sex <- as.factor(model_data$Sex)
model_data$Category <- as.factor(model_data$Category)

# Split into training and test sets
set.seed(123)
train_index <- createDataPartition(model_data[[target]], p = 0.8, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Train Random Forest model
rf_model <- randomForest(NHISS_Ratio ~ ., data = train_data, importance = TRUE, ntree = 500)

# Predict
predictions <- predict(rf_model, newdata = test_data)

# Plot predicted vs actual
ggplot(data = NULL, aes(x = test_data[[target]], y = predictions)) +
    geom_point(color = "steelblue", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
    labs(title = "Random Forest Predictions vs Actual (TP2)",
         x = "Actual NHISS_Ratio",
         y = "Predicted NHISS_Ratio") +
    theme_minimal()

# Optional: View variable importance
importance(rf_model)
varImpPlot(rf_model)


#### Random Forrest TP2 for gene expression data
# Filter data
tp2_data <- dataFINALmean %>% filter(Timepoint == "TP2")

# Define features and target
features <- c("SampleID", "Subpopulation","Age", "Sex", "Category", "mean_capZ",
              "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression")
target <- "NHISS_Ratio"

# Prepare data and drop rows with any missing values
model_data <- tp2_data %>%
    select(all_of(c(features, target))) %>%
    na.omit()

# Convert categorical variables to factors
model_data$Sex <- as.factor(model_data$Sex)
model_data$Category <- as.factor(model_data$Category)

# Split into training and test sets
set.seed(123)
train_index <- createDataPartition(model_data[[target]], p = 0.8, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Train Random Forest model
rf_model <- randomForest(NHISS_Ratio ~ ., data = train_data, importance = TRUE, ntree = 500)

# Predict
predictions <- predict(rf_model, newdata = test_data)

# Plot predicted vs actual
ggplot(data = NULL, aes(x = test_data[[target]], y = predictions)) +
    geom_point(color = "steelblue", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
    labs(title = "Random Forest Predictions vs Actual (TP2)",
         x = "Actual NHISS_Ratio",
         y = "Predicted NHISS_Ratio") +
    theme_minimal()

# Optional: View variable importance
importance(rf_model)
varImpPlot(rf_model)
