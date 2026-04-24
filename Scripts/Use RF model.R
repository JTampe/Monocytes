rd_data <- matrix_meanCapZ %>% filter(Timepoint == "TP2")

selected_row <- rf_data %>%
    drop_na() %>%
    slice_sample(n = 1)

input_row <- selected_row %>%
    select(all_of(features)) 

#input_row

predict(rf_model, newdata = input_row)

selected_row$NHISS_Diff

run_rf_cv_analysis <- function(data, features, target, prefix = "RF_CV_Model", plot_folder = "RF_Plots", imp_folder = "RF_Importance_Plots")
    
predict(rf_model, newdata = input_row)