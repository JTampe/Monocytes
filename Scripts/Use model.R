
selected_row <- rf_data %>%
    drop_na() %>%
    slice_sample(n = 1)

input_row <- selected_row %>%
    select(all_of(features)) 

#input_row

predict(rf_model, newdata = input_row)

selected_row$NHISS_Diff
