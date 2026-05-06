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

gene_panel <- c("CD33", "CD169", "TNFA", "TGFB1", "CD38")


result_important <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2"),
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", 
                 as.vector(gene_panel)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_panel"
)

## split by Subpopulation
# all
result_important <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2" & Subpopulation == "all"),
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", 
                 as.vector(genes_tp2)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_all"
)
# all
result_important <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2" & Subpopulation == "all"),
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", 
                 as.vector(genes_tp2)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_all"
)
# classical
result_important <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2" & Subpopulation == "classical"),
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", 
                 as.vector(genes_tp2)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_classical"
)
# intermediate
result_important <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2" & Subpopulation == "intermediate"),
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", 
                 as.vector(genes_tp2)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_intermediate"
)
# nonclassical
result_important <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2" & Subpopulation == "nonclassical"),
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", 
                 as.vector(genes_tp2)),
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_nonclassical"
)

#### single features:
result_important <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2"),
    features = "NHISS",
    target = "NHISS_Diff",
    prefix = "NHISS_Diff_NHISS"
)


### for NHISS_End instead of Diff ----------------------------------
result_End <- run_rf_analysis(
    data = matrix_meanCapZ %>% filter(Timepoint == "TP2"),
    features = c("Subpopulation", "Age", "Sex", "NHISS", "mRS", "Barthel", "MoCA", "HADS_Anxiety", "HADS_Depression", as.vector(genes_tp2)),
    target = "NHISS_End",
    prefix = "NHISS_End"
)
