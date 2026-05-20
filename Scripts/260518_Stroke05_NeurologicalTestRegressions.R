# *Neurological Tests Regressions* --------------------------------------------------------
#### NIHSS Regression --------------------------------------------------------
LinReg_NIHSS_capZ <- perform_linear_regression_correlation(dataFINALmean, "NIHSS", "NIHSS")

LinReg_NIHSS_sigGenes <- LinReg_NIHSS_capZ %>% filter(LinReg_NIHSS_capZ$p.value <0.05)
write.csv(LinReg_NIHSS_sigGenes, "LinReg_Results_capZ/LinReg_NIHSS_sigGenes.csv", row.names = FALSE)

#### NIHSS_End Regression* --------------------------------------------------------
# Correlation with "final" NIHSS
dataFINALmean$NIHSS_End <- NA
# Assign values from NIHSS where Timepoint is TP4
dataFINALmean$NIHSS_End[dataFINALmean$Timepoint == "TP4"] <- dataFINALmean$NIHSS[dataFINALmean$Timepoint == "TP4"]
dataFINALmean$NIHSS_End[dataFINALmean$Timepoint == "TP3"] <- dataFINALmean$NIHSS[dataFINALmean$Timepoint == "TP4"]
dataFINALmean$NIHSS_End[dataFINALmean$Timepoint == "TP2"] <- dataFINALmean$NIHSS[dataFINALmean$Timepoint == "TP4"]
dataFINALmean$NIHSS_End[dataFINALmean$Timepoint == "TP1"] <- dataFINALmean$NIHSS[dataFINALmean$Timepoint == "TP4"]

LinReg_NIHSS_End_capZ <- perform_linear_regression_correlation(dataFINALmean, "NIHSS_End", "NIHSS_End")

#Focus only on TP1 & TP2
LinReg_NIHSS_End_capZ <-LinReg_NIHSS_End_capZ %>% filter(LinReg_NIHSS_End_capZ$Timepoint == c("TP1", "TP2"))

LinReg_NIHSS_End_sigGenes <- LinReg_NIHSS_End_capZ %>% filter(LinReg_NIHSS_End_capZ$p.value <0.05)
write.csv(LinReg_NIHSS_End_sigGenes, "LinReg_Results_capZ/LinReg_NIHSS_End_sigGenes.csv", row.names = FALSE)

#### NIHSS_Diff Regression* --------------------------------------------------------
# Create a new column NIHSS_Diff initialized to NA
dataFINALmean$NIHSS_Diff <- NA

# Calculate the difference for each SampleID
dataFINALmean <- dataFINALmean %>%
    group_by(SampleID, Subpopulation, Gene) %>%
    mutate(
        NIHSS_Diff = if (all(c("TP1", "TP4") %in% Timepoint)) {
            NIHSS[Timepoint == "TP1"] - NIHSS[Timepoint == "TP4"]
        } else {
            NA
        }
    ) %>%
    ungroup()

LinReg_NIHSS_Diff_capZ <- perform_linear_regression_correlation(dataFINALmean, "NIHSS_Diff", "NIHSS_Diff")

#Focus only on TP1 & TP2
LinReg_NIHSS_Diff_capZ <-LinReg_NIHSS_Diff_capZ %>% filter(LinReg_NIHSS_Diff_capZ$Timepoint == c("TP1", "TP2"))

LinReg_NIHSS_Diff_sigGenes <- LinReg_NIHSS_Diff_capZ %>% filter(LinReg_NIHSS_Diff_capZ$p.value <0.05)
write.csv(LinReg_NIHSS_Diff_sigGenes, "LinReg_Results_capZ/LinReg_NIHSS_Diff_sigGenes.csv", row.names = FALSE)

#### NIHSS_Ratio Regression* --------------------------------------------------------
# Create a new column NIHSS_Diff initialized to NA
dataFINALmean$NIHSS_Ratio <- NA

# Calculate the ratio for each SampleID
dataFINALmean <- dataFINALmean %>%
    group_by(SampleID, Subpopulation, Gene) %>%
    mutate(
        NIHSS_Ratio = if (all(c("TP1", "TP4") %in% Timepoint)) {
            (NIHSS[Timepoint == "TP1"] - NIHSS[Timepoint == "TP4"])/NIHSS[Timepoint == "TP1"]
        } else {
            NA
        }
    ) %>%
    ungroup()

LinReg_NIHSS_Ratio_capZ <- perform_linear_regression_correlation(dataFINALmean, "NIHSS_Ratio", "NIHSS_Ratio")

#Focus only on TP1 & TP2
LinReg_NIHSS_Ratio_capZ <-LinReg_NIHSS_Ratio_capZ %>% filter(LinReg_NIHSS_Ratio_capZ$Timepoint == c("TP1", "TP2"))
#write.csv(LinReg_NIHSS_Ratio_capZ, "LinReg_Results_capZ/LinReg_NIHSS_Ratio_capZ.csv", row.names = FALSE)

LinReg_NIHSS_Ratio_sigGenes <- LinReg_NIHSS_Ratio_capZ %>% filter(LinReg_NIHSS_Ratio_capZ$p.value <0.05)
write.csv(LinReg_NIHSS_Ratio_sigGenes, "LinReg_Results_capZ/LinReg_NIHSS_Ratio_sigGenes.csv", row.names = FALSE)


# #### mRS Regression --------------------------------------------------------
# LinReg_mRS_capZ <- perform_linear_regression_correlation(dataFINALmean, "mRS", "mRS")
# 
# LinReg_mRS_sigGenes <- LinReg_mRS_capZ %>% filter(LinReg_mRS_capZ$p.value <0.05)
# write.csv(LinReg_mRS_sigGenes, "LinReg_Results_capZ/LinReg_mRS_sigGenes.csv", row.names = FALSE)
# 
# #### Barthel Regression --------------------------------------------------------
# LinReg_Barthel_capZ <- perform_linear_regression_correlation(dataFINALmean, "Barthel", "Barthel")
# 
# LinReg_Barthel_sigGenes <- LinReg_Barthel_capZ %>% filter(LinReg_Barthel_capZ$p.value <0.05)
# write.csv(LinReg_Barthel_sigGenes, "LinReg_Results_capZ/LinReg_Barthel_sigGenes.csv", row.names = FALSE)
# 
# #### MoCA Regression --------------------------------------------------------
# LinReg_MoCA_capZ <- perform_linear_regression_correlation(dataFINALmean, "MoCA", "MoCA")
# 
# LinReg_MoCA_sigGenes <- LinReg_MoCA_capZ %>% filter(LinReg_MoCA_capZ$p.value <0.05)
# write.csv(LinReg_MoCA_sigGenes, "LinReg_Results_capZ/LinReg_MoCA_sigGenes.csv", row.names = FALSE)
# 
# #### HADS_Anxiety Regression --------------------------------------------------------
# LinReg_HADS_Anxiety_capZ <- perform_linear_regression_correlation(dataFINALmean, "HADS_Anxiety", "HADS_Anxiety")
# 
# LinReg_HADS_Anxiety_sigGenes <- LinReg_HADS_Anxiety_capZ %>% filter(LinReg_HADS_Anxiety_capZ$p.value <0.05)
# write.csv(LinReg_HADS_Anxiety_sigGenes, "LinReg_Results_capZ/LinReg_HADS_Anxiety_sigGenes.csv", row.names = FALSE)
# 
# #### HADS_Depression Regression --------------------------------------------------------
# LinReg_HADS_Depression_capZ <- perform_linear_regression_correlation(dataFINALmean, "HADS_Depression", "HADS_Depression")
# 
# LinReg_HADS_Depression_sigGenes <- LinReg_HADS_Depression_capZ %>% filter(LinReg_HADS_Depression_capZ$p.value <0.05)
# write.csv(LinReg_HADS_Depression_sigGenes, "LinReg_Results_capZ/LinReg_HADS_Depression_sigGenes.csv", row.names = FALSE)
# 
# #### SPAN Regression --------------------------
# 
# # Create the SPAN variable as per Almekhlafi et. al 2014:
# dataFINALmean <- dataFINALmean %>%
#     mutate(SPAN = Age + NIHSS)
# 
# LinReg_SPAN_capZ <- perform_linear_regression_correlation(dataFINALmean, "SPAN", "SPAN")
# 
# LinReg_SPAN_sigGenes <- LinReg_SPAN_capZ %>% filter(LinReg_SPAN_capZ$p.value <0.05)
# write.csv(LinReg_SPAN_sigGenes, "LinReg_Results_capZ/LinReg_SPAN_sigGenes.csv", row.names = FALSE)
