# *Sig Regression - Gene Overview* ----------------------------------------------------
# with a mean for all Timepoints
Gene_Overview_NHISS <- LinReg_NHISS_capZ %>%
    filter(p.value <= 0.05) %>% 
    filter(Subgroup == "All") %>% 
    group_by(Subpopulation, Gene)%>%
    summarise(NHISS_Correlation = mean(Estimate, na.rm = TRUE))

Gene_Overview_NHISS_End <- LinReg_NHISS_End_capZ %>%
    filter(p.value <= 0.05) %>%  
    filter(Subgroup == "All") %>%
    group_by(Subpopulation, Gene)%>%
    summarise(NHISS_End_Correlation = mean(Estimate, na.rm = TRUE))

Gene_Overview_NHISS_Diff <- LinReg_NHISS_Diff_capZ %>%
    filter(p.value <= 0.05) %>%  
    filter(Subgroup == "All") %>%
    group_by(Subpopulation, Gene)%>%
    summarise(NHISS_Diff_Correlation = mean(Estimate, na.rm = TRUE))

Gene_Overview_NHISS_Ratio <- LinReg_NHISS_Ratio_capZ %>%
    filter(p.value <= 0.05) %>%  
    filter(Subgroup == "All") %>%
    group_by(Subpopulation, Gene)%>%
    summarise(NHISS_Ratio_Correlation = mean(Estimate, na.rm = TRUE))

# Gene_Overview_Age <- Age_correlation_capZ %>%
#     filter(p.value <= 0.05) %>% 
#     group_by(Subpopulation, Gene)%>%
#     summarise(Age_Correlation = mean(Estimate, na.rm = TRUE))
# Gene_Overview_SPAN <- LinReg_SPAN_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene)%>%
#     summarise(SPAN_Correlation = mean(Estimate, na.rm = TRUE))
# 
# Gene_Overview_mRS <- LinReg_mRS_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene)%>%
#     summarise(mRS_Correlation = mean(Estimate, na.rm = TRUE))
# 
# Gene_Overview_Barthel <- LinReg_Barthel_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene)%>%
#     summarise(Barthel_Correlation = mean(Estimate, na.rm = TRUE))
# 
# Gene_Overview_MoCA <- LinReg_MoCA_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene)%>%
#     summarise(MoCA_Correlation = mean(Estimate, na.rm = TRUE))
# 
# Gene_Overview_HADS_Anxiety <- LinReg_HADS_Anxiety_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene)%>%
#     summarise(HADS_Anxiety_Correlation = mean(Estimate, na.rm = TRUE))
# 
# Gene_Overview_HADS_Depression <- LinReg_HADS_Depression_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene)%>%
#     summarise(HADS_Depression_Correlation = mean(Estimate, na.rm = TRUE))

Gene_Overview_woTP <- merge(Gene_Overview_NHISS,
                            Gene_Overview_NHISS_End,
                            by = c("Gene", "Subpopulation"), 
                            all = TRUE)
Gene_Overview_woTP <- merge(Gene_Overview_woTP,
                            Gene_Overview_NHISS_Diff, 
                            by = c("Gene", "Subpopulation"), 
                            all = TRUE)
Gene_Overview_woTP <- merge(Gene_Overview_woTP,
                            Gene_Overview_NHISS_Ratio, 
                            by = c("Gene", "Subpopulation"), 
                            all = TRUE)

# Gene_Overview_woTP <- merge(Gene_Overview_woTP,
#                             Gene_Overview_mRS, 
#                             by = c("Gene", "Subpopulation"), 
#                             all = TRUE)
# Gene_Overview_woTP <- merge(Gene_Overview_woTP,
#                             Gene_Overview_Barthel, 
#                             by = c("Gene", "Subpopulation"), 
#                             all = TRUE)
# Gene_Overview_woTP <- merge(Gene_Overview_woTP,
#                             Gene_Overview_MoCA, 
#                             by = c("Gene", "Subpopulation"), 
#                             all = TRUE)
# Gene_Overview_woTP <- merge(Gene_Overview_woTP,
#                             Gene_Overview_HADS_Anxiety, 
#                             by = c("Gene", "Subpopulation"), 
#                             all = TRUE)
# Gene_Overview_woTP <- merge(Gene_Overview_woTP,
#                             Gene_Overview_HADS_Depression, 
#                             by = c("Gene", "Subpopulation"), 
#                             all = TRUE)
# # Aging at the end!
# Gene_Overview_woTP <- merge(Gene_Overview_woTP,
#                             Gene_Overview_Age,
#                             by = c("Gene", "Subpopulation"), 
#                             all = TRUE)
# Gene_Overview_woTP <- merge(Gene_Overview_woTP,
#                             Gene_Overview_SPAN, 
#                             by = c("Gene", "Subpopulation"), 
#                             all = TRUE)

# table including the details for the timepoints

Gene_Overview_NHISS <- LinReg_NHISS_capZ %>%
    filter(p.value <= 0.05) %>%  
    filter(Subgroup == "All") %>%
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarise(NHISS_Correlation = Estimate)
Gene_Overview_NHISS_End <- LinReg_NHISS_End_capZ %>%
    filter(p.value <= 0.05) %>%  
    filter(Subgroup == "All") %>%
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarise(NHISS_End_Correlation = Estimate)
Gene_Overview_NHISS_Diff <- LinReg_NHISS_Diff_capZ %>%
    filter(p.value <= 0.05) %>%  
    filter(Subgroup == "All") %>%
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarise(NHISS_Diff_Correlation = Estimate)
Gene_Overview_NHISS_Ratio <- LinReg_NHISS_Ratio_capZ %>%
    filter(p.value <= 0.05) %>%  
    filter(Subgroup == "All") %>%
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarise(NHISS_Ratio_Correlation = Estimate)

# Gene_Overview_Age <- Age_correlation_capZ %>%
#     filter(p.value <= 0.05) %>% 
#     group_by(Subpopulation, Gene, Timepoint) %>%
#     summarise(Age_Correlation = Estimate)
# Gene_Overview_SPAN <- LinReg_SPAN_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene, Timepoint) %>%
#     summarise(SPAN_Correlation = Estimate)
# 
# Gene_Overview_mRS <- LinReg_mRS_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene, Timepoint) %>%
#     summarise(mRS_Correlation = Estimate)
# Gene_Overview_Barthel <- LinReg_Barthel_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene, Timepoint) %>%
#     summarise(Barthel_Correlation = Estimate)
# Gene_Overview_MoCA <- LinReg_MoCA_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene, Timepoint) %>%
#     summarise(MoCA_Correlation = Estimate)
# Gene_Overview_HADS_Anxiety <- LinReg_HADS_Anxiety_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene, Timepoint) %>%
#     summarise(HADS_Anxiety_Correlation = Estimate)
# Gene_Overview_HADS_Depression <- LinReg_HADS_Depression_capZ %>%
#     filter(p.value <= 0.05) %>%  
#     filter(Subgroup == "All") %>%
#     group_by(Subpopulation, Gene, Timepoint) %>%
#     summarise(HADS_Depression_Correlation = Estimate)


Gene_Overview <- merge(Gene_Overview_NHISS,
                       Gene_Overview_NHISS_End, 
                       by = c("Gene", "Subpopulation", "Timepoint"), 
                       all = TRUE)
Gene_Overview <- merge(Gene_Overview,
                       Gene_Overview_NHISS_Diff, 
                       by = c("Gene", "Subpopulation", "Timepoint"), 
                       all = TRUE)
Gene_Overview <- merge(Gene_Overview,
                       Gene_Overview_NHISS_Ratio, 
                       by = c("Gene", "Subpopulation", "Timepoint"), 
                       all = TRUE)
# Gene_Overview <- merge(Gene_Overview,
#                        Gene_Overview_mRS, 
#                        by = c("Gene", "Subpopulation", "Timepoint"), 
#                        all = TRUE)
# Gene_Overview <- merge(Gene_Overview,
#                        Gene_Overview_Barthel, 
#                        by = c("Gene", "Subpopulation", "Timepoint"), 
#                        all = TRUE)
# Gene_Overview <- merge(Gene_Overview,
#                        Gene_Overview_MoCA, 
#                        by = c("Gene", "Subpopulation", "Timepoint"), 
#                        all = TRUE)
# Gene_Overview <- merge(Gene_Overview,
#                        Gene_Overview_HADS_Anxiety, 
#                        by = c("Gene", "Subpopulation", "Timepoint"), 
#                        all = TRUE)
# Gene_Overview <- merge(Gene_Overview,
#                        Gene_Overview_HADS_Depression, 
#                        by = c("Gene", "Subpopulation", "Timepoint"), 
#                        all = TRUE)
# # put aging at the end
# Gene_Overview <- merge(Gene_Overview,
#                        Gene_Overview_Age, 
#                        by = c("Gene", "Subpopulation", "Timepoint"), 
#                        all = TRUE)
# Gene_Overview <- merge(Gene_Overview,
#                        Gene_Overview_SPAN, 
#                        by = c("Gene", "Subpopulation", "Timepoint"), 
#                        all = TRUE)

# Replace the Estimated with the arrows!    
Gene_Overview_woTP_arrows <- Gene_Overview_woTP %>%
    mutate(across(c(NHISS_Correlation, NHISS_End_Correlation, NHISS_Diff_Correlation, NHISS_Ratio_Correlation, 
                    #Age_Correlation, SPAN_Correlation, mRS_Correlation, Barthel_Correlation, MoCA_Correlation, HADS_Anxiety_Correlation, HADS_Depression_Correlation, 
    ), 
    ~ ifelse(is.na(.), "-", ifelse(. > 0, "↑", "↓"))))
Gene_Overview_arrows <- Gene_Overview %>%
    mutate(across(c(NHISS_Correlation, NHISS_End_Correlation, NHISS_Diff_Correlation, NHISS_Ratio_Correlation, 
                    #Age_Correlation, SPAN_Correlation, mRS_Correlation, Barthel_Correlation, MoCA_Correlation, HADS_Anxiety_Correlation, HADS_Depression_Correlation, 
    ), 
    ~ ifelse(is.na(.), "-", ifelse(. > 0, "↑", "↓"))))

# Save the files:
file_name_s <- "LinReg_Results_capZ/Gene_Overview_woTP.csv"
write.csv(Gene_Overview_woTP, file_name_s, row.names = FALSE)
file_name_s <- "LinReg_Results_capZ/Gene_Overview.csv"
write.csv(Gene_Overview, file_name_s, row.names = FALSE)

file_name_s_arrows <- "LinReg_Results_capZ/Gene_Overview_woTP_arrows.csv"
write.csv(Gene_Overview_woTP_arrows, file_name_s_arrows, row.names = FALSE)
file_name_s_arrows <- "LinReg_Results_capZ/Gene_Overview_arrows.csv"
write.csv(Gene_Overview_arrows, file_name_s_arrows, row.names = FALSE)
