# *Demographics of Ctr & Patients* -----------------------------------------------------------------------------

# for FACS
metadataP_FACS <- metadataP %>% filter(!SampleID %in% Unmatched_TP0_FACS)
Demographics_FACS <- demographics_N_Age_Sex(metadataP_FACS)
write.csv(Demographics_FACS, file = "Demographics_FACS.csv", row.names = FALSE)
# for Gene expr.
metadataP_Gene <- metadataP %>% filter(SampleID %in% unique(data_mean_matched$SampleID))
Demographics_Gene <- demographics_N_Age_Sex(metadataP_Gene)
write.csv(Demographics_Gene, file = "Demographics_Gene.csv", row.names = FALSE)

#### NHISS Recovery----------------------------------------------------------------------------------------------------

NHISS_Recovery <- merge(Metadata_NeuroTest, metadataP, by = "SampleID", all.x = TRUE)

# Remove Patient 15
NHISS_Recovery <- NHISS_Recovery %>% filter(!(SampleID == "Patient15"))
NHISS_Recovery$Timepoint <- factor(NHISS_Recovery$Timepoint, levels = c("TP0", "TP1", "TP2", "TP3", "TP4", "TP5"))

ggline(NHISS_Recovery, 
       x = "Timepoint", 
       y = "NHISS", 
       color = "Recovery", 
       add = "mean_se",  # Add mean and standard error
       size = 1.2) +  
    labs(title = paste("Mean ± SEM of NHISS by Recovery across Timepoints")) +
    scale_color_manual(values = c("Good" = "darkgreen", "Bad" = "#700606")) +
    theme_bw() +
    scale_x_discrete(labels = c(
        "TP0" = "Control",
        "TP1" = "24 hours",
        "TP2" = "3-5 days",
        "TP3" = "1 month",
        "TP4" = "3 months",
        "TP5" = ">1 year"
    )) +
    theme(
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black")
    ) +
    stat_compare_means(method = "wilcox.test",  # or "t.test"
                       aes(group = Recovery),
                       label = "p.signif",  
                       size = 5) 

ggsave(filename = "Patient_NHISS_Recovery.png")
