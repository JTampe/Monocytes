# *Demographics of Ctr & Patients* -----------------------------------------------------------------------------

# for FACS
metadataP_FACS <- metadataP %>% filter(!SampleID %in% Unmatched_TP0_FACS)
Demographics_FACS <- demographics_N_Age_Sex(metadataP_FACS)
write.csv(Demographics_FACS, file = "Demographics_FACS.csv", row.names = FALSE)
# for Gene expr.
metadataP_Gene <- metadataP %>% filter(SampleID %in% unique(data_mean_matched$SampleID))
Demographics_Gene <- demographics_N_Age_Sex(metadataP_Gene)
write.csv(Demographics_Gene, file = "Demographics_Gene.csv", row.names = FALSE)

