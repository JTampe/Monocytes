## Set Up of Environment
# adjust it to all the other packages
if (requireNamespace("thematic")) 
  thematic::thematic_rmd(font = "auto")

today <- Sys.Date()

output_location <- paste(today,"_Stroke_Results", sep="")

setwd("/Users/ju5263ta/Github/Monocytes/Data/")
if (!dir.exists(output_location)) {
    dir.create(output_location, recursive = TRUE)
  }
getwd ()

rm(list=ls())

library(ggplot2)
library(tidyverse)
library(stats)
library(skimr)
library(sjPlot)
library(readxl)
library(thematic)
library(knitr)
library(lme4)
library(ggpubr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(stringr)
library(colorRamp2)
library(rstatix)
library(scales)
#library(emmeans)
library(pwr)
library(lmtest)
library(writexl)
library(openxlsx)
#library(nondetects)
library(skimr)
library(visdat) 
library(mice)
library(lmerTest)
library(beepr)
#library(Seurat)

thematic_on(bg = "white", fg = "black", accent = "blue")

### *Functions* ---------------------------------------------------------------------
my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
my_grey_scale <- c("grey", "black", "black", "black", "black") 

## Create Folder
createFolder <- function(folderName) {
  if (!dir.exists(folderName)) {
    dir.create(folderName, recursive = TRUE)
  }
}

## Split data
createSeparatedDataFiles <- function(data, categoryColumn, categoryFolder) {
  uniqueCategories <- unique(data[[categoryColumn]])
  
  for (category in uniqueCategories) {
    categoryData <- data %>% filter(data[[categoryColumn]] == category)
    categoryName <- as.character(category)
    filename <- file.path(categoryFolder, paste0(categoryName, ".txt"))
    write.table(categoryData, filename, sep = "\t", row.names = FALSE)
    assign(categoryName, categoryData, envir = .GlobalEnv)
  }
  
  return(uniqueCategories)
}

## Function to perform correlation and linear regression by Timepoint
calculate_linreg_by_timepoint <- function(data, cells_colname) {
    
    # Initialize an empty list to store results for each Timepoint
    results_list <- list()
    
    # Loop through each unique Timepoint and perform the correlation and regression
    for (timepoint in unique(data$Timepoint)) {
        
        # Subset data for the current Timepoint
        data_subset <- subset(data, Timepoint == timepoint)
        
        # Pearson correlation
        cor_test <- cor.test(data_subset[[cells_colname]], data_subset$Age)
        
        # Linear regression model
        lm_result <- lm(data_subset[[cells_colname]] ~ data_subset$Age)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        intercept <- summary(lm_result)$coefficients[1, 1]
        
        # Calculate Fold Induction
        Max <- max(data_subset$Age, na.rm = TRUE)
        Min <- min(data_subset$Age, na.rm = TRUE)
        FinalValue <- coef_x * Max + intercept
        InitialValue <- coef_x * Min + intercept
        Fold <- FinalValue / InitialValue
        
        # Store the results in a data frame for this Timepoint
        result <- data.frame(Cells = cells_colname,
                             Timepoint = timepoint,
                             N = length(data_subset[[cells_colname]]),
                             p.value = cor_test$p.value,
                             Coefficient = coef_x,
                             Down95 = ci[1],
                             Up95 = ci[2],
                             Intercept = intercept,
                             Fold = Fold,
                             Max_age = Max,
                             Min_age = Min)
        
        # Append to the results list
        results_list[[timepoint]] <- result
    }
    
    # Combine all results into a single data frame
    final_results <- do.call(rbind, results_list)
    
    return(final_results)
}

getMeTheShitSamples = function( data, gene, cutoff ){
    
    data$Sample[which(data$Gene == gene & data$CT_Value > cutoff
                      #, na.rm = TRUE
    )]
}

# summary data.frame for dotplot
calculate_summary <- function(data, combo_col, gene_col) {
    # Calculate mean_Z
    mean_summary <- data %>%
        group_by(!!sym(combo_col), !!sym(gene_col)) %>%
        summarize(
            mean_Z = mean(mean_Z, na.rm = TRUE),  
            count = n()  
        ) %>%
        filter(count > 1) 
    
    # Calculate sd_mean_Z
    sd_summary <- data %>%
        group_by(!!sym(combo_col), !!sym(gene_col)) %>%
        summarize(
            sd_mean_Z = sd(mean_Z, na.rm = TRUE)  # Calculate standard deviation of mean_Z
        ) 
    
    # Combine both summaries
    final_summary <- mean_summary %>%
        left_join(sd_summary, by = c(combo_col, gene_col)) %>%
        mutate(
            Dot_Size = scales::rescale(sd_mean_Z, to = c(5, 1))  # Rescale SD to range 5 (low SD) to 1 (high SD)
        )
    
    return(final_summary)
}
check_sample_counts <- function(data, timepoint_col = "Timepoint", sampleid_col = "SampleID") {
    # List unique SampleID values for each timepoint
    unique_samples_tp0 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP0"])
    unique_samples_tp1 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP1"])
    unique_samples_tp2 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP2"])
    unique_samples_tp3 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP3"])
    unique_samples_tp4 <- unique(data[[sampleid_col]][data[[timepoint_col]] == "TP4"])
    
    # Store lengths of unique samples for each timepoint
    lengths <- c(
        TP0 = length(unique_samples_tp0),
        TP1 = length(unique_samples_tp1),
        TP2 = length(unique_samples_tp2),
        TP3 = length(unique_samples_tp3),
        TP4 = length(unique_samples_tp4)
    )
    
    # Check if all timepoints have the same number of unique SampleIDs
    all_equal_lengths <- all(lengths == lengths[1])
    
    # Print result of the comparison
    if (all_equal_lengths) {
        cat("All timepoints have the same number of unique SampleIDs.\n")
    } else {
        cat("Not all timepoints have the same number of unique SampleIDs.\n")
    }
    
    # Optionally return the lengths for further inspection
    return(lengths)
}

# get significance stars
get_significance <- function(p_value) {
    ifelse(is.na(p_value), "NA",  # Handle NA values
           ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                         ifelse(p_value < 0.05, "*", "ns"))))
}

# extrTP p-values
get_p_value <- function(tp1, tp2, tukey_comparisons, tukey_result) {
    comparison <- paste(tp1, tp2, sep = "-")
    reverse_comparison <- paste(tp2, tp1, sep = "-")
    
    if (comparison %in% tukey_comparisons) {
        return(tukey_result$Timepoint[comparison, "p adj"])
    } else if (reverse_comparison %in% tukey_comparisons) {
        return(tukey_result$Timepoint[reverse_comparison, "p adj"])
    } else {
        return(NA)
    }
}

my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
# FACS ANOVA
automate_anova_extraction <- function(results_folder, plots_folder, results_name, plot_save, dataset, loop_vars, timepoint_col = "Timepoint") {
    # Create directories if they don't exist
    if (!dir.exists(results_folder)) {
        dir.create(results_folder)
    }
    
    if (!dir.exists(plots_folder)) {
        dir.create(plots_folder)
    }
    
    # Initialize an empty dataframe to store results
    results <- data.frame(Variable = character(), P_Value = numeric(),
                          TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                          TP3_P_Value = numeric(), TP4_P_Value = numeric())
    
    # Loop through the specified variables (columns) in the dataset
    for (var_name in loop_vars) {
        # Subset the dataset for the variable and timepoints
        df_test <- dataset[, c(var_name, timepoint_col)]
        
        # Filter out rows with missing values
        df_test <- df_test %>% filter(!is.na(get(var_name)), !is.na(get(timepoint_col)))
        
        # Ensure Timepoint is treated as a factor
        df_test[[timepoint_col]] <- factor(df_test[[timepoint_col]], levels = c("TP0", "TP1", "TP2", "TP3", "TP4"))
        
        # Check if we have sufficient data for all timepoints
        if (length(unique(df_test[[timepoint_col]])) < 2) {
            cat("Not enough timepoints with data for:", var_name, "\n")
            next
        }
        
        # Skip if not enough data
        if (nrow(df_test) < 3) {
            cat("Not enough data to run Wilcoxon test for:", var_name, "\n")
            next
        }
        
        # Initialize timepoint p-values as NA
        tp1_p_value <- tp2_p_value <- tp3_p_value <- tp4_p_value <- NA
        
        # Perform Wilcoxon test comparing each timepoint to TP0
        try({
            if ("TP1" %in% df_test[[timepoint_col]]) {
                tp1_p_value <- wilcox.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP1"], 
                                           df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], na.rm = TRUE)$p.value
            }
            if ("TP2" %in% df_test[[timepoint_col]]) {
                tp2_p_value <- wilcox.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP2"], 
                                           df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], na.rm = TRUE)$p.value
            }
            if ("TP3" %in% df_test[[timepoint_col]]) {
                tp3_p_value <- wilcox.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP3"], 
                                           df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], na.rm = TRUE)$p.value
            }
            if ("TP4" %in% df_test[[timepoint_col]]) {
                tp4_p_value <- wilcox.test(df_test[[var_name]][df_test[[timepoint_col]] == "TP4"], 
                                           df_test[[var_name]][df_test[[timepoint_col]] == "TP0"], na.rm = TRUE)$p.value
            }
        }, silent = TRUE)
        
        # Save the Wilcoxon p-values to the results dataframe
        results <- rbind(results, data.frame(Variable = var_name, 
                                             P_Value = min(tp1_p_value, tp2_p_value, tp3_p_value, tp4_p_value, na.rm = TRUE),
                                             TP1_P_Value = tp1_p_value, TP2_P_Value = tp2_p_value,
                                             TP3_P_Value = tp3_p_value, TP4_P_Value = tp4_p_value))
        
        # Generate and save the plots if required
        if (plot_save == "y") {
            median_TP0 <- median(df_test[df_test[[timepoint_col]] == "TP0", var_name], na.rm = TRUE)
            my_colors <- c("darkgreen", "orange", "red", "magenta", "purple")
            
            # Boxplot generation
            box_plot <- ggboxplot(df_test, x = timepoint_col, y = var_name, color = timepoint_col, 
                                  add = "jitter", legend = "none", 
                                  ylab = paste(var_name, "percentage"), width = 0.8, 
                                  add.params = list(size = 1, alpha = 0.5)) +
                geom_hline(yintercept = median_TP0, linetype = 2) +
                stat_compare_means(label = "p.signif", 
                                   method = "wilcox", 
                                   ref.group = "TP0", 
                                   hide.ns = TRUE, 
                                   label.y = max(df_test[[var_name]])) +
                scale_color_manual(values = my_colors)
            
            # Save boxplot
            boxplot_file_name <- file.path(plots_folder, paste(var_name, "_Boxplot_Wilcox.png", sep = ""))
            ggsave(filename = boxplot_file_name, plot = box_plot)
            
            # Barplot generation with SEM
            mean_values <- df_test %>%
                group_by(get(timepoint_col)) %>%
                summarise(mean_value = mean(get(var_name), na.rm = TRUE),
                          sd_value = sd(get(var_name), na.rm = TRUE),
                          n = n()) %>%
                mutate(SEM = sd_value / sqrt(n)) %>%
                rename(Timepoint = `get(timepoint_col)`)
            
            #my_grey_scale <- c("grey", "black", "black", "black", "black") 
            
            # Barplot with SEM and significance
            bar_plot <- ggbarplot(mean_values, x = "Timepoint", y = "mean_value", fill = "Timepoint",
                                  ylab = paste(var_name, "mean value"), 
                                  add = "mean_se", width = 0.8) +
                scale_fill_manual(values = my_colors) +
                geom_errorbar(aes(ymin = mean_value - SEM, ymax = mean_value + SEM), width = 0.2) +
                stat_compare_means(data = df_test, aes(x = get(timepoint_col), y = get(var_name)),
                                   method = "wilcox", ref.group = "TP0", hide.ns = TRUE,
                                   label = "p.signif", label.y = max(mean_values$mean_value) + 0.1 * max(mean_values$mean_value))
            
            # Save barplot
            barplot_file_name <- file.path(plots_folder, paste(var_name, "_Barplot_Wilcox_SEM.png", sep = ""))
            ggsave(filename = barplot_file_name, plot = bar_plot)
            
        }
    }
    
    # Write the results to a CSV file
    results_file <- file.path(paste(results_name, "_Results.csv", sep = ""))
    write.csv(results, file = results_file, row.names = FALSE)
    
    # Return the results dataframe
    return(results)
}

## *Data* ---------------------------------------------------------------------

#### Loading of the Data -----------------------------------------
setwd("/Users/ju5263ta/Github/Monocytes/rawData_Stroke")
getwd ()
# make sure that before you check the files, that not an empty column has been added to the end of the file, remove and export as .csv otherwise.
files <- system( "ls *B2M.csv", intern=T)

# if you just have one data set:
i=1
raw_data <- read.csv(text=readLines(files[i])[-(1:11)], header = T, sep=',', dec = '.',
                     row.names = 1
                     )

# delete Call.1 column, and resent column names
if ( length(colnames(raw_data)) > 13) {
    raw_data <- raw_data[, !colnames(raw_data) %in% "Call.1"]
    colnames(raw_data) <- c("Name","Type","rConc","Name.1","Type.1","Reference","Value",
                            "Quality","Call","Threshold","Value.1","Quality.1","Call.1")
}
raw_data$Dataset<-files[i]

# for more the 1 data set:
i = 2
for (i in (2:length(files))) {
    y = raw_data
    x <- read.csv(text=readLines(files[i])[-(1:11)], header = T, sep=',', dec = '.',row.names = 1)
    
    if (length(colnames(x))>13) {
        x <- x[, !colnames(x) %in% "Call.1"]
        colnames(x) <- c("Name","Type","rConc","Name.1","Type.1","Reference","Value",
                         "Quality","Call","Threshold","Value.1","Quality.1","Call.1")
    }
    x$Dataset<-files[i]
    y <- rbind(y, x)
    i = i +1
    raw_data = y
} 

#### trimming of data ------------------------------------------------------
data <- cbind.data.frame(Sample=raw_data$Name,
                         Gene=raw_data$Name.1,
                         Value=raw_data$Value,
                         dCT_Value=raw_data$Value.1,
                         Call=raw_data$Call.1)

#### evaluate PASS FAIL------------------------------------------------
data$Value<-ifelse(data$Call=='Flag', yes=35, no=data$Value) #999 
data$Value<-ifelse(data$Call=='mFlag', yes=NA, no=data$Value) #999 
data$dCT_Value<-ifelse(data$Call=='Flag', yes=35, no=data$dCT_Value) #999
data$dCT_Value<-ifelse(data$Call=='mFlag', yes=NA, no=data$dCT_Value) #999 instead of NA

# sort through quality and give NA to everything that Failed (=999) or is above CT 35
data$Value<-ifelse(data$Value>=35, yes=35, no=data$Value) 
data$dCT_Value<-ifelse(data$dCT_Value>=35, yes=35, no=data$dCT_Value) 
#data$dCT_Value<-ifelse(data$dCT_Value>=35, yes=35, no=data$dCT_Value) 

# Trim your data
dataTRIM <-cbind.data.frame(Sample=data$Sample,
                            Gene1=data$Gene,
                            CT_Value=data$Value,
                            dCT_Value=data$dCT_Value,
                            Comment = NA)

# Information on the data ---------------------------------------------------------------------
# adjust the indicators for the genes list here (for the heatmap) - line 138:

dataTRIM$TechnicalReplicate <- raw_data$Name.1
info_genes = t(data.frame(lapply( dataTRIM$TechnicalReplicate, function(x){ 
    ret = unlist( stringr::str_split( x, "[_\\s\\.]")); 
    if ( length(ret) < 2){
        ret = c(ret, rep(0, 2-length(ret) ))
    } 
    ret  } ) ))
colnames(info_genes) = c("Gene", "TechnicalReplicate")

info_genes <- as.data.frame(info_genes)

dataTRIM$TechnicalReplicate <- NULL
dataTRIM = cbind(dataTRIM, info_genes)

# List of Genes
gene_names <- unique(dataTRIM$Gene)

# List of Genes
genes <- unique(dataTRIM$Gene)
genes

# Number of Genes
no_genes <- length(genes)
no_genes

# List of Samples
samples <- unique(dataTRIM$Sample)

# Split technical and (pseudo-)biological replicates....:

#without the B1/B2 segregation, we would have each read of the gene normalized against the mean of the two biological repeats of B2M. This is in theory correct but we lose the information on the mean/SD while doing so through the machine.
#without the T1/T2 separation, we have the plate as a whole normalized against all 4 reads... with it we can split the whole plate into two halves, each standing for its own technical read.

info_rows <- t(data.frame(
  sapply(dataTRIM$Sample, function(x) {
    ret <- unlist(stringr::str_split(x, "[_\\s\\.]"))
    if (length(ret) < 5) {
      ret <- c(ret, rep(NA, 5 - length(ret)))
    }
    ret[1:5]  # Select only the first five elements
  })
))
colnames(info_rows) = c("Subpopulation", "Cell_Count", "SampleID","Timepoint","BiologicalReplicate")

# NRT: no Polymerase, technical control 
# NTC: no Template, neg-Sample control   
# 10xLC:Linarity control
info_rows <- as.data.frame(info_rows)

# Label Subpopulations:
# remove the C in NC to avoid problems:
info_rows$Subpopulation <- gsub("NC","N",as.character(info_rows$Subpopulation))
info_rows$Subpopulation[grep("N",info_rows$Subpopulation)]  <- "nonclassical"
info_rows$Subpopulation[grep("M",info_rows$Subpopulation)]  <- "all"
info_rows$Subpopulation[grep("C",info_rows$Subpopulation)]  <- "classical"
info_rows$Subpopulation[grep("I",info_rows$Subpopulation)]  <- "intermediate"

dataTRIM = cbind(dataTRIM, info_rows)

## remove unnessessary info (for now)
#It's easier to make the imputation on less columnes/variables

# Remove Extra sorts
dataTRIM <- dataTRIM %>% filter(BiologicalReplicate != "Extra")
dataTRIM <- dataTRIM %>% filter(Subpopulation != "noRT")
dataTRIM <- dataTRIM %>% filter(Cell_Count != "200")
# remove all Timepoint 5
dataTRIM <- dataTRIM %>% filter(Timepoint != "TP5")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test01")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test02")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test03")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test04")

dataTRIM$Cell_Count <- NULL
#### Evaluate Housekeeping genes -----------
ShitHK <-table ( c( 
    getMeTheShitSamples(dataTRIM, "ACTB", 30),
    getMeTheShitSamples(dataTRIM, "B2M", 30),
    getMeTheShitSamples(dataTRIM, "GAPDH", 30)
)
)
hist(ShitHK)

failedSAMPLES <- names(ShitHK) [ which(ShitHK > 2) ]

# Set CT & dCT NA for the Sampels with failed HK gene
# Set data$Value to NA for the samples that are in the failedSAMPLES vector
dataTRIM$CT_Value <- ifelse(dataTRIM$Sample %in% failedSAMPLES, yes = NA, no = dataTRIM$CT_Value)
dataTRIM$dCT_Value <- ifelse(dataTRIM$Sample %in% failedSAMPLES, yes = NA, no = dataTRIM$dCT_Value)
dataTRIM$Comment <- ifelse(dataTRIM$Sample %in% failedSAMPLES, yes = "More than 2 houskeeing gene reads failed", no = dataTRIM$Comment)

#### Expand data --------------------------------------------------------
unique_SampleID <- unique(dataTRIM$SampleID[!grepl("^Ctr", dataTRIM$SampleID)])
unique_Gene <- unique(dataTRIM$Gene)
unique_Subpopulation <- unique(dataTRIM$Subpopulation)
unique_Timepoint <- unique(dataTRIM$Timepoint[!grepl("TP0", dataTRIM$Timepoint)])
unique_TechnicalReplicate <- unique(dataTRIM$TechnicalReplicate)
unique_BiologicalReplicate <- unique(dataTRIM$BiologicalReplicate)

expanded_df <- expand.grid(SampleID = unique_SampleID,
                           Gene = unique_Gene,
                           Subpopulation = unique_Subpopulation,
                           Timepoint = unique_Timepoint,
                           TechnicalReplicate = unique_TechnicalReplicate,
                           BiologicalReplicate = unique_BiologicalReplicate)

dataEXPAND <- merge(expanded_df, dataTRIM, by = c("SampleID", "Gene", "Subpopulation", "Timepoint", "TechnicalReplicate", "BiologicalReplicate"), all.x = TRUE)

dataEXPAND$SampleID <- as.character(dataEXPAND$SampleID)
dataEXPAND$Gene <- as.character(dataEXPAND$Gene)
dataEXPAND$Subpopulation <- as.character(dataEXPAND$Subpopulation)
dataEXPAND$Timepoint <- as.character(dataEXPAND$Timepoint)
dataEXPAND$TechnicalReplicate <- as.character(dataEXPAND$TechnicalReplicate)
dataEXPAND$BiologicalReplicate <- as.character(dataEXPAND$BiologicalReplicate)

dataCTR <- dataTRIM[grep("^Ctr", dataTRIM$SampleID), ]

dataEXPAND <- rbind(dataEXPAND, dataCTR)

dataEXPAND$Comment <- ifelse(is.na(dataEXPAND$Sample), yes = "No Fludigm run, imputed Values",  no = dataEXPAND$Comment)

#### Metadata & FACSdata---------------------------------------------------------------
metadataP <- read_excel("Metadata_PatientID.xlsx")
Metadata_GeneID <- read_excel("Metadata_GeneID.xlsx")
Plate_info <- read_excel("241014_PlateID.xlsx")

#### match the right controls
Matched_TP0_Gene <- c("Ctr01", "Ctr02", "Ctr08", "Ctr09", "Ctr11", "Ctr14", "Ctr16", "Ctr17", "Ctr18", "Ctr25", "Ctr27", "Ctr29", "Ctr31","Ctr41", "Ctr43")

Unmatched_TP0_FACS <- c("Ctr16", "Ctr21", "Ctr23", "Ctr30", "Ctr39", "Ctr44", 
                        "Ctr15", "Ctr33")


FACSdata <- as.data.frame(read_excel("241008_StrokeControlMFI_combined.xls", # 241008_StrokeControlMFI_combined.xls for updated one
                                     col_types = c("text", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric")))

# split Sample ID into all info
FACS_info_rows <- t(data.frame(sapply(FACSdata$...1, function(x) {
    ret <- unlist(stringr::str_split(x, "[_\\s\\.]"))
    if (length(ret) < 4) {
        ret <- c(ret, rep(NA, 4 - length(ret)))
    }
    ret[1:4]  # Select only the first five elements
})
))

colnames(FACS_info_rows) = c("Type", "Plate", "Sorting Map","Sorted Sample")
FACS_info_rows <- as.data.frame(FACS_info_rows)
FACS_info_rows$PlateID <- paste(FACS_info_rows$Plate,FACS_info_rows$`Sorted Sample` )
FACSdata <- cbind(FACSdata, FACS_info_rows)

# Match PlateID with the right sample with the right data
FACSdata <- merge(FACSdata, Plate_info, by = "PlateID")
# remove the 5th timepoint
FACSdata <- FACSdata %>% filter(Timepoint != "TP5")

# add the metadata of Patients
FACSdata <- merge(FACSdata, metadataP, by = "SampleID", all.x = TRUE)

# Combine the data frames based on common values in the "id" column
# --> here i loose all the test because they are not in the metadata
dataAll <- merge(dataEXPAND, metadataP, by = "SampleID", all.x = TRUE)

#relabel Male & Female
dataAll$Sex[grep("M",dataAll$Sex)]  <- "Male"
dataAll$Sex[grep("F",dataAll$Sex)]  <- "Female"

FACSdata$Sex[grep("M",FACSdata$Sex)]  <- "Male"
FACSdata$Sex[grep("F",FACSdata$Sex)]  <- "Female"

dataAll <- dataAll %>% filter(!is.na(Age), !is.na(SampleID), !is.na(Category), !is.na(Sex))

#### failed Genes --------------------------------------------------------
fails <- dataAll %>%
    group_by(Subpopulation, Gene, Sex) %>%
    summarize(
        N = n(),
        Count_Above_35_CT = sum(CT_Value >= 35, na.rm = TRUE) + sum(is.na(CT_Value)),  # Include NAs in the count
    )

fails$Per_Above_35_CT = fails$Count_Above_35_CT/fails$N * 100

failed_genes <- fails %>% filter(Per_Above_35_CT >= 82.7) %>% select(Gene)
wrong_threshold <- c("CLEC7A", "S100A8", "S100A9")
exclude_genes  <- unique(c("Xeno",failed_genes$Gene, wrong_threshold))
paste("The gene",unique(failed_genes$Gene), "was excluded, as more than 80% of the reads failed.")
paste("The gene",wrong_threshold, "was excluded, as the global threshold setting did not allow proper evaluation.")

# remove the failed genes:
dataWORKING <- dataAll %>%  filter(!Gene %in% exclude_genes)

#### failed Samples --------------------------------------------------------
fails_Sample <- dataWORKING %>%
    group_by(Sample) %>%
    summarize(
        N = n(),
        Count_Above_35_CT <- sum(CT_Value >= 35, na.rm = TRUE),
        Per_Above_35_CT = mean(CT_Value >= 35, na.rm = TRUE) * 100,
        Count_Above_35_dCT = sum(dCT_Value >= 35, na.rm = TRUE),
        Per_Above_35_dCT = mean(dCT_Value >= 35, na.rm = TRUE) * 100
    )

failed_samples <- fails_Sample %>% filter(Per_Above_35_CT >= 30) %>% select(Sample)
failed_samples <- failed_samples$Sample
paste("The Sample ",failed_samples, "was excluded, as more than 30% of the reads failed.")

# Set all Values from failed Samples (> 30%) NA:
dataWORKING$CT_Value <- ifelse(dataWORKING$Sample %in% failed_samples, yes = NA, no = dataWORKING$CT_Value)
dataWORKING$dCT_Value <- ifelse(dataWORKING$Sample %in% failed_samples, yes = NA, no = dataWORKING$dCT_Value)
dataWORKING$Comment <- ifelse(dataWORKING$Sample %in% failed_samples, yes = "More 30% of CT values were 35 and above, imputed Values", no = dataWORKING$Comment)
unique(dataWORKING$Comment)

# Filter the data for the desired genes
filtered_data <- dataWORKING %>%
    filter(Gene %in% c("B2M", "ACTB", "GAPDH"))

# Calculate the CV for each group
variance_table <- filtered_data %>%
    group_by(Gene, Timepoint, Subpopulation) %>%
    summarize(CV = sd(CT_Value, na.rm = TRUE) / mean(CT_Value, na.rm = TRUE)) %>%
    ungroup()

# Display the variance table
print(variance_table)

# Plotting the variances
ggplot(variance_table, aes(x = Timepoint, y = CV, fill = Subpopulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Gene, scales = "fixed") +
    labs(title = "Coefficient of Variation Across Timepoints and Subpopulations",
         x = "Timepoint", y = "Coefficient of Variation (CV)") +
    theme_minimal() +
    scale_fill_manual(values = c("darkgreen", "orange", "red", "magenta"))


# *Imputation* --------------------------------------------------------------------------------------------
#### Prepare Matrix for imputation ---------------------------------------------------------------------
# visualize the missing values:
vis_miss(dataWORKING)
#hist(dataWORKING$dCT_Value, main = "Histogram of Original Data")

initial2<-mice(dataWORKING, maxit=0, print=F)
# Check
initial2$method
# exclude form the imputation
initial2$method["Sample"]<-"none"
initial2$method["Gene1"]<-"none"
initial2$method["TechnicalReplicate"]<-"none"
initial2$method["BiologicalReplicate"]<-"none"
initial2$method["Comment"]<-"none"

#Check again
initial2$method

# convert to factors
dataWORKING$Sex <- as.factor(dataWORKING$Sex)
dataWORKING$Subpopulation <- as.factor(dataWORKING$Subpopulation)
dataWORKING$Gene <- as.factor(dataWORKING$Gene)
dataWORKING$Category <- as.factor(dataWORKING$Category) # check the effect of ordered
dataWORKING$Timepoint <- factor(dataWORKING$Timepoint, ordered = TRUE)

# Matrix of predictors:
initial2$predictorMatrix

#### Impute Gene expression --------------------------------------------------------------------------------------------
data_imputed<-mice(dataWORKING, m=20, maxit=10, seed=1234, meth=initial2$method, pred=initial2$predictorMatrix)

# extract the imputed data into a “normal” data frame, run the analyses on each imputation separately and do model diagnostics.
finaldata <- complete(data_imputed, "long")
#View(finaldata)
table(finaldata$.imp)

data_imputed_mean <- finaldata  %>%
    group_by(SampleID, Gene, Subpopulation,Timepoint, TechnicalReplicate, BiologicalReplicate, Sample, Gene1, Age, Sex, Category) %>%    
    summarize(CT = mean(CT_Value, na.rm = TRUE),
             # CT_sd = sd(CT_Value, na.rm = TRUE), 
             # CT_N = sum(!is.na(CT_Value)), 
              dCT_B2M = mean(dCT_Value, na.rm = TRUE), 
             # dCT_B2M_sd = sd(dCT_Value, na.rm = TRUE),
              .groups = "drop")

dataFINAL <- as.data.frame(data_imputed_mean)

# Variance after imputation
# Filter the data for the desired genes
filtered_data <- dataWORKING %>%
    filter(Gene %in% c("B2M", "ACTB", "GAPDH"))

# Calculate the CV for each group
variance_table <- filtered_data %>%
    group_by(Gene, Timepoint, Subpopulation) %>%
    summarize(CV = sd(CT_Value, na.rm = TRUE) / mean(CT_Value, na.rm = TRUE)) %>%
    ungroup()

# Display the variance table
print(variance_table)

# Plotting the variances
ggplot(variance_table, aes(x = Timepoint, y = CV, fill = Subpopulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Gene, scales = "fixed") +
    labs(title = "Coefficient of Variation Across Timepoints and Subpopulations",
         x = "Timepoint", y = "Coefficient of Variation (CV)") +
    theme_minimal() +
    scale_fill_manual(values = c("darkgreen", "orange", "red", "magenta"))

# *Normalize the data vs the  Geometric mean of each Sample* ----------------

# Calculate the CT median or mean
dataFINAL$CTmedian <- NA
dataFINAL$CTmean <- NA
dataFINAL$CTgmean <- NA  # For geometric mean

dataFINAL <- dataFINAL %>%
    group_by(Sample) %>%
    mutate(
       CTmedian = median(CT[CT < 35], na.rm = TRUE),
       CTmean = mean(CT[CT < 35], na.rm = TRUE),
       CTgmean = exp(mean(log(CT[CT < 35]), na.rm = TRUE))  # Geometric mean
    )
# normalize all CT value to the Media CT value of each Sample
dataFINAL$dCT_median <- dataFINAL$CT - dataFINAL$CTmedian
dataFINAL$dCT_mean <- dataFINAL$CT - dataFINAL$CTmean
dataFINAL$dCT_gmean <- dataFINAL$CT - dataFINAL$CTgmean

# check the variance again after the imputation
filtered_data <- dataFINAL %>%
    filter(Gene %in% c("B2M", "ACTB", "GAPDH"))

# Calculate the CV for each group
variance_table <- filtered_data %>%
    group_by(Gene, Timepoint, Subpopulation) %>%
    summarize(CV = sd(dCT_median, na.rm = TRUE) / mean(dCT_median, na.rm = TRUE)) %>%
    ungroup()

# Display the variance table
print(variance_table)

# Plotting the variances
ggplot(variance_table, aes(x = Timepoint, y = CV, fill = Subpopulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Gene, scales = "fixed") +
    labs(title = "Coefficient of Variation Across Timepoints and Subpopulations",
         x = "Timepoint", y = "Coefficient of Variation (CV)") +
    theme_minimal() +
    scale_fill_manual(values = c("darkgreen", "orange", "red", "magenta"))

# *Zscore: Normalize the data vs Gene expression ----------------
# as such: apply( data_summary , 1 , function(x) { (x["mean_Z"]- x["MeanForGene"])/ x["sd_for_gene"] } )
dataFINAL <- dataFINAL %>%
    group_by(Gene) %>%
    mutate(
        mean_dCT_Gene = mean(dCT_median, na.rm = TRUE),  # Mean of mean_rE by Gene
        sd_dCT_Gene = sd(dCT_median, na.rm = TRUE),      # Standard deviation of mean_rE by Gene
        z_score_dCT = (dCT_median - mean_dCT_Gene) / sd_dCT_Gene,  # Z-score calculation
        inverted_z_score = -z_score_dCT,
        capped_z_score = pmin(pmax(inverted_z_score, -3), 3) 
        ) %>%
    ungroup() 

#### Calculate rE ----------------

# change here to median if you want
dataFINAL$rE_Value <- 2^((-1)*dataFINAL$dCT_median)
#dataFINAL$rE_Value <- 2^((-1)*dataFINAL$dCT_median)

dataFINAL$rE_Value <- ifelse(dataFINAL$CT>=35, yes=0, no=dataFINAL$rE_Value)
dataFINAL$rE_Value <- ifelse(dataFINAL$dCT_median >=35, yes=0, no=dataFINAL$rE_Value) 

# Check but if the g-median is above 30 the run was probably not successfull
dataFINAL$rE_Value <- ifelse(dataFINAL$CTgmean>=30, yes=0, no=dataFINAL$rE_Value)
    
#### Gene Info  ----------------
## Gene Metadata & other edits:
Gene_Info <- data.frame(Gene = Metadata_GeneID$Gene, GeneID = Metadata_GeneID$GeneID, Gene_Group = Metadata_GeneID$Gene_Group)

dataFINAL <- merge(dataFINAL, Gene_Info, by = "Gene")

# List of Patients (2x 15))
patients <- unique(dataFINAL$SampleID)
no_patients <- length(patients)

#### Days Post-Stroke data ----------------
# DaysPS <- read_excel("Metadata_DaysPS.xls")
dataFINAL$DaysPS <- as.numeric(NA)
dataFINAL$DaysPS[grep("TP0",dataFINAL$Timepoint)]  <- 0
dataFINAL$DaysPS[grep("TP1",dataFINAL$Timepoint)]  <- 1
dataFINAL$DaysPS[grep("TP2",dataFINAL$Timepoint)]  <- 4
dataFINAL$DaysPS[grep("TP3",dataFINAL$Timepoint)]  <- 30
dataFINAL$DaysPS[grep("TP4",dataFINAL$Timepoint)]  <- 90
#dataFINAL$DaysPS[grep("TP5",dataFINAL$Timepoint)]  <- 360


## MEAN values for the ANOVA & so ----------------
dataFINALmean <- dataFINAL %>%
        group_by(SampleID, Sex, Age, Category, Timepoint, DaysPS,
                 Subpopulation, Gene, GeneID, Gene_Group) %>%
        summarise(mean_CT = mean(CT), sd_CT = sd(CT),
                  mean_dCT = mean(dCT_median), sd_dCT = sd(dCT_median),
                  mean_Z = mean(inverted_z_score), sd_Z = sd(inverted_z_score),
                  mean_capZ = mean(capped_z_score), sd_capZ = sd(capped_z_score),
                  mean_rE = mean(rE_Value), sd_rE = sd(rE_Value)
                  )%>%
  ungroup()

#### File Location ----------------
today <- Sys.Date()
output_location <- paste(today,"_Stroke_Results", sep="")
setwd(paste("/Users/ju5263ta/Github/Monocytes/Data/",output_location, "/",sep=""))
getwd()

write_xlsx(fails, "Fail_rates.xlsx")
write_xlsx(dataFINAL, "dataFINAL.xlsx")
write_xlsx(dataFINALmean, "dataFINALmean.xlsx")

# *Dotplot of expressed Genes* -----------------------------------------------------------------------------
createFolder("Dotplots_Gene_summary")

# Create a new column for Timepoint and Subpopulation combination
dataFINALmean <- dataFINALmean %>%
    mutate(Sample_Combo = paste(Timepoint, Subpopulation, sep = "_"))

dataFINALmean <- dataFINALmean %>%
    mutate(Sample_Combo2 = paste(Subpopulation, Timepoint, sep = "_"))

# Create a combined column for Gene and GeneID for better sorting
dataFINALmean <- dataFINALmean %>%
    mutate(Gene_Combined = paste(GeneID, Gene, sep = "_"))  # Combine GeneID and Gene

# Using the function for data_summary and data_summary2
data_summary <- calculate_summary(dataFINALmean, "Sample_Combo", "Gene_Combined")
data_summary2 <- calculate_summary(dataFINALmean, "Sample_Combo2", "Gene_Combined")

data_summary <- data_summary %>%
    mutate(Gene_Combined = factor(Gene_Combined, levels = rev(unique(Gene_Combined))))

data_summary2 <- data_summary2 %>%
    mutate(Gene_Combined = factor(Gene_Combined, levels = rev(unique(Gene_Combined))))

# Create a linear gradient function for the original plots
linear_gradient <- function() {
    # Define the colors for the gradient
    colors <- c("blue", "purple", "red")
    
    # Set the breakpoints for the gradient (-1.2 to 1.2)
    scale_color_gradientn(colors = colors,
                          limits = c(-1.2, 1.2),  # Set the limits from -1.2 to 1.2
                          guide = "colorbar",
                          na.value = "grey50")  # Color for NA values
}

ggplot(data_summary, aes(x = Sample_Combo, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_Z)) +  # Use Dot_Size for scaling
    scale_size(range = c(1, 5)) +  # Size range reflects 1-5 scaling
    linear_gradient() +  # Use the custom linear gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of qPCR Data (high Expr. in red; low Expr. in blue)",
        color = "Gene Expr. \n(Zscored)",
        size = "high SD (small)\nlow SD (big)"  # Use "\n" for line break
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "white"),  # Set the plot background to white
        panel.background = element_rect(fill = "white"),  # Set the panel background to white
        panel.grid.major = element_line(color = "grey", size = 0.5),  # Optional: adjust grid lines
        panel.grid.minor = element_line(color = "lightgrey", size = 0.25)  # Optional: adjust minor grid lines
    )
ggsave(filename = "Dotplots_Gene_summary/Exp_by_TP.png")

ggplot(data_summary2, aes(x = Sample_Combo2, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_Z)) +  # Use Dot_Size for scaling
    scale_size(range = c(1, 5)) +  # Size range reflects 1-5 scaling
    linear_gradient() +  # Use the custom linear gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of qPCR Data (high Expr. in red; low Expr. in blue)",
        color = "Gene Expr. \n(Zscored)",
        size = "high SD (small)\nlow SD (big)"  # Use "\n" for line break
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "white"),  # Set the plot background to white
        panel.background = element_rect(fill = "white"),  # Set the panel background to white
        panel.grid.major = element_line(color = "grey", size = 0.5),  # Optional: adjust grid lines
        panel.grid.minor = element_line(color = "lightgrey", size = 0.25)  # Optional: adjust minor grid lines
    )
ggsave(filename = "Dotplots_Gene_summary/Exp_by_Suptype.png")

# Create a linear gradient function for the original plots
linear_gradient <- function() {
    # Define the colors for the gradient
    colors <- c("white", "grey", "black")
    
    # Set the breakpoints for the gradient (-1.2 to 1.2)
    scale_color_gradientn(colors = colors,
                          limits = c(-1.2, 1.2),  # Set the limits from -1.2 to 1.2
                          guide = "colorbar",
                          na.value = "red")  # Color for NA values
}

ggplot(data_summary, aes(x = Sample_Combo, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_Z)) +  # Use Dot_Size for scaling
    scale_size(range = c(1, 5)) +  # Size range reflects 1-5 scaling
    linear_gradient() +  # Use the custom linear gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of qPCR Data (high Expr. in red; low Expr. in blue)",
        color = "Gene Expr. \n(Zscored)",
        size = "high SD (small)\nlow SD (big)"  # Use "\n" for line break
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "white"),  # Set the plot background to white
        panel.background = element_rect(fill = "white"),  # Set the panel background to white
        panel.grid.major = element_line(color = "grey", size = 0.5),  # Optional: adjust grid lines
        panel.grid.minor = element_line(color = "lightgrey", size = 0.25)  # Optional: adjust minor grid lines
    )
ggsave(filename = "Dotplots_Gene_summary/Exp_by_TP_grey.png")

ggplot(data_summary2, aes(x = Sample_Combo2, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_Z)) +  # Use Dot_Size for scaling
    scale_size(range = c(1, 5)) +  # Size range reflects 1-5 scaling
    linear_gradient() +  # Use the custom linear gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of qPCR Data (high Expr. in red; low Expr. in blue)",
        color = "Gene Expr. \n(Zscored)",
        size = "high SD (small)\nlow SD (big)"  # Use "\n" for line break
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "white"),  # Set the plot background to white
        panel.background = element_rect(fill = "white"),  # Set the panel background to white
        panel.grid.major = element_line(color = "grey", size = 0.5),  # Optional: adjust grid lines
        panel.grid.minor = element_line(color = "lightgrey", size = 0.25)  # Optional: adjust minor grid lines
    )
ggsave(filename = "Dotplots_Gene_summary/Exp_by_Suptype_grey.png")



## LinReg FACS outliers-----------------------------------------
createFolder("LinReg_FACS_Plots")

# Define custom colors for the Timepoints
my_colors <- c("TP0" = "darkgreen", "TP1" = "orange", "TP2" = "red", "TP3" = "magenta", "TP4" = "purple")

# Initialize data frames to store results
Age_correlation_FACS <- data.frame(Cells = character(), Timepoint = character(), Sex = character(), N = numeric(), 
                                   p.value = numeric(), Coefficient = numeric(), Down95 = numeric(), 
                                   Up95 = numeric(), Intercept = numeric(), Fold = numeric(), 
                                   Max_age = numeric(), Min_age = numeric())

outliers_FACS <- FACSdata[0,]
outliers_FACS$Outlier <- rep("NA", nrow(outliers_FACS))

# Create folder to save plots if it doesn't exist
output_folder <- "LinReg_FACS_Plots"
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}

# remove Patient40 as I do not have age and Sex from them
FACSdata <- FACSdata[!is.na(FACSdata$Age), ]

# Loop over the relevant columns (assumed to be columns 4 to 21)
for (i in 4:16) {
    cells_colname <- colnames(FACSdata)[i]
    
    #### All ---------
    rows_all <- calculate_linreg_by_timepoint(FACSdata, cells_colname)
    rows_all$Sex <- "All"  # Add the Sex label for all samples
    Age_correlation_FACS <- rbind(Age_correlation_FACS, rows_all)
    
    # Plot regression for "All" samples with custom colors, size, and theme
    p_all <- ggplot(FACSdata, aes(x = Age, y = .data[[cells_colname]], color = Timepoint)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", aes(color = Timepoint), se = FALSE, size = 1) +  # Regression line in corresponding color
        ggtitle(paste("Linear Regression - All Samples (", cells_colname, ")", sep = "")) +
        xlab("Age") +
        ylab(cells_colname) +
        theme(text = element_text(size = 14)) +
        theme(
            panel.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "white", color = "black"),
            text = element_text(color = "black"),
            panel.grid.major = element_line(color = "gray", size = 0.2),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()
        ) +
        scale_color_manual(values = my_colors) +  # Apply custom color scale for points and lines
        scale_fill_manual(values = my_colors)  # Apply custom fill scale for confidence intervals
    
    # Add confidence interval only if p-value is <= 0.05
    for (timepoint in unique(FACSdata$Timepoint)) {
        if (rows_all$p.value[rows_all$Timepoint == timepoint] <= 0.05) {
            p_all <- p_all + geom_smooth(data = subset(FACSdata, Timepoint == timepoint), 
                                         method = "lm", 
                                         aes(fill = Timepoint), 
                                         se = TRUE, 
                                         alpha = 0.25, 
                                         color = NA)  # CI without color for transparency
        }
    }
    
    # Save the plot for "All" samples
    ggsave(filename = file.path(output_folder, paste0("LinReg_All_", cells_colname, ".png")), plot = p_all)
    
    ### Male ----
    FACSdata_male <- FACSdata[FACSdata$Sex == "Male", ]
    if (nrow(FACSdata_male) > 0) {
        rows_male <- calculate_linreg_by_timepoint(FACSdata_male, cells_colname)
        rows_male$Sex <- "Male"  # Add the Sex label for male samples
        Age_correlation_FACS <- rbind(Age_correlation_FACS, rows_male)
        
        # Plot regression for "Male" samples with custom colors, size, and theme
        p_male <- ggplot(FACSdata_male, aes(x = Age, y = .data[[cells_colname]], color = Timepoint)) +
            geom_point(size = 2) +
            geom_smooth(method = "lm", aes(color = Timepoint), se = FALSE, size = 1) +  # Regression line in corresponding color
            ggtitle(paste("Linear Regression - Male (", cells_colname, ")", sep = "")) +
            xlab("Age") +
            ylab(cells_colname) +
            theme(text = element_text(size = 14)) +
            theme(
                panel.background = element_rect(fill = "white", color = "black"),
                plot.background = element_rect(fill = "white", color = "black"),
                text = element_text(color = "black"),
                panel.grid.major = element_line(color = "gray", size = 0.2),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_manual(values = my_colors) +  # Apply custom color scale for points and lines
            scale_fill_manual(values = my_colors)  # Apply custom fill scale for confidence intervals
        
        # Add confidence interval only if p-value is <= 0.05
        for (timepoint in unique(FACSdata_male$Timepoint)) {
            if (rows_male$p.value[rows_male$Timepoint == timepoint] <= 0.05) {
                p_male <- p_male + geom_smooth(data = subset(FACSdata_male, Timepoint == timepoint), 
                                               method = "lm", 
                                               aes(fill = Timepoint), 
                                               se = TRUE, 
                                               alpha = 0.25, 
                                               color = NA)  # CI without color for transparency
            }
        }
        
        # Save the plot for "Male" samples
        ggsave(filename = file.path(output_folder, paste0("LinReg_Male_", cells_colname, ".png")), plot = p_male)
    }
    
    ### Female ---------
    FACSdata_female <- FACSdata[FACSdata$Sex == "Female", ]
    if (nrow(FACSdata_female) > 0) {
        rows_female <- calculate_linreg_by_timepoint(FACSdata_female, cells_colname)
        rows_female$Sex <- "Female"  # Add the Sex label for female samples
        Age_correlation_FACS <- rbind(Age_correlation_FACS, rows_female)
        
        # Plot regression for "Female" samples with custom colors, size, and theme
        p_female <- ggplot(FACSdata_female, aes(x = Age, y = .data[[cells_colname]], color = Timepoint)) +
            geom_point(size = 2) +
            geom_smooth(method = "lm", aes(color = Timepoint), se = FALSE, size = 1) +  # Regression line in corresponding color
            ggtitle(paste("Linear Regression - Female (", cells_colname, ")", sep = "")) +
            xlab("Age") +
            ylab(cells_colname) +
            theme(text = element_text(size = 14)) +
            theme(
                panel.background = element_rect(fill = "white", color = "black"),
                plot.background = element_rect(fill = "white", color = "black"),
                text = element_text(color = "black"),
                panel.grid.major = element_line(color = "gray", size = 0.2),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank()
            ) +
            scale_color_manual(values = my_colors) +  # Apply custom color scale for points and lines
            scale_fill_manual(values = my_colors)  # Apply custom fill scale for confidence intervals
        
        # Add confidence interval only if p-value is <= 0.05
        for (timepoint in unique(FACSdata_female$Timepoint)) {
            if (rows_female$p.value[rows_female$Timepoint == timepoint] <= 0.05) {
                p_female <- p_female + geom_smooth(data = subset(FACSdata_female, Timepoint == timepoint), 
                                                   method = "lm", 
                                                   aes(fill = Timepoint), 
                                                   se = TRUE, 
                                                   alpha = 0.25, 
                                                   color = NA)  # CI without color for transparency
            }
        }
        
        # Save the plot for "Female" samples
        ggsave(filename = file.path(output_folder, paste0("LinReg_Female_", cells_colname, ".png")), plot = p_female)
    }
}

# Save results to CSV
write.csv(Age_correlation_FACS, "Age_correlation_FACS.csv", row.names = FALSE)

# *1WAY ANOVA* FACS----------------------------------------------------------------------------------------------------
# with wilcox
ANOVA_FACSall <- automate_anova_extraction(output_location, "Wilcox_FACS_Plots", 
                                                  "Wilcox_FACS", "y", FACSdata, colnames(FACSdata)[4:16], 
                                                  "Timepoint")
# paired t.test
FACSdata_matched <- FACSdata %>%
    filter(!SampleID %in% Unmatched_TP0_FACS)

paired_Ttest_FACS <- automate_ttest_extraction(output_location, "paired_Ttest_FACS_Plots",
                          "paired_Ttest_FACS", "y", FACSdata_matched, colnames(FACSdata_matched)[4:16],
                          "Timepoint")


# *Gene Expr. statistics* ----------------------------------------------------------------------------------------------------
# As they are non-normal disributed, dependent with similar variance
# Adjust to yes if you want plots to be saved or n if not
plot_save <- "y"

# remove unmatched CTR for ANOVA too 
data_rE_mean_matched <- dataFINALmean %>%
    filter(!(Timepoint == "TP0" & !SampleID %in% Matched_TP0_Gene))
# Check if all groups have the same size:
check_sample_counts(data_rE_mean_matched)

#### Wilcox test for rE of Genes --------------------------------------------------
folder <- "Wilcox_rE_Plots"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric(),
                      TP1_P_Value = numeric(), TP2_P_Value = numeric(), TP3_P_Value = numeric(), TP4_P_Value = numeric())

j = 1
for (j in 1:length(unique(data_rE_mean_matched$Gene))) {
    gene <- data_rE_mean_matched$Gene[j]
    df_gene <- filter(data_rE_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Category), !is.na(mean_rE), !is.na(Subpopulation), !is.na(Timepoint))
    i = 1
    for (i in 1:length(unique(df_gene$Subpopulation))) {
        sup <- df_gene$Subpopulation[i]
        df <- filter(df_gene, Subpopulation == sup)
        if (nrow(df) <= 0 || gene %in% exclude_genes) {
            cat("Skipping:", gene, "for subpopulation:", sup, "due to insufficient data or in excluded genes list\n")
            next  # Skip to the next iteration of the loop
        }
        
        df <- df[!df$mean_rE %in% boxplot.stats(df$mean_rE)$out, ]
        df$Timepoint <- factor(as.vector(df$Timepoint))
        baseline <- filter(df, Timepoint == "TP0")
        median_TP0 <- median(baseline$mean_rE)
        
        # Perform Wilcoxon test for each comparison with TP0
        timepoints <- c("TP1", "TP2", "TP3", "TP4")
        tp_p_values <- c(NA, NA, NA, NA)  # Initialize p-values for TP1 to TP4
        
        # Loop through each timepoint and compare to TP0
        for (tp in seq_along(timepoints)) {
            comparison_tp <- timepoints[tp]
            
            # Check if there are enough data for both TP0 and the current timepoint
            if (nrow(filter(df, Timepoint == comparison_tp)) > 1 && nrow(baseline) > 1) {
                # Perform Wilcoxon test between TP0 and the current timepoint
                test_result <- wilcox.test(df$mean_rE[df$Timepoint == "TP0"], 
                                           df$mean_rE[df$Timepoint == comparison_tp],
                                           paired = FALSE, exact = FALSE)
                tp_p_values[tp] <- test_result$p.value  # Store the p-value
            }
        }
        
        # Save results to data frame
        results <- rbind(results, data.frame(Gene = gene, Subpopulation = sup, 
                                             P_Value = min(tp_p_values, na.rm = TRUE),  # Store the smallest p-value as general
                                             TP1_P_Value = tp_p_values[1], 
                                             TP2_P_Value = tp_p_values[2], 
                                             TP3_P_Value = tp_p_values[3], 
                                             TP4_P_Value = tp_p_values[4]))
        
        if (plot_save == "y") {
            my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
            
            ggboxplot(df, 
                      x = "Timepoint", 
                      y = "mean_rE", 
                      color = "Timepoint", 
                      add = "jitter", 
                      legend = "none", 
                      ylab = paste(gene, "expression relative to the geometric Mean"), 
                      width = 0.8, 
                      add.params = list(size = 1, alpha = 0.5)) +  
                geom_hline(yintercept = median_TP0, linetype = 2) + 
                #stat_compare_means(method = "wilcox", label.y = max(df$mean_rE)) +        
                stat_compare_means(label = "p.signif", 
                                   method = "wilcox", 
                                   ref.group = "TP0", 
                                   hide.ns = TRUE, 
                                   label.y = max(df$mean_rE)) +
                scale_color_manual(values = my_colors)
            
            file_name <- file.path(folder, paste(gene,"_", sup, "_rE_Wilcoxon.png", sep = ""))
            ggsave(filename = file_name)
        }
    }
}
Wilcox_rE <- results
Wilcox_rE$Significance <- sapply(Wilcox_rE$P_Value, get_significance)
Wilcox_rE$TP1_Significance <- sapply(Wilcox_rE$TP1_P_Value, get_significance)
Wilcox_rE$TP2_Significance <- sapply(Wilcox_rE$TP2_P_Value, get_significance)
Wilcox_rE$TP3_Significance <- sapply(Wilcox_rE$TP3_P_Value, get_significance)
Wilcox_rE$TP4_Significance <- sapply(Wilcox_rE$TP4_P_Value, get_significance)
write.csv(Wilcox_rE, file = "Wilcox_Results_rE.csv", row.names = FALSE)

#### ANOVA for Zscore with Baseline --------------------------------------------------
# ANOVA repeated mesures as we have dependent, normally distributed values with similar varaiance
folder <- "ANOVA_Zscore_Plots"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric(),
                      TP1_P_Value = numeric(), TP2_P_Value = numeric(), TP3_P_Value = numeric(), TP4_P_Value = numeric())
j = 1
for (j in 1:length(unique(data_rE_mean_matched$Gene))) {
    gene <- data_rE_mean_matched$Gene[j]
    df_gene <- filter(data_rE_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Category), !is.na(mean_Z), !is.na(Subpopulation), !is.na(Timepoint))
    i=1
    for (i in 1:length(unique(df_gene$Subpopulation))) {
        sup <- df_gene$Subpopulation[i]
        df <- filter(df_gene, Subpopulation == sup)
        if (nrow(df) <= 0 || gene %in% exclude_genes) {
            cat("Skipping:", gene, "for subpopulation:", sup, "due to insufficient data or in excluded genes list\n")
            next  # Skip to the next iteration of the loop
        }
        
        df <- df[!df$mean_Z %in% boxplot.stats(df$mean_Z)$out, ]
        df$Timepoint <- factor(as.vector(df$Timepoint))
        baseline <- filter(df, Timepoint == "TP0")
        median_TP0 <- median(baseline$mean_Z)
        
        # Perform ANOVA
        anova_result <- aov(mean_Z ~ Timepoint, data = df)
        
        # Perform Tukey's post-hoc test
        tukey_result <- TukeyHSD(anova_result, "Timepoint")
        
        # Extract p-values for each timepoint comparison against TP0
        tp1_p_value <- tukey_result$Timepoint["TP1-TP0", "p adj"]
        tp2_p_value <- tukey_result$Timepoint["TP2-TP0", "p adj"]
        tp3_p_value <- tukey_result$Timepoint["TP3-TP0", "p adj"]
        tp4_p_value <- tukey_result$Timepoint["TP4-TP0", "p adj"]
        
        # Save results to data frame
        results <- rbind(results, data.frame(Gene = gene, Subpopulation = sup, 
                                             P_Value = summary(anova_result)[[1]]$"Pr(>F)"[1],
                                             TP1_P_Value = tp1_p_value, TP2_P_Value = tp2_p_value, 
                                             TP3_P_Value = tp3_p_value, TP4_P_Value = tp4_p_value))
        
        if (plot_save == "y") {
            my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
            
            ggboxplot(df, 
                      x = "Timepoint", 
                      y = "mean_Z", 
                      color = "Timepoint", 
                      add = "jitter", 
                      legend = "none", 
                      ylab = paste(gene, "expression relative to B2M"), 
                      width = 0.8, 
                      add.params = list(size = 1, alpha = 0.5)) +  
                geom_hline(yintercept = median_TP0, linetype = 2) + 
                stat_compare_means(method = "anova", label.y = max(df$mean_Z)) +        
                stat_compare_means(label = "p.signif", 
                                   method = "wilcox", 
                                   ref.group = "TP0", 
                                   hide.ns = TRUE, 
                                   label.y = max(df$mean_Z)) +
                scale_color_manual(values = my_colors)
            
            file_name <- file.path(folder, paste(gene,"_", sup, "_Zscore_ANOVAwilcox.png", sep = ""))
            ggsave(filename = file_name)
            
            # Timeline plot
                title_facet <- paste(gene, "expression in", sup, "monocytes (Zscored dCT)", sep=" ")
                
            ggline(subset(df), 
                       x = "Timepoint", 
                       y = "mean_Z",  
                       add = "mean_se",  # Add mean and standard error
                       size = 1.2) +  
                    labs(y = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                         title = title_facet) +
                    theme_bw() +  # Apply white background
                    theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
                          panel.grid.minor = element_blank(), 
                          panel.border = element_blank(),
                          axis.line = element_line(color = "black")) +
                    stat_compare_means(label = "p.signif", 
                                       method = "anova", #wilcox before
                                       ref.group = "TP0", 
                                       hide.ns = TRUE, 
                                       label.y = max(df$mean_Z))
                
                    file_name_facet <- paste0(folder, "/", title_facet, "_Mean_SEM.png")
                    ggsave(filename = file_name_facet)
        }
    }
}
ANOVA_Zscore <- results
ANOVA_Zscore$Significance <- sapply(ANOVA_Zscore$P_Value, get_significance)
ANOVA_Zscore$TP1_Significance <- sapply(ANOVA_Zscore$TP1_P_Value, get_significance)
ANOVA_Zscore$TP2_Significance <- sapply(ANOVA_Zscore$TP2_P_Value, get_significance)
ANOVA_Zscore$TP3_Significance <- sapply(ANOVA_Zscore$TP3_P_Value, get_significance)
ANOVA_Zscore$TP4_Significance <- sapply(ANOVA_Zscore$TP4_P_Value, get_significance)
write.csv(ANOVA_Zscore, file = "ANOVA_Results_Zscore.csv", row.names = FALSE)

#### *Male vs Female Analysis (Mean ± SEM)* -------------------------------------
folder <- "Mean_SEM_Sex_Plots"
createFolder(folder)  # Ensure the folder exists
results <- data.frame(Sex = character(), Gene = character(), Subpopulation = character(), 
                      P_Value = numeric(), TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                      TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (i in 1:length(unique(data_rE_mean_matched$Gene))) {
    gene <- unique(data_rE_mean_matched$Gene)[i]
    df_gene <- filter(data_rE_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_Z), !is.na(Subpopulation), !is.na(Sex))
    
    for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- unique(df_gene$Subpopulation)[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
        
        # Remove outliers
        df_subpop_clean <- df_subpop[!df_subpop$mean_Z %in% boxplot.stats(df_subpop$mean_Z)$out, ]
        
        if (nrow(df_subpop_clean) <= 0) next  # Skip if no rows remain
        
        tryCatch({
            title_facet <- paste(gene, "expression in", sup, "monocytes (Zscored dCT)", sep=" ")
            
            # Line plot comparing Male vs Female with mean ± SEM across Timepoints
            ggline(subset(df_subpop_clean, !is.na(Sex)), 
                   x = "Timepoint", 
                   y = "mean_Z", 
                   color = "Sex", 
                   add = "mean_se",  # Add mean and standard error
                   size = 1.2) +  
                labs(y = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                     title = title_facet) +
                scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2")) +
                theme_bw() +  # Apply white background
                theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
                      panel.grid.minor = element_blank(), 
                      panel.border = element_blank(),
                      axis.line = element_line(color = "black")) +
                stat_compare_means(method = "wilcox.test",  # or "t.test" 
                                   aes(group = Sex),
                                   label = "p.signif",  
                                   size = 5) 
            
            if (plot_save == "y") {
                file_name_facet <- paste0(folder, "/", title_facet, "_Mean_SEM.png")
                ggsave(filename = file_name_facet)
            }
            
            # Initialize a list to store p-values for timepoints
            p_values_list <- list(TP0 = NA, TP1 = NA, TP2 = NA, TP3 = NA, TP4 = NA)
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            timepoints <- c("TP0", "TP1", "TP2", "TP3", "TP4")
            for (tp in timepoints) {
                male_values <- df_subpop_clean$mean_Z[df_subpop_clean$Timepoint == tp & df_subpop_clean$Sex == "Male"]
                female_values <- df_subpop_clean$mean_Z[df_subpop_clean$Timepoint == tp & df_subpop_clean$Sex == "Female"]
                
                # Perform Wilcoxon test for the comparison of Male vs Female at each timepoint
                if (length(male_values) > 0 && length(female_values) > 0) {
                    wilcox_result <- wilcox.test(male_values, female_values)
                    p_values_list[[tp]] <- wilcox_result$p.value
                }
            }
            
            # Add results for Male vs Female for all timepoints as a single row
            results <- rbind(results, data.frame(WilcoxComparison = "Male vs Female", Gene = gene, Subpopulation = sup, 
                                                 TP0_P_Value = p_values_list[["TP0"]],
                                                 TP1_P_Value = p_values_list[["TP1"]],
                                                 TP2_P_Value = p_values_list[["TP2"]],
                                                 TP3_P_Value = p_values_list[["TP3"]],
                                                 TP4_P_Value = p_values_list[["TP4"]]))
            
            
        }, error = function(e) {
            cat("Error processing", gene, ":", e$message, "\n")
        })
    }
}
Timecourse_Wilcox_Sex_Z <- results
Timecourse_Wilcox_Sex_Z$TP0_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP0_P_Value, get_significance)
Timecourse_Wilcox_Sex_Z$TP1_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP1_P_Value, get_significance)
Timecourse_Wilcox_Sex_Z$TP2_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP2_P_Value, get_significance)
Timecourse_Wilcox_Sex_Z$TP3_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP3_P_Value, get_significance)
Timecourse_Wilcox_Sex_Z$TP4_Significance <- sapply(Timecourse_Wilcox_Sex_Z$TP4_P_Value, get_significance)
write.csv(Timecourse_Wilcox_Sex_Z , file = "Timecourse_Wilcox_Sex_Z_Results.csv", row.names = FALSE)
    
 #### Healthy vs MILD vs MODERATE  (Mean ± SEM)* -------------------------------------
folder <- "Mean_SEM_Category_Plots"
createFolder(folder)  # Ensure the folder exists
results <- data.frame(Category = character(), Gene = character(), Subpopulation = character(), 
                      P_Value = numeric(), TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                      TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (i in 1:length(unique(data_rE_mean_matched$Gene))) {
    gene <- unique(data_rE_mean_matched$Gene)[i]
    df_gene <- filter(data_rE_mean_matched, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_Z), !is.na(Subpopulation), !is.na(Category))
    
    for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- unique(df_gene$Subpopulation)[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
        
        # Remove outliers
        df_subpop_clean <- df_subpop[!df_subpop$mean_Z %in% boxplot.stats(df_subpop$mean_Z)$out, ]
        
        if (nrow(df_subpop_clean) <= 0) next  # Skip if no rows remain
        
        tryCatch({
            title_facet <- paste(gene, "expression in", sup, "monocytes (Zscored dCT)", sep=" ")
            
            # Line plot comparing Male vs Female with mean ± SEM across Timepoints
            ggline(subset(df_subpop_clean, !is.na(Category)), 
                   x = "Timepoint", 
                   y = "mean_Z", 
                   color = "Category", 
                   add = "mean_se",  # Add mean and standard error
                   size = 1.2) +  
                labs(y = paste(gene, "expression (Z-scored from CT - CT Sample Median)"), 
                     title = title_facet) +
                scale_color_manual(values = c("MINOR" = "#A67C00", "MODERATE" = "#700606")) +
                theme_bw() +  # Apply white background
                theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid', color = "gray80"), 
                      panel.grid.minor = element_blank(), 
                      panel.border = element_blank(),
                      axis.line = element_line(color = "black")) +
                stat_compare_means(method = "wilcox.test",  # or "t.test" 
                                   aes(group = Category),
                                   label = "p.signif",  
                                   size = 5) 
            
            if (plot_save == "y") {
                file_name_facet <- paste0(folder, "/", title_facet, "_Mean_SEM.png")
                ggsave(filename = file_name_facet)
            }
           
            # Initialize a list to store p-values for timepoints
            p_values_list <- list(TP0 = NA, TP1 = NA, TP2 = NA, TP3 = NA, TP4 = NA)
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            timepoints <- c("TP0", "TP1", "TP2", "TP3", "TP4")
            for (tp in timepoints) {
                minor_values <- df_subpop_clean$mean_Z[df_subpop_clean$Timepoint == tp & df_subpop_clean$Category == "MINOR"]
                moderate_values <- df_subpop_clean$mean_Z[df_subpop_clean$Timepoint == tp & df_subpop_clean$Category == "MODERATE"]
                
                # Perform Wilcoxon test for the comparison of Minor vs Moderate at each timepoint
                if (length(minor_values) > 0 && length(moderate_values) > 0) {
                    wilcox_result <- wilcox.test(minor_values, moderate_values)
                    p_values_list[[tp]] <- wilcox_result$p.value
                }
            }
            
            # Add results for Male vs Female for all timepoints as a single row
            results <- rbind(results, data.frame(WilcoxComparison = "Minor vs Moderate", Gene = gene, Subpopulation = sup, 
                                                 TP0_P_Value = p_values_list[["TP0"]],
                                                 TP1_P_Value = p_values_list[["TP1"]],
                                                 TP2_P_Value = p_values_list[["TP2"]],
                                                 TP3_P_Value = p_values_list[["TP3"]],
                                                 TP4_P_Value = p_values_list[["TP4"]]))
            
        }, error = function(e) {
            cat("Error processing", gene, ":", e$message, "\n")
        })
    }
}

Timecourse_Wilcox_Category_Z <- results
Timecourse_Wilcox_Category_Z$TP0_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP0_P_Value, get_significance)
Timecourse_Wilcox_Category_Z$TP1_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP1_P_Value, get_significance)
Timecourse_Wilcox_Category_Z$TP2_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP2_P_Value, get_significance)
Timecourse_Wilcox_Category_Z$TP3_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP3_P_Value, get_significance)
Timecourse_Wilcox_Category_Z$TP4_Significance <- sapply(Timecourse_Wilcox_Category_Z$TP4_P_Value, get_significance)
write.csv(Timecourse_Wilcox_Category_Z , file = "Timecourse_Wilcox_Category_Z_Results.csv", row.names = FALSE)

# *Regression* --------------------------------------------------------
#### Pearson Gene vs Age --------------------------------------------------------

createFolder("LinReg_Results")

# Pearson's Correlation uses linear relationship to correlate the data
Age_correlation_rE <- data.frame(Subpopulation = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

Age_correlation_Sex_rE <- data.frame(Subpopulation = character(), Sex = character(), Timepoint = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

for (gene in unique(dataFINALmean$Gene)) {
    df <- dataFINALmean %>% filter(Gene == gene)
    df <- df %>% filter(!is.na(Timepoint), !is.na(mean_rE), !is.na(Subpopulation), !is.na(Sex))
  for (subpop in unique(df$Subpopulation)) {
    df_subpop <- df %>% filter(Subpopulation == subpop)
    for (tp in unique(df_subpop$Timepoint)) {
      df_tp <- df_subpop %>% filter(Timepoint== tp)
    if (nrow(df_tp) <= 0) {
      next
    }
    tryCatch({
      x <- cor.test(df_tp$mean_rE, df_tp$Age)
      row <- data.frame(
        Subpopulation = subpop,
        Timepoint = tp,
        Gene = gene,
        GeneID = df_tp$GeneID[1],
        Gene_Group = df_tp$Gene_Group[1],
        N = nrow(df_tp),
        p.value = x$p.value,
        Estimate = x$estimate,
        Est_CI_Lower = x$conf.int[1],
        Est_CI_Upper = x$conf.int[2],
        Coefficient = NA,
        Coef_CI_Lower = NA,
        Coef_CI_Upper = NA  
      )
      
      # Linear regression to get the coefficient and its CI
      lm_result <- lm(mean_rE ~ Age, data = df_tp)
      coef_x <- summary(lm_result)$coefficients[2, 1]
      ci <- confint(lm_result)[2, ]
      
      # Update the row with the coefficient and its CI
      row$Coefficient <- coef_x
      row$Coef_CI_Lower <- ci[1]
      row$Coef_CI_Upper <- ci[2]
      
      Age_correlation_rE <- rbind(Age_correlation_rE, row)
    }, error = function(e) {
      cat("Error in Pearson's Correlation for gene", genes[i], "& Subpopulation with Age", subpop, "- skipping this comparison.\n")
    })
      for (sex in unique(df_tp$Sex)) {
      df_sex <- df_tp %>% filter(Sex == sex)
    if (nrow(df_sex) <= 0) {
      next
    }
    tryCatch({
        x <- cor.test(df_sex$mean_rE, df_sex$Age)
        row <- data.frame(
          Subpopulation = subpop,
          Sex = sex,
          Timepoint = tp,
          Gene = gene,
          GeneID = df_sex$GeneID[1],
          Gene_Group = df_sex$Gene_Group[1],
          N = nrow(df_sex),
          p.value = x$p.value,
          Estimate = x$estimate,
          Est_CI_Lower = x$conf.int[1],
          Est_CI_Upper = x$conf.int[2],
          Coefficient = NA, 
          Coef_CI_Lower = NA, 
          Coef_CI_Upper = NA  
        )
        
        # Linear regression to get the coefficient and its CI
        lm_result_sex <- lm(mean_rE ~ Age, data = df_sex)
        coef_x <- summary(lm_result_sex)$coefficients[2, 1]
        ci <- confint(lm_result_sex)[2, ]
        
        # Update the row with the coefficient and its CI
        row$Coefficient <- coef_x
        row$Coef_CI_Lower <- ci[1]
        row$Coef_CI_Upper <- ci[2]
        
        Age_correlation_Sex_rE <- rbind(Age_correlation_Sex_rE, row)
      }, error = function(e) {
        cat("Error for ", gene, ", ", tp, sex, "&", subpop, "- skipping this comparison.\n")
      }) 
      }      
    }
  }
  }

file_name <- "LinReg_Results/Pearson's Correlation_Age_rE.csv"
write.csv(Age_correlation_rE, file_name, row.names = FALSE)

file_name_s <- "LinReg_Results/Pearson's Correlation_Age_Sex_rE.csv"
write.csv(Age_correlation_Sex_rE, file_name_s, row.names = FALSE)

#### Extract significant LinRegs --------------------------------------------------------

Age_correlation_rE$Sex <- "All"
All_Age_Pearson_rE <- rbind(Age_correlation_Sex_rE, Age_correlation_rE)
Age_sigGenes_Pearson <- All_Age_Pearson_rE %>% filter(All_Age_Pearson_rE$p.value <0.05)
write.csv(Age_sigGenes_Pearson, "LinReg_Results/Age_sigGenes_Pearson.csv", row.names = FALSE)

#### Plot significant LinRegs (F/M) --------------------------------------------------------
createFolder("LinReg_Sex")

i <- 1
j <- nrow(Age_sigGenes_Pearson)
for (i in 1:j) {
    sub <- Age_sigGenes_Pearson$Subpopulation[i]
    sig <- Age_sigGenes_Pearson$Gene[i]
    tp <- Age_sigGenes_Pearson$Timepoint[i]
    df <- filter(dataFINALmean, Gene == sig)
    df <- filter(df, Subpopulation == sub)
    df <- filter(df, Timepoint == tp)
    s <- Age_sigGenes_Pearson$Sex[i]
    if (s=="All"){
        titel <- paste(sig, " expression in ", sub, " monocytes at ", tp," (LinReg, rE, combined)", sep="")
file_name <- paste("LinReg_Sex/",titel, ".png", sep = "")
    ggplot(data = df, aes(x = Age, y = mean_rE)) +
  geom_smooth(method = "glm", color = "black") +
  geom_point(aes(color = Sex), size = 2) +
  ggtitle(titel) +
  xlab("Age") +
  ylab("Gene expression relative to B2M") +
  scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2")) +
  theme(text = element_text(size=14)) +
    theme(
        panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = "black"),
        text = element_text(color = "black"),
        panel.grid.major = element_line(color = "gray", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
    ggsave(filename=file_name)
    #dev.off()
    i = i +1
    }
    
    else {
    titel <- paste(sig, " expression in ",sub," monocytes at ", tp," (",s,") (LinReg, rE)", sep="")
file_name <- paste("LinReg_Sex/",titel, ".png", sep = "")
#png(filename=file_name, height=750, width=750)
ggplot(data = df, aes(x = Age, y = mean_rE, colour = Sex)) +
    geom_smooth(method = glm) + 
    geom_point(size = 2) +
    ggtitle(titel) +
    xlab("Age") + 
    ylab("Gene expression relative to B2M") +
    scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2")) +
    theme_bw(base_size = 14) +
    theme(
        panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = "black"),
        text = element_text(color = "black"),
        panel.grid.major = element_line(color = "gray", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
    ggsave(filename=file_name)
    #dev.off()
    }
}

## Plot significant LinRegs (TPs) --------------------------------------------------------

createFolder("LinReg_TP")

j <- nrow(Age_sigGenes_Pearson)
for (i in 1:j) {
    sub <- Age_sigGenes_Pearson$Subpopulation[i]
    sig <- Age_sigGenes_Pearson$Gene[i]
    tp <- Age_sigGenes_Pearson$Timepoint[i]
    df <- filter(dataFINALmean, Gene == sig)
    df <- filter(df, Subpopulation == sub)
    df <- filter(df, Timepoint == tp)
    s <- Age_sigGenes_Pearson$Sex[i]
    if (s=="All"){
        titel <- paste(sig, " expression in ", sub, " monocytes at ", tp," (LinReg, rE, combined)", sep="")
file_name <- paste("LinReg_TP/",titel, ".png", sep = "")
    ggplot(data = df, aes(x = Age, y = mean_rE, colour = Timepoint)) +
  geom_smooth(method = "glm") +
  geom_point(size = 2) +
  ggtitle(titel) +
  xlab("Age") +
  ylab("Gene expression relative to B2M") +
  scale_color_manual(values = my_colors) +
  theme(text = element_text(size=14)) +
    theme(
        panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = "black"),
        text = element_text(color = "black"),
        panel.grid.major = element_line(color = "gray", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
    ggsave(filename=file_name)
    #dev.off()
    }
    
    else {
        df <- filter(df, Sex == s)
    titel <- paste(sig, " expression in ",sub," monocytes at ", tp," (",s,") (LinReg, rE)", sep="")
file_name <- paste("LinReg_TP/",titel, ".png", sep = "")
#png(filename=file_name, height=750, width=750)
ggplot(data = df, aes(x = Age, y = mean_rE, colour = Timepoint)) +
    geom_smooth(method = glm) + 
    geom_point(size = 2) +
    ggtitle(titel) +
    xlab("Age") + 
    ylab("Gene expression relative to B2M") +
    scale_color_manual(values = my_colors) +
    theme_bw(base_size = 14) +
    theme(
        panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = "black"),
        text = element_text(color = "black"),
        panel.grid.major = element_line(color = "gray", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
    ggsave(filename=file_name)
    #dev.off()
    }
}

## Plot significant LinRegs (all TPs) --------------------------------------------------------

j <- nrow(Age_sigGenes_Pearson)
for (i in 1:j) {
    sub <- Age_sigGenes_Pearson$Subpopulation[i]
    sig <- Age_sigGenes_Pearson$Gene[i]
    tp <- Age_sigGenes_Pearson$Timepoint[i]
    df <- filter(dataFINALmean, Gene == sig)
    df <- filter(df, Subpopulation == sub)
    s <- Age_sigGenes_Pearson$Sex[i]
    if (s=="All"){
        titel <- paste(sig, " expression in ", sub, " monocytes at ", tp," (LinReg, rE, combined)", sep="")
file_name <- paste("LinReg_TP_all/",titel, ".png", sep = "")
    ggplot(data = df, aes(x = Age, y = mean_rE, colour = Timepoint)) +
  geom_smooth(method = "glm") +
  geom_point(size = 2) +
  ggtitle(titel) +
  xlab("Age") +
  ylab("Gene expression relative to B2M") +
  scale_color_manual(values = my_colors) +
  theme(text = element_text(size=14)) +
    theme(
        panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = "black"),
        text = element_text(color = "black"),
        panel.grid.major = element_line(color = "gray", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
    ggsave(filename=file_name)
    #dev.off()
    }
    
    else {
        df <- filter(df, Sex == s)
    titel <- paste(sig, " expression in ",sub," monocytes at ", tp," (",s,") (LinReg, rE)", sep="")
file_name <- paste("LinReg_TP_all/",titel, ".png", sep = "")
#png(filename=file_name, height=750, width=750)
ggplot(data = df, aes(x = Age, y = mean_rE, colour = Timepoint)) +
    geom_smooth(method = glm) + 
    geom_point(size = 2) +
    ggtitle(titel) +
    xlab("Age") + 
    ylab("Gene expression relative to B2M") +
    scale_color_manual(values = my_colors) +
    theme_bw(base_size = 14) +
    theme(
        panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = "black"),
        text = element_text(color = "black"),
        panel.grid.major = element_line(color = "gray", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
    ggsave(filename=file_name)
    #dev.off()
    }
}
