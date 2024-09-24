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

thematic_on(bg = "white", fg = "black", accent = "blue")

### *Functions* ---------------------------------------------------------------------

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

getMeTheShitSamples = function( data, gene, cutoff ){
    
    data$Sample[which(data$Gene == gene & data$CT_Value > cutoff
                      #, na.rm = TRUE
    )]
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

## *Data* ---------------------------------------------------------------------

## Loading of the Data -----------------------------------------
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

## trimming of data ------------------------------------------------------
data <- cbind.data.frame(Sample=raw_data$Name,
                         Gene=raw_data$Name.1,
                         Value=raw_data$Value,
                         dCT_Value=raw_data$Value.1,
                         Call=raw_data$Call.1)

# evaluate PASS FAIL------------------------------------------------
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

## Information on the data ---------------------------------------------------------------------
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

## Expand data --------------------------------------------------------
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

## Metadata ---------------------------------------------------------------
metadataP <- read_excel("Metadata_PatientID.xlsx")
Metadata_GeneID <- read_excel("Metadata_GeneID.xlsx")
Plate_info <- read_excel("240410_PlateID.xlsx")

# Combine the data frames based on common values in the "id" column
# --> here i loose all the test because they are not in the metadata
dataAll <- merge(dataEXPAND, metadataP, by = "SampleID", all.x = TRUE)

#relable Male & Female
dataAll$Sex[grep("M",dataAll$Sex)]  <- "Male"
dataAll$Sex[grep("F",dataAll$Sex)]  <- "Female"

dataAll <- dataAll %>% filter(!is.na(Age), !is.na(SampleID), !is.na(Category), !is.na(Sex))

## failed Genes --------------------------------------------------------
fails <- dataAll %>%
  group_by(Subpopulation, Gene, Sex) %>%
  summarize(
    N = n(),
    Count_Above_35_CT <- sum(CT_Value >= 35, na.rm = TRUE),
    Per_Above_35_CT = mean(CT_Value >= 35, na.rm = TRUE) * 100,
    Count_Above_35_dCT = sum(dCT_Value >= 35, na.rm = TRUE),
    Per_Above_35_dCT = mean(dCT_Value >= 35, na.rm = TRUE) * 100
  )

failed_genes <- fails %>% filter(Per_Above_35_CT >= 81) %>% select(Gene)
wrong_threshold <- c("CLEC7A", "S100A8", "S100A9")
exclude_genes  <- unique(c("Xeno",failed_genes$Gene, wrong_threshold))
paste("The gene",unique(failed_genes$Gene), "was excluded, as more than 80% of the reads failed.")
paste("The gene",wrong_threshold, "was excluded, as the global threshold setting did not allow proper evaluation.")

# remove the failed genes:
dataWORKING <- dataAll %>%  filter(!Gene %in% exclude_genes)

## failed Samples --------------------------------------------------------
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



## Prepare Matrix for imputation ---------------------------------------------------------------------
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

# Imputation --------------------------------------------------------------------------------------------
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

# *Normalize the data vs the median CT_Value if each Sample* ----------------

# Calculate the CT median
dataFINAL$CTmedian <- NA
dataFINAL <- dataFINAL %>%
    group_by(Sample) %>%
    mutate(CTmedian = median(CT, na.rm = TRUE))

# normalize all CT value to the Media CT value of each Sample
dataFINAL$dCT_median <- dataFINAL$CT - dataFINAL$CTmedian

# check the variance again after the imputation
filtered_data <- dataFINAL %>%
    filter(Gene %in% c("B2M", "ACTB", "GAPDH"))

# Calculate the CV for each group
variance_table <- filtered_data %>%
    group_by(Gene, Timepoint, Subpopulation) %>%
    summarize(CV = sd(CT, na.rm = TRUE) / mean(CT, na.rm = TRUE)) %>%
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


# calculated the relative expression
dataFINAL$rE_Value <- 2^((-1)*dataFINAL$dCT_median)

dataFINAL$rE_Value <- ifelse(dataFINAL$CT>=35, yes=0, no=dataFINAL$rE_Value)
dataFINAL$rE_Value <- ifelse(dataFINAL$dCT_median>=35, yes=0, no=dataFINAL$rE_Value) 

# CHeck but if the median is above 30 the run was probably not successfull
dataFINAL$rE_Value <- ifelse(dataFINAL$CTmedian>=30, yes=0, no=dataFINAL$rE_Value)
    
# *FACS data* ----------------------------------------------------------------------------------------------------
FACSdata <- read_excel("240410_FACSdata_StrokeControl.xls")
FACSdata <- as.data.frame(FACSdata, col_types = c("text", 
    "text", "numeric", "text", "numeric", 
    "text", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric"))

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

# Gene Info  ----------------
## Gene Metadata & other edits:
Gene_Info <- data.frame(Gene = Metadata_GeneID$Gene, GeneID = Metadata_GeneID$GeneID, Gene_Group = Metadata_GeneID$Gene_Group)

dataFINAL <- merge(dataFINAL, Gene_Info, by = "Gene")

# List of Patients (2x 15))
patients <- unique(dataFINAL$SampleID)
no_patients <- length(patients)

# Days Post-Stroke data ----------------
# DaysPS <- read_excel("Metadata_DaysPS.xls")
dataFINAL$DaysPS <- as.numeric(NA)
dataFINAL$DaysPS[grep("TP1",dataFINAL$Timepoint)]  <- 0
dataFINAL$DaysPS[grep("TP1",dataFINAL$Timepoint)]  <- 1
dataFINAL$DaysPS[grep("TP2",dataFINAL$Timepoint)]  <- 4
dataFINAL$DaysPS[grep("TP3",dataFINAL$Timepoint)]  <- 30
dataFINAL$DaysPS[grep("TP4",dataFINAL$Timepoint)]  <- 90
dataFINAL$DaysPS[grep("TP5",dataFINAL$Timepoint)]  <- 360
#unique(dataFINAL$DaysPS)

## create mean values for the ANOVA & so
dataFINALmean <- dataFINAL %>%
        group_by(SampleID, Sex, Age, Category, Timepoint, DaysPS,
                 Subpopulation, Gene, GeneID) %>%
        summarise(mean_CT = mean(CT), sd_CT = sd(CT),
                  mean_dCT = mean(dCT_median), sd_dCT = sd(dCT_median),
                  mean_rE = mean(rE_Value), sd_rE = sd(rE_Value)
                  )%>%
  ungroup()

## Change working directory so that everything is getting saved in the right place ----------------
today <- Sys.Date()
output_location <- paste(today,"_Stroke_Results", sep="")
setwd(paste("/Users/ju5263ta/Github/Monocytes/Data/",output_location, "/",sep=""))
getwd()

write_xlsx(fails, "Fail_rates.xlsx")
write_xlsx(dataFINAL, "dataFINAL.xlsx")
write_xlsx(dataFINALmean, "dataFINALmean.xlsx")

# *ANOVA* ----------------------------------------------------------------------------------------------------
# Adjust to yes if you want plots to be saved or n if not
anova_plot_save <- "y"

# rE with Baseline - add adjusted p-Value ? --------------------------------------------------
createFolder("ANOVA_results")
folder <- "ANOVAwilcox_rE exOutliers"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric(),
                      TP1_P_Value = numeric(), TP2_P_Value = numeric(), TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (j in 1:length(unique(dataFINALmean$Gene))) {
    gene <- dataFINALmean$Gene[j]
    df_gene <- filter(dataFINALmean, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Category), !is.na(mean_rE), !is.na(Subpopulation), !is.na(Timepoint))
    
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
        
        # Perform ANOVA
        anova_result <- aov(mean_rE ~ Timepoint, data = df)
        
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
        
        if (anova_plot_save == "y") {
        my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 
        
        ggboxplot(df, 
                  x = "Timepoint", 
                  y = "mean_rE", 
                  color = "Timepoint", 
                  add = "jitter", 
                  legend = "none", 
                  ylab = paste(gene, "expression relative to B2M"), 
                  width = 0.8, 
                  add.params = list(size = 1, alpha = 0.5)) +  
            geom_hline(yintercept = median_TP0, linetype = 2) + 
            stat_compare_means(method = "anova", label.y = max(df$mean_rE)) +        
            stat_compare_means(label = "p.signif", 
                               method = "wilcox", 
                               ref.group = "TP0", 
                               hide.ns = TRUE, 
                               label.y = max(df$mean_rE)) +
            scale_color_manual(values = my_colors)
        
        file_name <- file.path(folder, paste(gene,"_", sup, "_rE_ANOVAwilcox.png", sep = ""))
        ggsave(filename = file_name)
        }
    }
}
ANOVA_rE <- results
ANOVA_rE$TP1_Significance <- sapply(ANOVA_rE$TP1_P_Value, get_significance)
ANOVA_rE$TP2_Significance <- sapply(ANOVA_rE$TP2_P_Value, get_significance)
ANOVA_rE$TP3_Significance <- sapply(ANOVA_rE$TP3_P_Value, get_significance)
ANOVA_rE$TP4_Significance <- sapply(ANOVA_rE$TP4_P_Value, get_significance)
write.csv(ANOVA_rE, file = "ANOVA_results/ANOVA_Results_rE.csv", row.names = FALSE)

# dCT with Baseline - with Outlier exclusion --------------------------------------------------
folder <- "ANOVAwilcox_dCT exOutliers"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric(),
                      TP1_P_Value = numeric(), TP2_P_Value = numeric(), TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (j in 1:length(unique(dataFINALmean$Gene))) {
    gene <- dataFINALmean$Gene[j]
    df_gene <- filter(dataFINALmean, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Category), !is.na(mean_dCT), !is.na(Subpopulation), !is.na(Timepoint))
    
    for (i in 1:length(unique(df_gene$Subpopulation))) {
        sup <- df_gene$Subpopulation[i]
        df <- filter(df_gene, Subpopulation == sup)
        if (nrow(df) <= 0 || gene %in% exclude_genes) {
            cat("Skipping:", gene, "for subpopulation:", sup, "due to insufficient data or in excluded genes list\n")
            next  # Skip to the next iteration of the loop
        }
        
        df <- df[!df$mean_dCT %in% boxplot.stats(df$mean_dCT)$out, ]
        df$Timepoint <- factor(as.vector(df$Timepoint))
        baseline <- filter(df, Timepoint == "TP0")
        median_TP0 <- median(baseline$mean_dCT)
        
        # Perform ANOVA
        anova_result <- aov(mean_dCT ~ Timepoint, data = df)
        
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
        
        if (anova_plot_save == "y") {
        # Plot the boxplot and comparison without outliers
        ggboxplot(df, 
                  x = "Timepoint", 
                  y = "mean_dCT", 
                  color = "Timepoint", 
                  add = "jitter", 
                  legend = "none", 
                  ylab = paste(gene, "expression [dCT to B2M]"), 
                  width = 0.8, 
                  add.params = list(size = 1, alpha = 0.5)) +  
            geom_hline(yintercept = median_TP0, linetype = 2) + 
            stat_compare_means(method = "anova", label.y = max(df$mean_dCT)) +        
            stat_compare_means(label = "p.signif", 
                               method = "wilcox", 
                               ref.group = "TP0", 
                               hide.ns = TRUE, 
                               label.y = max(df$mean_dCT)) +
            scale_color_manual(values = my_colors)
        
        file_name <- file.path(folder, paste(gene,"_", sup, "_dCT_ANOVAwilcox.png", sep = ""))
        ggsave(filename = file_name)
        }
    }
}

ANOVA_dCT <- results
ANOVA_dCT$TP1_Significance <- sapply(ANOVA_dCT$TP1_P_Value, get_significance)
ANOVA_dCT$TP2_Significance <- sapply(ANOVA_dCT$TP2_P_Value, get_significance)
ANOVA_dCT$TP3_Significance <- sapply(ANOVA_dCT$TP3_P_Value, get_significance)
ANOVA_dCT$TP4_Significance <- sapply(ANOVA_dCT$TP4_P_Value, get_significance)
write.csv(ANOVA_dCT, file  = "ANOVA_results//ANOVA_Results_dCT.csv", row.names = FALSE) 


## Male vs Female Analysis --------------------------------------------------
## you can try to adjjust the loop that it also only plots the facet comparisons upon request
folder <- "ANOVA_rE_Sex"
createFolder(folder)  # Ensure the folder exists
results <- data.frame(Sex = character(), Gene = character(), Subpopulation = character(), 
                      P_Value = numeric(), TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                      TP3_P_Value = numeric(), TP4_P_Value = numeric())

for (i in 1:length(unique(dataFINALmean$Gene))) {
    gene <- unique(dataFINALmean$Gene)[i]
    df_gene <- filter(dataFINALmean, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_rE), !is.na(Subpopulation), !is.na(Sex))
    
    for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- unique(df_gene$Subpopulation)[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
        
        # Remove outliers
        df_subpop_clean <- df_subpop[!df_subpop$mean_rE %in% boxplot.stats(df_subpop$mean_rE)$out, ]
        
        if (nrow(df_subpop_clean) <= 0) next  # Skip if no rows remain
        
        tryCatch({
            title_facet <- paste(gene, "expression in", sup, "monocytes (LinReg, rE)", sep=" ")
            
            # Faceted boxplot: comparing Male vs Female across Timepoints
            ggboxplot(df_subpop_clean, 
                      x = "Sex", 
                      y = "mean_rE", 
                      color = "Sex",  
                      title = title_facet,
                      add = "jitter", 
                      facet.by = "Timepoint",
                      legend = "none", 
                      ylab = paste(gene, "expression relative to B2M"), 
                      width = 0.8, 
                      add.params = list(size = 1, alpha = 0.5)) +  
                stat_compare_means(method = "anova", label.y = max(df_subpop_clean$mean_rE)) +        
                stat_compare_means(label = "p.signif", 
                                   method = "wilcox", 
                                   ref.group = "Male", 
                                   hide.ns = TRUE, 
                                   label.y = max(df_subpop_clean$mean_rE)) +
                scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2"))
            
            if (anova_plot_save == "y") {
            file_name_facet <- paste0(folder, "/", title_facet, "_Facets.png")
            ggsave(filename = file_name_facet)
            }
            
            # Male analysis
            male <- filter(df_subpop_clean, Sex == "Male")
            title_male <- paste(gene, "expression in", sup, "(Male) monocytes (LinReg, rE)", sep = "")
            male_clean <- male[!male$mean_rE %in% boxplot.stats(male$mean_rE)$out, ]
            
            baseline_male <- filter(male_clean, Timepoint == "TP0")
            median_TP0_male <- median(baseline_male$mean_rE)
            
            # ANOVA for Male
            anova_m <- aov(mean_rE ~ Timepoint, data = male_clean)
            p_value_m <- summary(anova_m)[[1]]$"Pr(>F)"[1]
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            timepoints <- c("TP1", "TP2", "TP3", "TP4")
            p_values_tp_m <- sapply(timepoints, function(tp) {
                if (sum(male_clean$Timepoint == tp) > 0) {
                    tryCatch({
                        wilcox.test(male_clean$mean_rE[male_clean$Timepoint == "TP0"],
                                    male_clean$mean_rE[male_clean$Timepoint == tp])$p.value
                    }, error = function(e) NA)
                } else {
                    NA
                }
            })
            
            results <- rbind(results, data.frame(Sex = "Male", Gene = gene, Subpopulation = sup, 
                                                 P_Value = p_value_m, 
                                                 TP1_P_Value = p_values_tp_m[1], 
                                                 TP2_P_Value = p_values_tp_m[2], 
                                                 TP3_P_Value = p_values_tp_m[3], 
                                                 TP4_P_Value = p_values_tp_m[4]))
            
            if (anova_plot_save == "y") {
                ggboxplot(male_clean, 
                          x = "Timepoint", 
                          y = "mean_rE", 
                          color = "Timepoint", 
                          title = title_male,
                          add = "jitter", 
                          legend = "none", 
                          ylab = paste(gene, "expression relative to B2M"), 
                          width = 0.6, 
                          add.params = list(size = 1, alpha = 0.5)) +  
                    geom_hline(yintercept = median_TP0_male, linetype = 2) + 
                    stat_compare_means(method = "anova", label.y = max(male_clean$mean_rE)) +        
                    stat_compare_means(label = "p.signif", 
                                       method = "wilcox", 
                                       ref.group = "TP0", 
                                       hide.ns = TRUE, 
                                       label.y = max(male_clean$mean_rE)) +
                    scale_color_manual(values = my_colors)
                
                file_name_male <- paste0(folder, "/", title_male, ".png")
                ggsave(filename = file_name_male)
            }
            
            # Female analysis
            female <- filter(df_subpop_clean, Sex == "Female")
            title_female <- paste(gene, "expression in", sup, "(Female) monocytes (LinReg, rE)", sep = "")
            female_clean <- female[!female$mean_rE %in% boxplot.stats(female$mean_rE)$out, ]
            
            baseline_female <- filter(female_clean, Timepoint == "TP0")
            median_TP0_female <- median(baseline_female$mean_rE)
            
            # ANOVA for Female
            anova_f <- aov(mean_rE ~ Timepoint, data = female_clean)
            p_value_f <- summary(anova_f)[[1]]$"Pr(>F)"[1]
            
            # Extract p-values for each timepoint comparison (pairwise Wilcoxon)
            p_values_tp_f <- sapply(timepoints, function(tp) {
                if (sum(female_clean$Timepoint == tp) > 0) {
                    tryCatch({
                        wilcox.test(female_clean$mean_rE[female_clean$Timepoint == "TP0"],
                                    female_clean$mean_rE[female_clean$Timepoint == tp])$p.value
                    }, error = function(e) NA)
                } else {
                    NA
                }
            })
            
            results <- rbind(results, data.frame(Sex = "Female", Gene = gene, Subpopulation = sup, 
                                                 P_Value = p_value_f, 
                                                 TP1_P_Value = p_values_tp_f[1], 
                                                 TP2_P_Value = p_values_tp_f[2], 
                                                 TP3_P_Value = p_values_tp_f[3], 
                                                 TP4_P_Value = p_values_tp_f[4]))
            
            if (anova_plot_save == "y") {
                ggboxplot(female_clean, 
                          x = "Timepoint", 
                          y = "mean_rE", 
                          color = "Timepoint", 
                          title = title_female,
                          add = "jitter", 
                          legend = "none", 
                          ylab = paste(gene, "expression relative to B2M"), 
                          width = 0.6, 
                          add.params = list(size = 1, alpha = 0.5)) +  
                    geom_hline(yintercept = median_TP0_female, linetype = 2) + 
                    stat_compare_means(method = "anova", label.y = max(female_clean$mean_rE)) +        
                    stat_compare_means(label = "p.signif", 
                                       method = "wilcox", 
                                       ref.group = "TP0", 
                                       hide.ns = TRUE, 
                                       label.y = max(female_clean$mean_rE)) +
                    scale_color_manual(values = my_colors)
                
                file_name_female <- paste0(folder, "/", title_female, ".png")
                ggsave(filename = file_name_female)
            }
        }, error = function(e) {
            cat("Error processing", gene, ":", e$message, "\n")
        })
    }
}
results$TP1_Significance <- sapply(results$TP1_P_Value, get_significance)
results$TP2_Significance <- sapply(results$TP2_P_Value, get_significance)
results$TP3_Significance <- sapply(results$TP3_P_Value, get_significance)
results$TP4_Significance <- sapply(results$TP4_P_Value, get_significance)
ANOVA_Sex <- results 
# Save ANOVA results with significance levels
write.csv(ANOVA_Sex, file = "ANOVA_results/ANOVA_Results_Sex.csv", row.names = FALSE)

## Healthy vs MILD vs MODERATE --------------------------------------------------------
folder <- "ANOVA_rE_Category"
createFolder(folder)  # Ensure the folder exists
results <- data.frame(Category = character(), Gene = character(), Subpopulation = character(), 
                      P_Value = numeric(), TP1_P_Value = numeric(), TP2_P_Value = numeric(), 
                      TP3_P_Value = numeric(), TP4_P_Value = numeric())

i <- 1
for (i in 1:length(unique(dataFINALmean$Gene))) {
    gene <- unique(dataFINALmean$Gene)[i]  # Use unique to avoid duplicates
    df_gene <- filter(dataFINALmean, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_rE), !is.na(Subpopulation), !is.na(Category))
    
    j <- 1
    for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- unique(df_gene$Subpopulation)[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
        
        # Keep TP0 data in both categories
        df_minor <- df_subpop %>% filter(Category != "MODERATE" | Timepoint == "TP0")
        df_moderate <- df_subpop %>% filter(Category != "MINOR" | Timepoint == "TP0")
        
        # Remove outliers for MINOR and MODERATE categories
        df_minor <- df_minor[!df_minor$mean_rE %in% boxplot.stats(df_minor$mean_rE)$out, ]
        df_moderate <- df_moderate[!df_moderate$mean_rE %in% boxplot.stats(df_moderate$mean_rE)$out, ]
        
        # Ensure Timepoint is a factor
        df_minor$Timepoint <- factor(as.vector(df_minor$Timepoint))
        df_moderate$Timepoint <- factor(as.vector(df_moderate$Timepoint))
        
        # Get baseline median for MINOR and MODERATE
        baseline_mi <- filter(df_minor, Timepoint == "TP0")
        baseline_mo <- filter(df_moderate, Timepoint == "TP0")
        
        median_TP0_mi <- median(baseline_mi$mean_rE, na.rm = TRUE)
        median_TP0_mo <- median(baseline_mo$mean_rE, na.rm = TRUE)
        
        # ANOVA
        if (nrow(df_minor) > 1) {  # Check if there are enough rows for ANOVA
            anova_mi <- aov(mean_rE ~ Timepoint, data = df_minor)
            p_value_mi <- summary(anova_mi)[[1]]$"Pr(>F)"[1]
        } else {
            p_value_mi <- NA  # Handle cases where ANOVA cannot be performed
        }
        
        if (nrow(df_moderate) > 1) {  # Check if there are enough rows for ANOVA
            anova_mo <- aov(mean_rE ~ Timepoint, data = df_moderate)
            p_value_mo <- summary(anova_mo)[[1]]$"Pr(>F)"[1]
        } else {
            p_value_mo <- NA  # Handle cases where ANOVA cannot be performed
        }
        
        # Store the pairwise Wilcoxon test results for each timepoint
        results <- rbind(results, data.frame(Category = "MINOR", Gene = gene, Subpopulation = sup, 
                                             P_Value = p_value_mi, TP1_P_Value = p_values_mi[1], TP2_P_Value = p_values_mi[2], 
                                             TP3_P_Value = p_values_mi[3], TP4_P_Value = p_values_mi[4]))
        
        results <- rbind(results, data.frame(Category = "MODERATE", Gene = gene, Subpopulation = sup, 
                                             P_Value = p_value_mo, TP1_P_Value = p_values_mo[1], TP2_P_Value = p_values_mo[2], 
                                             TP3_P_Value = p_values_mo[3], TP4_P_Value = p_values_mo[4]))
        # if we want the plots:
        if (anova_plot_save == "y") {
            # MINOR
            titel <- paste(gene, " expression in ",sup," monocytes (Minor) (LinReg, rE)", sep="")
            ggboxplot(df_minor, 
                      x = "Timepoint", 
                      y = "mean_rE", 
                      color = "Timepoint", 
                      title = titel,
                      add = "jitter", 
                      legend = "none", 
                      ylab = paste(gene, "expression relative to B2M"), 
                      width = 0.6, 
                      add.params = list(size = 1, alpha = 0.5)) +   
                geom_hline(yintercept = median_TP0_mi, linetype = 2) + 
                stat_compare_means(method = "anova", label.y = max(df_minor$mean_rE)) +        
                stat_compare_means(label = "p.signif", 
                                   method = "wilcox", 
                                   ref.group = "TP0", 
                                   hide.ns = TRUE, 
                                   label.y = max(df_minor$mean_rE)) +
                scale_color_manual(values = my_colors)
            # save
            file_name <- paste("ANOVA_rE_Category/",titel, ".png", sep = "")
            ggsave(filename=file_name)
            
            # MODERATE
            titel <- paste(gene, " expression in ",sup," monocytes (Moderate) (LinReg, rE)", sep="")
            
            # All together 
            ggboxplot(df_moderate, 
                      x = "Timepoint", 
                      y = "mean_rE", 
                      color = "Timepoint", 
                      title = titel,
                      add = "jitter", 
                      legend = "none", 
                      ylab = paste(gene, "expression relative to B2M"), 
                      width = 0.6, 
                      add.params = list(size = 1, alpha = 0.5)) +   
                geom_hline(yintercept = median_TP0_mo, linetype = 2) + 
                stat_compare_means(method = "anova", label.y = max(df_moderate$mean_rE)) +        
                stat_compare_means(label = "p.signif", 
                                   method = "wilcox", 
                                   ref.group = "TP0", 
                                   hide.ns = TRUE, 
                                   label.y = max(df_moderate$mean_rE)) +
                scale_color_manual(values = my_colors)
            # save
            file_name <- paste("ANOVA_rE_Category/",titel, ".png", sep = "")
            ggsave(filename=file_name)
        }  
    }
}

# Add significance levels
results$TP1_Significance <- sapply(results$TP1_P_Value, get_significance)
results$TP2_Significance <- sapply(results$TP2_P_Value, get_significance)
results$TP3_Significance <- sapply(results$TP3_P_Value, get_significance)
results$TP4_Significance <- sapply(results$TP4_P_Value, get_significance)
ANOVA_Category <- results
# Save ANOVA results with significance levels
write.csv(ANOVA_Category, file = "ANOVA_results/ANOVA_Results_Category.csv", row.names = FALSE)

# Combine all ANOVA results--------------------------------------------------------
# ANOVA_dCT$Sex <- "All"
# ANOVA_dCT$Category <- "All"
ANOVA_rE$Sex <- "All"
ANOVA_rE$Category <- "All"
ANOVA_Sex$Category <- "All"
ANOVA_Category$Sex <- "All"

ANOVA_All <- rbind(ANOVA_rE, ANOVA_Sex,ANOVA_Category)
ANOVA_Sig <- ANOVA_All %>% filter(P_Value <= 0.05)

write.csv(ANOVA_All, file = "ANOVA_results/ANOVA_Results_All.csv", row.names = FALSE)
write.csv(ANOVA_Sig, file = "ANOVA_results/ANOVA_Results_Sig.csv", row.names = FALSE)


# *Regression* --------------------------------------------------------
## Pearson Gene vs Age --------------------------------------------------------

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

## Extract significant LinRegs --------------------------------------------------------

Age_correlation_rE$Sex <- "All"
All_Age_Pearson_rE <- rbind(Age_correlation_Sex_rE, Age_correlation_rE)
Age_sigGenes_Pearson <- All_Age_Pearson_rE %>% filter(All_Age_Pearson_rE$p.value <0.05)
write.csv(Age_sigGenes_Pearson, "LinReg_Results/Age_sigGenes_Pearson.csv", row.names = FALSE)

## Plot significant LinRegs (F/M) --------------------------------------------------------
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
