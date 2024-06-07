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

## Calculate VIF for lmm_model_dCT
calculate_vif <- function(model) {
  X <- model.matrix(model)
  vif_values <- sapply(1:ncol(X), function(i) {
    x <- X[, i]
    lm_i <- lm(x ~ ., data = as.data.frame(X[, -i]))
    rsquared_i <- summary(lm_i)$r.squared
    if (is.nan(rsquared_i)) {
      return(Inf)
    } else {
      return(1 / (1 - rsquared_i))
    }
  })
  names(vif_values) <- colnames(X)
  return(vif_values)
}

### *Data* ---------------------------------------------------------------------
## Loading of the Data
setwd("/Users/ju5263ta/Github/Monocytes/Data/Fluidigm")
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

# Bring data in the right format.
names(raw_data)[names(raw_data) == 'Name'] <- 'Sample'
names(raw_data)[names(raw_data) == 'Name.1'] <- 'Gene'
names(raw_data)[names(raw_data) == 'Value.1'] <- 'dCT_Value'

# remove values from the run that had a bubble:
#raw_data <- raw_data[raw_data$Sample != "I_20_Ctr43_B1", ]

## only include runs with working B2M-----------------------------------
#extract failed runs
failed_T1 <- raw_data %>%
  filter(Gene == 'B2M_T1' & Value >= 30) %>%
  select(Sample) %>%
  distinct()
failed_T2 <- raw_data %>%
  filter(Gene == 'B2M_T2' & Value >= 30) %>%
  select(Sample) %>%
  distinct()
failed_B2M_T1 <- failed_T1$Sample
failed_B2M_T2 <- failed_T2$Sample
# Flag B2M as failed 
raw_data <- raw_data %>%
  mutate(Call.1 = if_else(Sample %in% failed_B2M_T1 & Gene == 'B2M_T1', 'mFlag', Call.1))
raw_data <- raw_data %>%
  mutate(Call.1 = if_else(Sample %in% failed_B2M_T2 & Gene == 'B2M_T2', 'mFlag', Call.1))
# Flag other genes as failed 
raw_data <- raw_data %>%
  mutate(Call.1 = if_else(Sample %in% failed_B2M_T1 & Reference == 'B2M_T1', 'mFlag', Call.1))
raw_data <- raw_data %>%
  mutate(Call.1 = if_else(Sample %in% failed_B2M_T2 & Reference == 'B2M_T2', 'mFlag', Call.1))

#mFlag_count <- sum(raw_data$Call == 'mFlag')
#print(mFlag_count)
## trimming of data ------------------------------------------------------
data <- cbind.data.frame(Sample=raw_data$Sample,
                         Gene=raw_data$Gene,
                         Value=raw_data$Value,
                         dCT_Value=raw_data$dCT_Value,
                         Call=raw_data$Call.1)

# evaluate PASS FAIL------------------------------------------------
data$Value<-ifelse(data$Call=='Flag', yes=35, no=data$Value) #999 
data$Value<-ifelse(data$Call=='mFlag', yes=NA, no=data$Value) #999 
data$dCT_Value<-ifelse(data$Call=='Flag', yes=35, no=data$dCT_Value) #999
data$dCT_Value<-ifelse(data$Call=='mFlag', yes=NA, no=data$dCT_Value) #999 instead of NA

# sort through quality and give NA to everything that Failed (=999) or is above CT 35
data$dCT_Value<-ifelse(data$Value>=35, yes=35, no=data$dCT_Value) 
data$dCT_Value<-ifelse(data$dCT_Value>=35, yes=35, no=data$dCT_Value) 

#data$Value<-ifelse(data$Value>=35, yes=35, no=data$Value) 
#remove empty rows (there should not be any...)
#data<-data[!data$Gene=="",]

# take the data
dataTRIM <-cbind.data.frame(Sample=data$Sample,
                            Gene=data$Gene,
                            CT_Value=data$Value,
                            dCT_Value=data$dCT_Value)

## Information on the data ---------------------------------------------------------------------
# adjust the indicators for the genes list here (for the heatmap) - line 138:

names(dataTRIM)[names(dataTRIM) == 'Gene'] <- 'Gene1'
dataTRIM$TechnicalReplicate <- raw_data$Gene
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
#dataTRIM <- dataTRIM %>% filter(Cell_Count != "200")
# remove all Timepoint 5
dataTRIM <- dataTRIM %>% filter(Timepoint != "TP5")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test01")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test02")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test03")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test04")

dataTRIM$Cell_Count <- NULL
#dataTRIM$Sample <- NULL
#dataTRIM$Gene1 <- NULL

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

## load all Metadata ---------------------------------------------------------------
setwd("/Users/ju5263ta/Github/Monocytes/Data/Metadata")
getwd ()
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

# *Imputation of missing data* ---------------------------------------------------------------------

## Failed Genes --------------------------------------------------------
fails <- dataAll %>%
  group_by(Subpopulation, Gene, Sex) %>%
  summarize(
    N = n(),
    Count_Above_35_CT <- sum(CT_Value >= 35, na.rm = TRUE),
    Per_Above_35_CT = mean(CT_Value >= 35, na.rm = TRUE) * 100,
    Count_Above_35_dCT = sum(dCT_Value >= 35, na.rm = TRUE),
    Per_Above_35_dCT = mean(dCT_Value >= 35, na.rm = TRUE) * 100
  )

failed_genes <- fails %>% filter(Per_Above_35_CT >= 96) %>%select(Gene)
wrong_threshold <- c("CLEC7A", "S100A8", "S100A9")
exclude_genes  <- unique(c("Xeno",failed_genes$Gene, wrong_threshold))

#Make data frame for the imputation with only working genes
dataWORKING <- dataAll %>%  filter(!Gene %in% exclude_genes)

## Prepare Matrix for imputation ---------------------------------------------------------------------
# visualize the missing values:
vis_miss(dataWORKING)
hist(dataWORKING$dCT_Value, main = "Histogram of Original Data")

initial2<-mice(dataWORKING, maxit=0, print=F)
# Check
initial2$method
# exclude form the imputation
initial2$method["SampleID"]<-"none"
initial2$method["Sample"]<-"none"
initial2$method["Gene1"]<-"none"
initial2$method["TechnicalReplicate"]<-"none"
initial2$method["BiologicalReplicate"]<-"none"
initial2$method["Category"]<-"none"
initial2$method["Age"]<-"none"
initial2$method["Sex"]<-"none"
initial2$method["Gene"]<-"none"

# We imputes missing values assuming that the data follows a normal (Gaussian) distribution.
initial2$method["CT_Value"]<-"norm"
initial2$method["dCT_Value"]<-"norm"

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
#data_imputed<-mice(dataWORKING, m=20, maxit=10, seed=1234, meth=initial2$method, pred=initial2$predictorMatrix)
#dataimputed100<-mice(data, m=20, maxit=100, seed=1234, meth=initial1$method, pred=initial1$predictorMatrix)
#data_imputed<-mice(data_forImputation, m=20, maxit=10, seed=1234, meth=initial1$method, pred=initial1$predictorMatrix)
###### data_imputed<-mice(dataWORKING, m=20, maxit=10, seed=1234, meth=initial2$method, pred=initial2$predictorMatrix)

### ALTERNATIVE
genesWORKING <- unique(dataWORKING$Gene)

# Initialize an empty list to store imputed dataframes
imputed_data <- list()

# Loop through each gene
for (gene in genesWORKING) {
  cat("Imputing for gene:", gene, "\n")
  
  # Subset data for the current gene
  gene_data <- filter(dataWORKING, Gene == gene)
  
  # Impute missing values for the current gene
  imputed_gene_data <- mice(gene_data, m = 20, maxit = 10, seed = 1234, meth = initial2$method, pred = initial2$predictorMatrix)
  
  completed_gene_data <- complete(imputed_gene_data, "long")
  completed_gene_data$Gene <- gene
  
  # Store imputed data in the list
  imputed_data[[gene]] <- completed_gene_data
}

data_imputed <- lapply(imputed_data, as.data.frame)

finaldata <- Reduce(function(x, y) merge(x, y, by = names(x), all = TRUE), imputed_data)

# QQ-plot of regression residuals and plot of residuals vs fitted values
# not needed, just for QC
# extract the imputed data into a “normal” data frame, run the analyses on each imputation separately and do model diagnostics.
#finaldata <- complete(data_imputed, "long")

#View(finaldata)
table(finaldata$.imp)

data_imputed_mean <- finaldata  %>%
  group_by(SampleID, Sample, Gene1, Gene, TechnicalReplicate, BiologicalReplicate, Timepoint, Category, Subpopulation, Age, Sex) %>%    
  summarize(CT_Mean = mean(CT_Value, na.rm = TRUE),
            CT_sd = sd(CT_Value, na.rm = TRUE), 
            CT_N = sum(!is.na(CT_Value)), 
            dCT_Mean = mean(dCT_Value, na.rm = TRUE), 
            dCT_sd = sd(dCT_Value, na.rm = TRUE),
            .groups = "drop")
dataFINAL <- as.data.frame(data_imputed_mean)

# calculated the relative expression
dataFINAL$rE_Value <- 2^((-1)*dataFINAL$dCT_Mean)
dataFINAL$rE_Value_sd <- abs(dataFINAL$rE_Value) * dataFINAL$dCT_sd

dataFINAL$rE_Value <- ifelse(dataFINAL$CT_Mean>=35, yes=0, no=dataFINAL$rE_Value)
dataFINAL$rE_Value_sd <- ifelse(dataFINAL$CT_Mean>=35, yes=0, no=dataFINAL$rE_Value) 
dataFINAL$rE_Value <- ifelse(dataFINAL$dCT_Mean>=35, yes=0, no=dataFINAL$rE_Value) 
dataFINAL$rE_Value_sd <-ifelse(dataFINAL$dCT_Mean>=35, yes=0, no=dataFINAL$rE_Value) 


## Effect of different variables ---------------------------------------------------------------------
# not needed, just for QC

## Effect of the different variables on the data or better on the dCT_Value!
# with the original  data
#lm(dCT_Mean~Gene+Age+Sex+Subpopulation+SampleID+Timepoint, data = dataFINAL) %>% sjPlot::tab_model()

# with the original  data but without individuals!
#lm(dCT_Mean~Gene+Age+Sex+Subpopulation+Timepoint, data = dataFINAL) %>% sjPlot::tab_model()
#lm(dCT_Mean~Gene+Age+Sex+Subpopulation+Timepoint, data = dataFINAL) %>% sjPlot::tab_model()

# with the imputed data
#lm(dCT_Mean~Gene+Age+Sex+Subpopulation+SampleID, data = data_imputed) %>% sjPlot::tab_model()

# similar but how to get the variables out:
# one is the actual analysis, here a linear regression with the function lm. 
#fit2<-with(data_imputed,lm(dCT_Value~Gene+Age+Sex+Subpopulation+SampleID))
#fit1<-with(data_imputed,lm(dCT_Value~Gene+Age+Sex+Subpopulation))

# The other stage is pooling of the estimates, which is done with function pool.
#res2<-pool(fit2)
#res1<-pool(fit1)
#summary(res2, conf.int=TRUE)
#summary(res1, conf.int=TRUE)

## run this instead of imputation --------------------------------------------------
#dataFINAL <- dataAll 
#dataFINAL$rE_Value <- 2^((-1)*dataFINAL$dCT_Value)
#file_name <- paste("dataFINAL_woImputation.xlsx", sep = "")
#write_xlsx(dataFINAL, file_name)



# *FACS data* ----------------------------------------------------------------------------------------------------

## Gene Metadata & other edits:
Gene_Info <- data.frame(Gene = Metadata_GeneID$Gene, GeneID = Metadata_GeneID$GeneID, Gene_Group = Metadata_GeneID$Gene_Group)

dataFINAL <- merge(dataFINAL, Gene_Info, by = "Gene")

# List of Patients
patients <- unique(dataFINAL$SampleID)
no_patients <- length(patients)

## import FACS data ----------------------------------------------------------------------------------------------------
setwd("/Users/ju5263ta/Github/Monocytes/Data/FACS")
getwd ()

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
#write_xlsx(FACS_info_rows, "FACS_info_rows.xlsx")
FACSdata <- cbind(FACSdata, FACS_info_rows)

# Match PlateID with the right sample with the right data
FACSdata <- merge(FACSdata, Plate_info, by = "PlateID")
# remove the 5th timepoint
FACSdata <- FACSdata %>% filter(Timepoint != "TP5")

## Days Post-Stroke data
# DaysPS <- read_excel("Metadata_DaysPS.xls")
dataFINAL$DaysPS <- as.numeric(NA)
dataFINAL$DaysPS[grep("TP1",dataFINAL$Timepoint)]  <- 1
dataFINAL$DaysPS[grep("TP2",dataFINAL$Timepoint)]  <- 4
dataFINAL$DaysPS[grep("TP3",dataFINAL$Timepoint)]  <- 30
dataFINAL$DaysPS[grep("TP4",dataFINAL$Timepoint)]  <- 90
dataFINAL$DaysPS[grep("TP5",dataFINAL$Timepoint)]  <- 360
#unique(dataFINAL$DaysPS)

## create mean values for the ANOVA & so
dataFINALmean <- dataFINAL %>%
        group_by(SampleID, Sex, Age, Category, Timepoint, DaysPS,
                 Subpopulation, Gene, Gene_Group, GeneID) %>%
        summarise(mean_CT = mean(CT_Mean), sd_CT = sd(CT_Mean),
                  mean_dCT = mean(dCT_Mean), sd_dCT = sd(dCT_Mean),
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


# *ANOVA* fix the p value ----------------------------------------------------------------------------------------------------

## dCT with Baseline - with Outlier exclusion --------------------------------------------------
createFolder("ANOVA_results")
folder <- "ANOVAwilcox_dCT exOutliers"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric())

for (j in 1:length(unique(dataFINALmean$Gene))) {
    gene <- dataFINALmean$Gene[j]
    df_gene <- filter(dataFINALmean, Gene == gene)
    df_gene <- df_gene %>% filter(!is.na(Category), !is.na(mean_dCT), !is.na(Subpopulation),  !is.na(Timepoint))
    
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
        #posthoc_result <- TukeyHSD(anova_result)
         
        # Extract p-value from ANOVA summary
        p_value <- summary(anova_result)[[1]]$"Pr(>F)"[1]
        
        # Save results to data frame
        results <- rbind(results, data.frame(Gene = gene, Subpopulation = sup, P_Value = p_value))
        #my_colors <- c("black", "red", "orange", "#b8860b", "green") 
        my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 

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

ANOVA_dCT <- results
write.csv(ANOVA_dCT, file  = "ANOVA_results//ANOVA_Results_dCT.csv", row.names = FALSE)

## rE with Baseline - add adjusted p-Value ? --------------------------------------------------

folder <- "ANOVAwilcox_rE exOutliers"
results <- data.frame(Gene = character(), Subpopulation = character(), P_Value = numeric())

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
        
        # Extract p-value from ANOVA summary
        p_value <- summary(anova_result)[[1]]$"Pr(>F)"[1]
        # Extract adjusted p-value (NA for ANOVA)
        adjusted_p_value <- NA
        
        # Save results to data frame
        results <- rbind(results, data.frame(Gene = gene, Subpopulation = sup, P_Value = p_value))
        #my_colors <- c("black", "red", "orange", "#b8860b", "green") 
        my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 

# Plot the boxplot and comparison without outliers
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
ANOVA_rE <- results
write.csv(ANOVA_rE, file = "ANOVA_results/ANOVA_Results_rE.csv", row.names = FALSE)

## Male vs Female --------------------------------------------------
folder <- "ANOVA_rE_Sex"
results <- data.frame(Sex = character(), Gene = character(), Subpopulation = character(), P_Value = numeric())

i <- 1
for (i in 1:length(unique(dataFINALmean$Gene))) {
    gene <- dataFINALmean$Gene[i]
    df_gene <- filter(dataFINALmean, Gene == gene)
  df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_rE), !is.na(Subpopulation), !is.na(Sex))
  j <- 1
  for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- df_gene$Subpopulation[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)

    df_subpop$Timepoint <- factor(as.vector(df_subpop$Timepoint))
    df_all <- df_subpop[!df_subpop$mean_rE  %in% boxplot.stats(df_subpop$mean_rE )$out, ]
    
    if (nrow(df_subpop) <= 0) {
      next
    }
    tryCatch({
    titel <- paste(gene, " expression in ",sup," monocytes (LinReg, rE)", sep="")

# Separeted by Facets        
ggboxplot(df_all, 
          x = "Sex", 
          y = "mean_rE", 
          color = "Sex",  
          title = titel,
          add = "jitter", 
          facet.by = "Timepoint",
          legend = "none", 
          ylab = paste(gene, "expression relative to B2M"), 
          width = 0.8, 
          add.params = list(size = 1, alpha = 0.5)) +  
  stat_compare_means(method = "anova", label.y = max(df_all$mean_rE)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "Male", 
                     hide.ns = TRUE, 
                     label.y = max(df_all$mean_rE)) +
  scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2"))

    file_name <- paste("ANOVA_rE_Sex/",titel, "_Facets.png", sep = "")
    ggsave(filename=file_name)

# Male plot
male <- df_subpop %>% filter(Sex == "Male")        
titel <- paste(gene, " expression in ",sup," (Male) monocytes (LinReg, rE)", sep="")
male <- male[!male$mean_rE  %in% boxplot.stats(male$mean_rE )$out, ]

baseline <- filter(male, Timepoint == "TP0")
median_TP0 <- median(baseline$mean_rE)

# ANOVA
anova_m <- aov(mean_rE ~ Timepoint, data = male)
p_value_m <- summary(anova_m)[[1]]$"Pr(>F)"[1]

results <- rbind(results, data.frame(Sex = male$Sex[1], Gene = gene, Subpopulation = sup, P_Value = p_value_m))

ggboxplot(male, 
          x = "Timepoint", 
          y = "mean_rE", 
          color = "Timepoint", 
          title = titel,
          add = "jitter", 
          legend = "none", 
          ylab = paste(gene, "expression relative to B2M"), 
          width = 0.6, 
          add.params = list(size = 1, alpha = 0.5)) +  
     geom_hline(yintercept = median_TP0, linetype = 2) + 
  stat_compare_means(method = "anova", label.y = max(male$mean_rE)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "TP0", 
                     hide.ns = TRUE, 
                     label.y = max(male$mean_rE)) +
  scale_color_manual(values = my_colors)

    file_name <- paste("ANOVA_rE_Sex/",titel, ".png", sep = "")
    ggsave(filename=file_name)
    
# Female plot
female <- df_subpop %>% filter(Sex == "Female")        
titel <- paste(gene, " expression in ",sup," (Female) monocytes (LinReg, rE)", sep="")
female <- female[!female$mean_rE  %in% boxplot.stats(female$mean_rE )$out, ]

baseline <- filter(female, Timepoint == "TP0")
median_TP0 <- median(baseline$mean_rE)

# ANOVA
anova_f <- aov(mean_rE ~ Timepoint, data = female)
p_value_f <- summary(anova_f)[[1]]$"Pr(>F)"[1]

results <- rbind(results, data.frame(Sex = female$Sex[1], Gene = gene, Subpopulation = sup, P_Value = p_value_f))

ggboxplot(female, 
          x = "Timepoint", 
          y = "mean_rE", 
          color = "Timepoint", 
          title = titel,
          add = "jitter", 
          legend = "none", 
          ylab = paste(gene, "expression relative to B2M"), 
          width = 0.6, 
          add.params = list(size = 1, alpha = 0.5)) +  
     geom_hline(yintercept = median_TP0, linetype = 2) + 
  stat_compare_means(method = "anova", label.y = max(female$mean_rE)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "TP0", 
                     hide.ns = TRUE, 
                     label.y = max(female$mean_rE)) +
  scale_color_manual(values = my_colors)

    file_name <- paste("ANOVA_rE_Sex/",titel, ".png", sep = "")
    ggsave(filename=file_name)
    }
  )
}
}

ANOVA_Sex <- results
write.csv(ANOVA_Sex, file = "ANOVA_results/ANOVA_Results_Sex.csv", row.names = FALSE)

## Healthy vs MILD vs MODERATE --------------------------------------------------------

folder <- "ANOVA_rE_Category"
results <- data.frame(Category = character(), Gene = character(), Subpopulation = character(), P_Value = numeric())

i <- 1
for (i in 1:length(unique(dataFINALmean$Gene))) {
    gene <- dataFINALmean$Gene[i]
    df_gene <- filter(dataFINALmean, Gene == gene)
   df_gene <- df_gene %>% filter(!is.na(Timepoint), !is.na(mean_rE), !is.na(Subpopulation), !is.na(Category))
  j <- 1
  for (j in 1:length(unique(df_gene$Subpopulation))) {
        sup <- df_gene$Subpopulation[j]
        df_subpop <- filter(df_gene, Subpopulation == sup)
    
    df_minor <- df_subpop %>% filter(Category != "MODERATE")
    df_moderate <- df_subpop %>% filter(Category != "MINOR")
    
    df_minor <- df_minor[!df_minor$mean_rE %in% boxplot.stats(df_minor$mean_rE)$out, ]
    df_moderate <- df_moderate[!df_moderate$mean_rE %in% boxplot.stats(df_moderate$mean_rE)$out, ]
    
    df_minor$Timepoint <- factor(as.vector(df_minor$Timepoint))
    df_moderate$Timepoint <- factor(as.vector(df_moderate$Timepoint))
    
    baseline_mi <- filter(df_minor, Timepoint == "TP0")
    baseline_mo <- filter(df_moderate, Timepoint == "TP0")
    median_TP0_mi <- median(baseline_mi$mean_rE)
    median_TP0_mo <- median(baseline_mo$mean_rE)
    
    # ANOVA
    anova_mi <- aov(mean_rE ~ Timepoint, data = df_minor)
    anova_mo <- aov(mean_rE ~ Timepoint, data = df_moderate)
    p_value_mi <- summary(anova_mi)[[1]]$"Pr(>F)"[1]
    p_value_mo <- summary(anova_mo)[[1]]$"Pr(>F)"[1]

    results <- rbind(results, data.frame(Category = "MINOR", Gene = gene, Subpopulation = sup, P_Value = p_value_mi))
    results <- rbind(results, data.frame(Category = "MODERATE", Gene = gene, Subpopulation = sup, P_Value = p_value_mo))
    
    if (nrow(df_minor) <= 0) {
      next
    }
    tryCatch({
    titel <- paste(gene, " expression in ",sup," monocytes (Minor) (LinReg, rE)", sep="")

# All together 
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

    file_name <- paste("ANOVA_rE_Category/",titel, ".png", sep = "")
    ggsave(filename=file_name)
    })

    if (nrow(df_minor) <= 0) {
      next
    }
    tryCatch({
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

    file_name <- paste("ANOVA_rE_Category/",titel, ".png", sep = "")
    ggsave(filename=file_name)
    }
  )
  }
}

ANOVA_Category <- results
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

# *Housekeeping gene evaluation* --------------------------------------------------------
HK_A <- filter(dataFINAL, Gene=="ACTB")
HK_B <- filter(dataFINAL, Gene=="B2M")
HK_G <- filter(dataFINAL, Gene=="GAPDH")
HK_A$Timepoint <- factor(as.vector(HK_A$Timepoint))
HK_B$Timepoint <- factor(as.vector(HK_B$Timepoint))
HK_G$Timepoint <- factor(as.vector(HK_G$Timepoint))
        
anova_A <- aov(CT_Mean ~ Timepoint, data = HK_A)
anova_B <- aov(CT_Mean ~ Timepoint, data = HK_B)
anova_G <- aov(CT_Mean ~ Timepoint, data = HK_G)
         
p_value <- data.frame(HKgene = c("ACTB", "B2M", "GAPDH"), P_Value = c(summary(anova_A)[[1]]$"Pr(>F)"[1],                                 summary(anova_B)[[1]]$"Pr(>F)"[1],                               summary(anova_G)[[1]]$"Pr(>F)"[1]))
        
my_colors <- c("darkgreen", "orange", "red", "magenta", "purple") 

# ACTB ----------------------------------------------------------------
ggboxplot(HK_A, 
          x = "Timepoint", 
          y = "CT_Mean", 
          color = "Timepoint", 
          add = "jitter", 
          legend = "none", 
          ylab = "ACTB expression [CT]", 
          width = 0.8, 
          add.params = list(size = 1, alpha = 0.5)) +  
  stat_compare_means(method = "anova", label.y = max(HK_A$CT_Mean)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "TP0", 
                     hide.ns = TRUE, 
                     label.y = max(HK_A$CT_Mean)) +
  scale_color_manual(values = my_colors)
ggsave(filename="HK_evaluation/ACTB.png")

HK_A_mean <- HK_A %>%
        group_by(SampleID, Timepoint, Subpopulation) %>%
        summarise(mean_CT = mean(CT_Mean), sd_CT = sd(CT_Mean))

ggboxplot(HK_A_mean, 
          x = "Timepoint", 
          y = "mean_CT", 
          color = "Timepoint", 
          add = "jitter", 
          legend = "none", 
          ylab = "ACTB expression [CT]",  
          width = 0.8, 
          add.params = list(size = 1, alpha = 0.5)) +  
  stat_compare_means(method = "anova", label.y = max(HK_A_mean$mean_CT)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "TP0", 
                     hide.ns = TRUE, 
                     label.y = max(HK_A_mean$mean_CT)) +
  scale_color_manual(values = my_colors)
ggsave(filename="HK_evaluation/ACTB_mean.png")

# B2M -----------------------------------------------------------------
ggboxplot(HK_B, 
          x = "Timepoint", 
          y = "CT_Mean", 
          color = "Timepoint", 
          add = "jitter", 
          legend = "none", 
          ylab = "B2M expression [CT]", 
          width = 0.8, 
          add.params = list(size = 1, alpha = 0.5)) +  
  stat_compare_means(method = "anova", label.y = max(HK_B$CT_Mean)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "TP0", 
                     hide.ns = TRUE, 
                     label.y = max(HK_B$CT_Mean)) +
  scale_color_manual(values = my_colors)
ggsave(filename="HK_evaluation/B2M.png")

HK_B_mean <- HK_B %>%
        group_by(SampleID, Timepoint, Subpopulation) %>%
        summarise(mean_CT = mean(CT_Mean), sd_CT = sd(CT_Mean))

ggboxplot(HK_B_mean, 
          x = "Timepoint", 
          y = "mean_CT", 
          color = "Timepoint", 
          add = "jitter", 
          legend = "none", 
          ylab = "B2M expression [CT]", 
          width = 0.8, 
          add.params = list(size = 1, alpha = 0.5)) +  
  stat_compare_means(method = "anova", label.y = max(HK_B_mean$mean_CT)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "TP0", 
                     hide.ns = TRUE, 
                     label.y = max(HK_B_mean$mean_CT)) +
  scale_color_manual(values = my_colors)
ggsave(filename="HK_evaluation/B2M_mean.png")

# GAPDH ----------------------------------------------------------------
ggboxplot(HK_G, 
          x = "Timepoint", 
          y = "CT_Mean", 
          color = "Timepoint", 
          add = "jitter", 
          legend = "none", 
          ylab = "GAPDH expression [CT]", 
          width = 0.8, 
          add.params = list(size = 1, alpha = 0.5)) +  
  stat_compare_means(method = "anova", label.y = max(HK_G$CT_Mean)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "TP0", 
                     hide.ns = TRUE, 
                     label.y = max(HK_G$CT_Mean)) +
  scale_color_manual(values = my_colors)
ggsave(filename="HK_evaluation/GAPDH.png")

HK_G_mean <- HK_G %>%
        group_by(SampleID, Timepoint, Subpopulation) %>%
        summarise(mean_CT = mean(CT_Mean), sd_CT = sd(CT_Mean))

ggboxplot(HK_G_mean, 
          x = "Timepoint", 
          y = "mean_CT", 
          color = "Timepoint", 
          add = "jitter", 
          legend = "none", 
          ylab = "GAPDH expression [CT]", 
          width = 0.8, 
          add.params = list(size = 1, alpha = 0.5)) +  
  stat_compare_means(method = "anova", label.y = max(HK_G_mean$mean_CT)) +        
  stat_compare_means(label = "p.signif", 
                     method = "wilcox", 
                     ref.group = "TP0", 
                     hide.ns = TRUE, 
                     label.y = max(HK_G_mean$mean_CT)) +
  scale_color_manual(values = my_colors)
ggsave(filename="HK_evaluation/GAPDH_mean.png")

anova_A_mean <- aov(mean_CT ~ Timepoint, data = HK_A_mean)
anova_B_mean <- aov(mean_CT ~ Timepoint, data = HK_B_mean)
anova_G_mean <- aov(mean_CT ~ Timepoint, data = HK_G_mean)
         
p_value$P_Value_mean <- c(summary(anova_A_mean)[[1]]$"Pr(>F)"[1],                                 summary(anova_B_mean)[[1]]$"Pr(>F)"[1],                               summary(anova_G_mean)[[1]]$"Pr(>F)"[1])

write.csv(p_value, file = "HK_evaluation/HK_pValues.csv", row.names = FALSE)

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
