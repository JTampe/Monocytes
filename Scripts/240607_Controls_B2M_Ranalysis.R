
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
#library(nondetects)

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("mice")
#BiocManager::install("skimr")
#BiocManager::install("visdat")
# Load the mice package
library(mice)
library(skimr)
library(visdat) 

thematic_on(bg = "white", fg = "black", accent = "blue")

### *Functions* ---------------------------------------------------------------------

## Create folder
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

## LinReg Plot 
generateLinRegPlot <- function(data, xvar, yvar, title, folder, filename) {
  file_name <- file.path(folder, paste(filename, ".png", sep = ""))
  
  ggplot(data = data, aes_string(x = xvar, y = yvar)) +
    geom_smooth(method = "glm", color = "black") +
    geom_point(aes(color = Sex), size = 2) +
    ggtitle(title) +
    xlab(xvar) +
    ylab(yvar) +
    scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2")) +
    theme(text = element_text(size = 14)) +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = "black"),
      text = element_text(color = "black"),
      panel.grid.major = element_line(color = "gray", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  ggsave(filename = file_name)
}

## LinReg Plot - sex
generateLinRegPlot_Sex <- function(data, xvar, yvar, title, folder, filename) {
  file_name <- file.path(folder, paste(filename, ".png", sep = ""))
  
  p <- ggplot(data = data, aes_string(x = xvar, y = yvar)) +
    geom_smooth(method = "glm", aes(color = Sex), size = 1) +
    geom_point(aes(color = Sex), size = 2) +
    ggtitle(title) +
    xlab(xvar) +
    ylab(yvar) +
    scale_color_manual(values = c("Male" = "#A67C00", "Female" = "#1D04C2")) +
    theme(text = element_text(size = 14)) +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = "black"),
      text = element_text(color = "black"),
      panel.grid.major = element_line(color = "gray", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  ggsave(filename = file_name, plot = p)
}

## Weighted Correlation
weighted_correlation <- function(x, y, sd_x, sd_y) {
  n <- length(x)
  if (n != length(y) || n != length(sd_x) || n != length(sd_y)) {
    stop("Input vectors must have the same length.")
  }
  
  cov_xy <- sum((x - mean(x)) * (y - mean(y)) * sd_x * sd_y) / sum(sd_x^2 * sd_y^2)
  var_x <- sum((x - mean(x))^2 * sd_x^2) / sum(sd_x^2)
  var_y <- sum((y - mean(y))^2 * sd_y^2) / sum(sd_y^2)
  
  weighted_cor <- cov_xy / sqrt(var_x * var_y)
  
  return(weighted_cor)
}

# *Data*
## Loading of the Data -----------------------------------------
setwd("/Users/ju5263ta/Github/Monocytes/Data/Fluidigm")
getwd ()
# make sure that before you check the files, that not an empty column has been added to the end of the file, remove and export as .csv otherwise.
files <- system( "ls *.csv", intern=T)

raw_data <- read_csv("Controls_B2M.csv", 
    skip = 11)
# delete Call.1 column, and resent column names
if ( length(colnames(raw_data)) > 13) {
    raw_data <- raw_data[, !colnames(raw_data) %in% "Call.1"]
    colnames(raw_data) <- c("ID","Name","Type","rConc","Name.1","Type.1","Reference","Value", "Quality","Call","Threshold","Value.1","Quality.1","Call.1")
}
i <- 1
raw_data$Dataset<-files[i]

# give the right labels and combine in one big file
names(raw_data)[names(raw_data) == 'Name'] <- 'Sample'
names(raw_data)[names(raw_data) == 'Name.1'] <- 'Gene'
names(raw_data)[names(raw_data) == 'Value.1'] <- 'dCT_Value'

# remove values from the run that had a bubble:
raw_data <- raw_data[raw_data$Sample != "I_20_Ctr43_B1", ]

# trimming of data -------------------------------------------------------
data <- cbind.data.frame(Sample=raw_data$Sample,
                         Gene=raw_data$Gene,
                         Value=raw_data$Value,
                         dCT_Value=raw_data$dCT_Value,
                         Call=raw_data$Call)

#evaluate PASS FAIL
data$Value<-ifelse(data$Call=='Flag', yes=999, no=data$Value) #999 instead of NA
data$Value<-ifelse(data$Call=='mFlag', yes=NA, no=data$Value) #999 instead of NA

data$dCT_Value<-ifelse(data$Call=='Flag', yes=35, no=data$dCT_Value) #999 instead of NA
data$dCT_Value<-ifelse(data$Call=='mFlag', yes=NA, no=data$dCT_Value) #999 instead of NA

# sort through quality and give and NA to everything that Failed (=999) or is above CT 35
data$dCT_Value<-ifelse(data$Value>=35, yes=35, no=data$dCT_Value) 

#data$Value<-ifelse(data$Value>=35, yes=35, no=data$Value) 
#remove empty rows (there should not be any...)
#data<-data[!data$Gene=="",]

# take the data
dataTRIM <-cbind.data.frame(Sample=data$Sample,
                            Gene=data$Gene,
                            CT_Value=data$Value,
                            dCT_Value=data$dCT_Value)

## Information on the data
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

#Here we also split technical and (pseudo-)biological replicates....:
#without the B1/B2 segregation, we would have each read of the gene normalized against the mean of the two biological repeats of ACTB. This is in theory correct but we lose the information on the mean/SD while doing so through the machine.
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

# Remove Extra sorts
dataTRIM <- dataTRIM %>% filter(BiologicalReplicate != "Extra")
dataTRIM <- dataTRIM %>% filter(Subpopulation != "noRT")
dataTRIM <- dataTRIM %>% filter(Cell_Count != "200")

dataTRIM <- dataTRIM %>% filter(SampleID != "Test01")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test02")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test03")
dataTRIM <- dataTRIM %>% filter(SampleID != "Test04")

dataTRIM$Cell_Count <- NULL
#dataTRIM$Sample <- NULL
#dataTRIM$Gene1 <- NULL

## load all Metadata ---------------------------------------------------------------
setwd("/Users/ju5263ta/Github/Monocytes/Data/Metadata")
getwd ()
metadata <- read_excel("Metadata_Controls.xlsx")
Metadata_GeneID <- read_excel("Metadata_GeneID.xlsx")

# Combine the data frames based on common values in the "id" column
# --> here i loose all the test because they are not in the metadata
dataAll <- merge(dataTRIM, metadata, by = "SampleID", all.x = TRUE)

#relable Male & Female
dataAll$Sex[grep("man",dataAll$Sex)]  <- "Male"
dataAll$Sex[grep("kvinna",dataAll$Sex)]  <- "Female"

## remove Outliers
dataAll <- dataAll[dataAll$SampleID != "Ctr33", ]
dataAll <- dataAll[dataAll$SampleID != "Ctr15", ]
dataAll <- dataAll[dataAll$SampleID != "Ctr23", ]
dataAll <- dataAll[dataAll$SampleID != "Ctr39", ]

## remove unnessessary info (for now)

# remove test better!
dataAll <- dataAll %>% filter(!is.na(Age))
# remove Xeno
dataAll <- dataAll %>% filter(!is.na(Age))

dataAll$Cell_Count <- NULL
dataAll$Patient_ID <- NULL
dataAll$AgeID <- NULL
dataAll$SexID <- NULL
dataAll$Gene_Group <- NULL
#dataAll$Sample <- NULL
#dataAll$Gene1 <- NULL

# *Imputation of missing data* --------------------------------------------

## Failed Genes - List the failed genes as Quality Control
fails <- dataAll %>%
  group_by(Subpopulation, Gene, Sex) %>%
  summarize(
    N = n(),
    Count_Above_35_CT = sum(CT_Value > 35, na.rm = TRUE),
    Per_Above_35_CT = mean(CT_Value > 35, na.rm = TRUE) * 100,
    Count_Above_35_dCT = sum(dCT_Value > 35, na.rm = TRUE),
    Per_Above_35_dCT = mean(dCT_Value > 35, na.rm = TRUE) * 100
  )

failed_genes <- fails %>% filter(Per_Above_35_CT >= 96) %>%select(Gene) 

# exclude failed genes and the positive control Xeno
exclude_genes  <- unique(c("Xeno",failed_genes$Gene))

## Impute data ---------------------------------

#Make data frame for the imputation with only working genes
dataWORKING <- dataAll %>%  filter(!Gene %in% exclude_genes)

# visualize the missing values:
vis_miss(dataWORKING)
# QC 75% missing 1.8%

initial2<-mice(dataWORKING, maxit=0, print=F)
# Check
initial2$method
# exclude form the imputation
initial2$method["Sample"]<-"none"
initial2$method["SampleID"]<-"none"
initial2$method["Gene1"]<-"none"
initial2$method["TechnicalReplicate"]<-"none"
initial2$method["BiologicalReplicate"]<-"none"

# We imputes missing values assuming that the data follows a normal (Gaussian) distribution.
initial2$method["CT_Value"]<-"norm"
initial2$method["dCT_Value"]<-"norm"

#Check again
initial2$method

# convert to factors
dataAll$Sex <- as.factor(dataAll$Sex)
dataAll$Subpopulation <- as.factor(dataAll$Subpopulation)
dataAll$Gene <- as.factor(dataAll$Gene)

# Matrix of predictors:
#initial1$predictorMatrix
initial2$predictorMatrix

# Imputation --------------------------------------------------------------
#dataimputed100<-mice(data, m=20, maxit=100, seed=1234, meth=initial1$method, pred=initial1$predictorMatrix)
#data_imputed<-mice(data_forImputation, m=20, maxit=10, seed=1234, meth=initial1$method, pred=initial1$predictorMatrix)
data_imputed<-mice(dataWORKING, m=20, maxit=10, seed=1234, meth=initial2$method, pred=initial2$predictorMatrix)



# QQ-plot of regression residuals and plot of residuals vs fitted values
# extract the imputed data into a “normal” data frame, run the analyses on each imputation separately and do model diagnostics.
finaldata <- complete(data_imputed, "long")
#View(finaldata)
table(finaldata$.imp)

data_imputed_mean <- finaldata  %>%
  group_by(SampleID, Sample, Gene1, Gene, TechnicalReplicate, BiologicalReplicate, Subpopulation, Age, Sex) %>%    
  summarize(CT_Mean = mean(CT_Value, na.rm = TRUE),
            CT_sd = sd(CT_Value, na.rm = TRUE), 
            CT_N = sum(!is.na(CT_Value)), 
            dCT_Mean = mean(dCT_Value, na.rm = TRUE), 
            dCT_sd = sd(dCT_Value, na.rm = TRUE),
            .groups = "drop")
dataFINAL <- as.data.frame(data_imputed_mean)

## Effect of different variables
## Effect of the different variables on the data or better on the dCT_Mean!

# with the original  data
lm(dCT_Mean~Gene+Age+Sex+Subpopulation+SampleID, data = dataFINAL) %>% sjPlot::tab_model()
# with the original  data but without individuals!
lm(dCT_Mean~Gene+Age+Sex+Subpopulation, data = dataFINAL) %>% sjPlot::tab_model()

# with the imputed data
lm(dCT_Mean~Gene+Age+Sex+Subpopulation+SampleID, data = data_imputed) %>% sjPlot::tab_model()

# similar but how to get the variables out:
# one is the actual analysis, here a linear regression with the function lm. 
fit2<-with(data_imputed,lm(dCT_Value~Gene+Age+Sex+Subpopulation+SampleID))
fit1<-with(data_imputed,lm(dCT_Value~Gene+Age+Sex+Subpopulation))
# The other stage is pooling of the estimates, which is done with function pool.
res2<-pool(fit2)
res1<-pool(fit1)
summary(res2, conf.int=TRUE)
summary(res1, conf.int=TRUE)

## Check imputation & Questions

# Plots for imputation--------------------------------------------------
# Check convergence
plot(data_imputed)
# kernel density plots for imputed values 
#    densityplot(data_imputed)

### which one to look at, shows the spread of the imputated data
densityplot(data_imputed, data=~dCT_Value | Subpopulation)
densityplot(data_imputed, data=~dCT_Value | Sex)

#choose imputation number 1
# better: choose the mean from 20 impoutations
check1<-lm(dCT_Mean~Gene+Age+Sex+Subpopulation+SampleID, data=finaldata, subset=.imp==1)
#plot(check1)
dataFINAL_imp <- filter(finaldata, .imp==1)

## different approach but oly works on the rE!
dataFINAL$rE_Value <- 2^((-1)*dataFINAL$dCT_Mean)
# with the original  data but without individuals!
glm(rE_Value~Gene+Age+Sex+Subpopulation, family = "gaussian", data = dataFINAL) %>% tab_model()

# Questions:
# how is the summary interpretated? estimate? SD error? p value?
# especially in this one! what does it mean for the one thats not shown of the levels in factors?
#lm(dCT_Mean~Gene+Age+Sex+Subpopulation, data = dataFINAL) %>% sjPlot::tab_model()
# Gene vs Gene1
# exclude failed genes / neg. controls?


### What kind of plot is this? & Also not working!
#   stripplot(data_imputed)

# *FACS & Gene Metadata

## Gene Metadata & other edits:
Gene_Info <- data.frame(Gene = Metadata_GeneID$Gene, GeneID = Metadata_GeneID$GeneID, Gene_Group = Metadata_GeneID$Gene_Group)

dataFINAL <- merge(dataFINAL, Gene_Info, by = "Gene")

# Assign the right Age Groups: --> 2 Groups
i = 1
j = length(dataFINAL$Age)
for (i in (1:j)) {
  if (!is.na(dataFINAL$Age[i]) & dataFINAL$Age[i] > 65) {
    dataFINAL$Age_Group[i] <- "Older"
  } else {
    dataFINAL$Age_Group[i] <- "Adult"
  }
}

dataFINAL$Group <- NA
dataFINAL$Group <- ifelse(grepl("Female",dataFINAL$Sex) & grepl("Older",dataFINAL$Age_Group), "Female_Older", dataFINAL$Group)
dataFINAL$Group <- ifelse(grepl("Female",dataFINAL$Sex) & grepl("Adult",dataFINAL$Age_Group), "Female_Adult", dataFINAL$Group)
dataFINAL$Group <- ifelse(grepl("Male",dataFINAL$Sex) & grepl("Older",dataFINAL$Age_Group), "Male_Older", dataFINAL$Group)
dataFINAL$Group <- ifelse(grepl("Male",dataFINAL$Sex) & grepl("Adult",dataFINAL$Age_Group), "Male_Adult", dataFINAL$Group)

# List of Patients
patients <- unique(dataFINAL$SampleID)
patients

no_patients <- length(patients)
no_patients

# #Assign the right Age Groups: --> 3 Groups
# i = 1
# j = length(dataFINAL$Age)
# for (i in (1:j))  {
#     if (dataFINAL$Age[i] > 65) {
#         dataFINAL$Age_Group[i] <- "Older"
#     } else if (dataFINAL$Age[i] < 35) {
#         dataFINAL$Age_Group[i] <- "Adult"
#     } else {
#         dataFINAL$Age_Group[i] <- "Middle Aged"
#     }
# }     

## import FACS data ----------------------------------------------------------------------------------------------------
setwd("/Users/ju5263ta/Github/Monocytes/Data/FACS")
getwd ()

FACSdata <- read_excel("231002_FACSdata.xls")
FACSdata <- as.data.frame(FACSdata, col_types = c("text", 
    "text", "numeric", "text", "numeric", 
    "text", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric"))

# *Calculations*. -------------------------
## Mean, SD & N calculation
dataFINALmean <- dataFINAL %>%
  group_by(SampleID, Gene, GeneID, Gene_Group, Subpopulation, Age, Age_Group, Sex, Group) %>%    
  summarize(
    CT_Mean = mean(CT_Mean, na.rm = TRUE),
    CT_sd = sqrt(sum((CT_sd)^2, na.rm = TRUE) / sum(!is.na(CT_sd))), 
    CT_N = sum(!is.na(CT_Mean)), 
    dCT_Mean = mean(dCT_Mean, na.rm = TRUE), 
    dCT_sd = sqrt(sum((dCT_sd)^2, na.rm = TRUE) / sum(!is.na(dCT_sd))),
    .groups = "drop")

dataFINALmean <- as.data.frame(dataFINALmean)
dataFINALmean$CT_Mean[is.nan(dataFINALmean$CT_Mean)] <- NA
dataFINALmean$dCT_Mean[is.nan(dataFINALmean$dCT_Mean)] <- NA

## relative Expression ---------------------------
dataFINAL$rE_Value <- 2^((-1)*dataFINAL$dCT_Mean)
dataFINAL$rE_Value_sd <- abs(dataFINAL$rE_Value) * dataFINAL$dCT_sd
dataFINALmean$rE_Value <- 2^((-1)*dataFINALmean$dCT_Mean)
dataFINALmean$rE_Value_sd <- abs(dataFINALmean$rE_Value) * dataFINALmean$dCT_sd

# remove all lines that lack the needed information
# all lines should have everything - nothing should be removed in an optimal scenario!
dataFINALmean <- dataFINALmean %>% filter(!is.na(Age_Group), !is.na(rE_Value),!is.na(Sex), !is.na(Gene), !is.na(Sex),!is.na(Subpopulation))

## Change working directory so that everything is getting saved in the right place ----------------
today <- Sys.Date()
output_location <- paste(today,"_Stroke_Results", sep="")
setwd(paste("/Users/ju5263ta/Github/Monocytes/Data/",output_location, "/",sep=""))
getwd()

write_xlsx(fails, "Fail_rates.xlsx")
write_xlsx(dataFINAL, "dataFINAL.xlsx")
write_xlsx(dataFINALmean, "dataFINALmean.xlsx")

# *Summary files*
##  Age Summary -------------------------
# making a different dataframe for each gene:
new_folder_name <- "Summary Matrixes/"
if (!dir.exists(new_folder_name)) {
  dir.create(new_folder_name, recursive = TRUE)
}

Summary_rE_Age <- dataFINALmean  %>%
  group_by(Subpopulation, Age_Group, GeneID, Gene) %>%
  summarize(mean_value = mean(rE_Value, na.rm = TRUE), sd = sd(rE_Value, na.rm = TRUE), N = sum(!is.na(rE_Value)))
write.csv(Summary_rE_Age, paste(new_folder_name,"Summary_rE_Age", sep=""), row.names = FALSE)

## Age Summary Matrix -----------------------------------------
Summary_rE_Age$Info <- paste(Summary_rE_Age$Subpopulation, Summary_rE_Age$Age_Group, sep="_")
Summary_rE_Age$Gene_Order <- paste(Summary_rE_Age$GeneID, Summary_rE_Age$Gene, sep="_")
#Summary_rE_Age <- Summary_rE_Age[order(Summary_rE_Age$GeneID), ]
#genes_ordered <- c(unique(Summary_rE_Age$Gene))
Summary_rE_Age_Matrix <- dcast(Summary_rE_Age, Gene_Order ~ Info, value.var = "mean_value")
#Summary_rE_Age_Matrix <- Summary_rE_Age_Matrix[genes_ordered, ]

write.table(Summary_rE_Age_Matrix, paste(new_folder_name,"Summary_rE_Age_Matrix.txt", sep=""), sep="\t", row.names=FALSE)

## Sex Summary-----------------------------------------
Summary_rE_Sex <- dataFINALmean  %>%
  group_by(Subpopulation, Age_Group, Sex, GeneID, Gene) %>%
  summarize(mean_value = mean(rE_Value, na.rm = TRUE), sd = sd(rE_Value, na.rm = TRUE), N = sum(!is.na(rE_Value)))
write.csv(Summary_rE_Sex, paste(new_folder_name,"Summary_rE_Sex", sep=""), row.names = FALSE)

## Sex Summary Matrix-----------------------------------------
Summary_rE_Sex$Info <- paste(Summary_rE_Sex$Subpopulation, Summary_rE_Sex$Sex, Summary_rE_Sex$Age_Group, sep="_")
Summary_rE_Sex$Info2 <- paste(Summary_rE_Sex$Sex, Summary_rE_Sex$Subpopulation, Summary_rE_Sex$Age_Group, sep="_")
Summary_rE_Sex$Gene_Order <- paste(Summary_rE_Sex$GeneID, Summary_rE_Sex$Gene, sep="_")
#Summary_rE_Sex <- Summary_rE_Sex[order(Summary_rE_Sex$GeneID), ]
#genes_ordered <- c(unique(Summary_rE_Sex$Gene))
Summary_rE_Sex_Matrix <- dcast(Summary_rE_Sex, Gene_Order ~ Info, value.var = "mean_value")
Summary_rE_Sex_Matrix2 <- dcast(Summary_rE_Sex, Gene_Order ~ Info2, value.var = "mean_value")
#Summary_rE_Sex_Matrix <- Summary_rE_Sex_Matrix[genes_ordered, ]

write.table(Summary_rE_Sex_Matrix, paste(new_folder_name,"Summary_rE_Sex_Matrix.txt", sep=""), sep="\t", row.names=FALSE)
write.table(Summary_rE_Sex_Matrix2, paste(new_folder_name,"Summary_rE_Sex_Matrix2.txt", sep=""), sep="\t", row.names=FALSE)


# *Splitting the data for future uses*-----------------------------------------
## Gene data-----------------------------------------
geneFolder <- "Data_Genes separated"
createFolder(geneFolder)
gene_list <- createSeparatedDataFiles(dataFINALmean, "Gene", "Data_Genes separated")
print(gene_list)

## Subtype data -----------------------------------------
subtypeFolder <- "Data_Subtype separated"
createFolder(subtypeFolder)
subtype_list <- createSeparatedDataFiles(dataFINALmean, "Subpopulation", subtypeFolder)
print(subtype_list)

## Male/Female data-----------------------------------------

sexFolder <- "Data_Sex separated"
createFolder(sexFolder)
sex_list <- createSeparatedDataFiles(dataFINALmean, "Sex", sexFolder)
print(sex_list)
## Age & Sex Group data-----------------------------------------
groupFolder <- "Data_Group separated"
createFolder(groupFolder)
group_list <- createSeparatedDataFiles(dataFINALmean, "Group", groupFolder)
print(group_list)

unique(Female_Adult$SampleID)
length(unique(Female_Adult$SampleID))

unique(Male_Adult$SampleID)
length(unique(Male_Adult$SampleID))

unique(Female_Older$SampleID)
length(unique(Female_Older$SampleID))

unique(Male_Older$SampleID)
length(unique(Male_Older$SampleID))

# *Statistics*-----------------------------------------
## Description of the samples-----------------------------------------
createFolder("Statistics_Ttest")

pvalue_age_metadata <- metadata %>%  pairwise_t_test(Age ~ Sex,paired = FALSE)
pvalue_age_metadata$p
# Medians:-----------------------------------------
medians <- aggregate(Age ~ Sex, data = metadata, FUN = median)
medians

# Age distributions:-----------------------------------------
plot_ages <- ggplot(metadata, aes(x = Age, fill = Sex)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 30) +
  labs(title = "Age Distribution", x = "Age", y = "Frequency") +
  scale_fill_manual(values = c("#1D04C2","#A67C00")) +
  facet_wrap(~Sex, scales = "free_y")+
    theme(text = element_text(size = 14)) +
    theme(
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = "black"),
    text = element_text(color = "black"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
    )
    
# Save the plots as .png files
ggsave("Statistics_Ttest/histogram_age distribution.png", plot_ages, width = 8, height = 6, units = "in")

## LinReg FACS outliers-----------------------------------------
createFolder("Statistics_LinReg")

Age_correlation_FACS <- data.frame(Cells = character(), Sex = character(), N = numeric(), p.value = numeric(), Coefficient = numeric(), Down95= numeric(), Up95= numeric())

outliers_FACS <- FACSdata [0,]
outliers_FACS$Outlier <- rep("NA", nrow(outliers_FACS))
#i <- 7
for (i in 7:19) {
      x <- cor.test(FACSdata[,i], FACSdata$Age)
      row <- data.frame(Cells = colnames(FACSdata[i]),
                        Sex = "All",
                        N=length(FACSdata[,i]),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA, 
                        Coef_CI_Lower = NA, 
                        Coef_CI_Upper = NA 
                        )
        # Linear regression to get the coefficient and its CI
        lm_result <- lm(FACSdata[,i] ~ FACSdata$Age)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        
        # Update the row with the coefficient and its CI
        row$Coefficient <- coef_x
        row$Coef_CI_Lower <- ci[1]
        row$Coef_CI_Upper <- ci[2]
        
        
         # check for outliners
      residuals <- lm_result$residuals # Extract residuals
      residual_sd <- sd(residuals) # Calculate standard deviation of residuals
      outlier_threshold <- 3 * residual_sd # Define a threshold for outliers (e.g., 3x SD)
      outliers <- which(abs(residuals) > outlier_threshold)  # Identify outliers
      add <- FACSdata[outliers,]
      if (nrow(add) <= 0) {
        next
      }
      tryCatch({
      add$Outlier = colnames(FACSdata[i])
      outliers_FACS <- rbind(outliers_FACS, add)# Extract values corresponding to outliers
      } )
      outliers_FACS$Sex <- "All"
      Age_correlation_FACS <- rbind(Age_correlation_FACS, row)
}

Sex_correlation_FACS <- data.frame(Cells = character(), Sex = character(), N= numeric(), p.value = numeric(), Coefficient = numeric(), Down95= numeric(), Up95= numeric())

i= 7
FACS_female <- FACSdata %>% filter(FACSdata$Sex == "Female")
for (i in 7:19) {
      x_f <- cor.test(FACS_female[,i], FACS_female$Age)
      row_f <- data.frame(Cells = colnames(FACS_female[i]),  
                          Sex = "Female", 
                          N=length(FACS_female[,i]),
                          p.value = x_f$p.value,
                          Estimate = x_f$estimate,
                          Est_CI_Lower = x_f$conf.int[1],
                          Est_CI_Upper = x_f$conf.int[2],
                          Coefficient = NA, 
                          Coef_CI_Lower = NA, 
                          Coef_CI_Upper = NA 
                          )
        # Linear regression to get the coefficient and its CI
        lm_result_f <- lm(FACS_female[,i] ~ FACS_female$Age)
        coef_x_f <- summary(lm_result_f)$coefficients[2, 1]
        ci_f <- confint(lm_result_f)[2, ]
        
        # Update the row with the coefficient and its CI
        row_f$Coefficient <- coef_x_f
        row_f$Coef_CI_Lower <- ci_f[1]
        row_f$Coef_CI_Upper <- ci_f[2]
        
         Sex_correlation_FACS <- rbind(Sex_correlation_FACS, row_f)
        
        # check for outliners
      residuals_f <- lm_result_f$residuals # Extract residuals
      residual_sd_f <- sd(residuals_f) # Calculate standard deviation of residuals
      outlier_threshold_f <- 3 * residual_sd_f # Define a threshold for outliers (e.g., 3x SD)
      outliers_f <- which(abs(residuals_f) > outlier_threshold_f)  # Identify outliers
      add <- FACS_female[outliers_f,]
      if (nrow(add) <= 0) {
        next
      }
      tryCatch({
      add$Outlier = colnames(FACSdata[i])
      outliers_FACS <- rbind(outliers_FACS, add)# Extract values corresponding to outliers
      } )
}

i= 7
FACS_male <- FACSdata %>% filter(FACSdata$Sex == "Male")
for (i in 7:19) {      
      x_m <- cor.test(FACS_male[,i], FACS_male$Age)
      row_m <- data.frame(Cells = colnames(FACS_male[i]),   
                          Sex = "Male", 
                          N=length(FACS_male[,i]), 
                          p.value = x_m$p.value,
                          Estimate = x_m$estimate,
                          Est_CI_Lower = x_m$conf.int[1],
                          Est_CI_Upper = x_m$conf.int[2],
                          Coefficient = NA, 
                          Coef_CI_Lower = NA, 
                          Coef_CI_Upper = NA 
                          )
      # Linear regression to get the coefficient and its CI
        lm_result_m <- lm(FACS_male[,i] ~ FACS_male$Age)
        coef_x_m <- summary(lm_result_m)$coefficients[2, 1]
        ci_m <- confint(lm_result_m)[2, ]
        
        # Update the row with the coefficient and its CI
        row_m$Coefficient <- coef_x_m
        row_m$Coef_CI_Lower <- ci_m[1]
        row_m$Coef_CI_Upper <- ci_m[2]
        
         Sex_correlation_FACS <- rbind(Sex_correlation_FACS, row_m)
         
        # check for outliners
      residuals <- lm_result_m$residuals # Extract residuals
      residual_sd <- sd(residuals) # Calculate standard deviation of residuals
      outlier_threshold <- 3 * residual_sd # Define a threshold for outliers (e.g., 3x SD)
      outliers <- which(abs(residuals) > outlier_threshold)  # Identify outliers
      add <- FACS_male[outliers,]
      if (nrow(add) <= 0) {
        next
      }
      tryCatch({
      add$Outlier = colnames(FACSdata[i])
      outliers_FACS <- rbind(outliers_FACS, add)# Extract values corresponding to outliers
      } )
}

Correlation_FACS <- rbind(Age_correlation_FACS, Sex_correlation_FACS)

file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_Age_FACSdata.csv", sep = "")
write.csv(Age_correlation_FACS, file_name, row.names = FALSE)
file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_Sex_FACSdata.csv", sep = "")
write.csv(Sex_correlation_FACS, file_name, row.names = FALSE)
file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_FACSdata.csv", sep = "")
write.csv(Correlation_FACS, file_name, row.names = FALSE)

file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_FACSdata_Outliers.csv", sep = "")
write.csv(outliers_FACS, file_name, row.names = FALSE)

## LinReg Age Outliers-----------------------------------------
# Pearson's Correlation uses linear relationship to correlate the data
Age_correlation_rE <- data.frame(Subpopulation = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

outliers_Age_rE <-ACTB[0, ] 

i <- 1
for (i in 2:length(gene_list)) {
  df <- get(gene_list[i])
  df <- df %>% filter(!is.na(Age_Group), !is.na(dCT_Mean), !is.na(Subpopulation), !is.na(Sex))
  j <- 1
  for (j in 1:length(subtype_list)) {
    subpop <- subtype_list[j]
    df_subpop <- df %>% filter(Subpopulation == subpop)
    if (nrow(df_subpop) <= 0) {
      next
    }
    tryCatch({
      x <- cor.test(df_subpop$rE_Value, df_subpop$Age)
      row <- data.frame(
        Subpopulation = subpop,
        Gene = df_subpop$Gene[1],
        GeneID = df_subpop$GeneID[1],
        Gene_Group = df_subpop$Gene_Group[1],
        N = nrow(df_subpop),
        p.value = x$p.value,
        Estimate = x$estimate,
        Est_CI_Lower = x$conf.int[1],
        Est_CI_Upper = x$conf.int[2],
        Coefficient = NA,
        Coef_CI_Lower = NA,
        Coef_CI_Upper = NA  
      )
      
      # Linear regression to get the coefficient and its CI
      lm_result <- lm(rE_Value ~ Age, data = df_subpop)
      coef_x <- summary(lm_result)$coefficients[2, 1]
      ci <- confint(lm_result)[2, ]
      
      # Update the row with the coefficient and its CI
      row$Coefficient <- coef_x
      row$Coef_CI_Lower <- ci[1]
      row$Coef_CI_Upper <- ci[2]
      
      # check for outliners
      residuals <- lm_result$residuals # Extract residuals
      residual_sd <- sd(residuals) # Calculate standard deviation of residuals
      outlier_threshold <- 3 * residual_sd # Define a threshold for outliers (e.g., 3x SD)
      outliers <- which(abs(residuals) > outlier_threshold)  # Identify outliers
      outliers_Age_rE <- rbind(outliers_Age_rE, df_subpop[outliers,])# Extract values corresponding to outliers
      
      Age_correlation_rE <- rbind(Age_correlation_rE, row)
    }, error = function(e) {
      cat("Error in Pearson's Correlation for gene", genes[i], "& Subpopulation with Age", subpop, "- skipping this comparison.\n")
    })
    j <- j + 1
  }
  i <- i + 1
}
file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_Age_rE.csv", sep = "")
write.csv(Age_correlation_rE, file_name, row.names = FALSE)

file_name_outliers <- paste("Statistics_LinReg/", "Pearson's Correlation_Age_Outliers.csv", sep = "")
write.csv(outliers_Age_rE, file_name_outliers, row.names = FALSE)

## LinReg Sex Outliers-----------------------------------------
# Pearson's Correlation uses linear relationship to correlate the data
Age_correlation_Sex_rE <- data.frame(Subpopulation = character(), Sex = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

outliers_Sex_rE <-ACTB[0, ] 

i <- 1
for (i in 2:length(gene_list)) {
  df <- get(gene_list[i])
  for (sex in unique(df$Sex)) {
    df_sex <- df %>% filter(Sex == sex)
    for (subpop in unique(df_sex$Subpopulation)) {
      df_subpop <- df_sex %>% filter(Subpopulation == subpop)
      if (nrow(df_subpop) <= 0) {
        next
      }
      tryCatch({
        x <- cor.test(df_subpop$rE_Value, df_subpop$Age)
        row <- data.frame(
          Subpopulation = subpop,
          Sex = sex,
          Gene = df_subpop$Gene[1],
          GeneID = df_subpop$GeneID[1],
          Gene_Group = df_subpop$Gene_Group[1],
          N = nrow(df_subpop),
          p.value = x$p.value,
          Estimate = x$estimate,
          Est_CI_Lower = x$conf.int[1],
          Est_CI_Upper = x$conf.int[2],
          Coefficient = NA, 
          Coef_CI_Lower = NA, 
          Coef_CI_Upper = NA  
        )
        
        # Linear regression to get the coefficient and its CI
        lm_result <- lm(rE_Value ~ Age, data = df_subpop)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        
        # Update the row with the coefficient and its CI
        row$Coefficient <- coef_x
        row$Coef_CI_Lower <- ci[1]
        row$Coef_CI_Upper <- ci[2]
        
      # check for outliners
      residuals <- lm_result$residuals # Extract residuals
      residual_sd <- sd(residuals) # Calculate standard deviation of residuals
      outlier_threshold <- 3 * residual_sd # Define a threshold for outliers (e.g., 3x SD)
      outliers <- which(abs(residuals) > outlier_threshold)  # Identify outliers
      outliers_Sex_rE <- rbind(outliers_Sex_rE, df_subpop[outliers,])# Extract values corresponding to outliers
      
        Age_correlation_Sex_rE <- rbind(Age_correlation_Sex_rE, row)
      }, error = function(e) {
        cat("Error in Pearson's Correlation for gene ", genes[i], "& Subpopulation & Sex with Age", subpop, "- skipping this comparison.\n")
      })
    }
  }
}

file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_Age_Sex_rE.csv", sep = "")
write.csv(Age_correlation_Sex_rE, file_name, row.names = FALSE)

file_name_outliers <- paste("Statistics_LinReg/", "Pearson's Correlation_Sex_Outliers.csv", sep = "")
write.csv(outliers_Sex_rE, file_name_outliers, row.names = FALSE)

## Age Corr. w/o Outliners-----------------------------------------
createFolder("Statistics_LinReg_woOutliers")

# Pearson's Correlation uses linear relationship to correlate the data
Age_correlation_rE_woOutliers <- data.frame(Subpopulation = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

i <- 1
for (i in 2:length(gene_list)) {
  df <- get(gene_list[i])
  df <- df %>% filter(!is.na(Age_Group), !is.na(dCT_Mean), !is.na(Subpopulation), !is.na(Sex))
  j <- 1
  for (j in 1:length(subtype_list)) {
    subpop <- subtype_list[j]
    df_subpop <- df %>% filter(Subpopulation == subpop)
    df_subpop <- anti_join(df_subpop, outliers_Age_rE, by = c("SampleID", "Gene", "Subpopulation", "rE_Value"))
    if (nrow(df_subpop) <= 0) {
      next
    }
    tryCatch({
      x <- cor.test(df_subpop$rE_Value, df_subpop$Age)
      row <- data.frame(
        Subpopulation = subpop,
        Gene = df_subpop$Gene[1],
        GeneID = df_subpop$GeneID[1],
        Gene_Group = df_subpop$Gene_Group[1],
        N = nrow(df_subpop),
        p.value = x$p.value,
        Estimate = x$estimate,
        Est_CI_Lower = x$conf.int[1],
        Est_CI_Upper = x$conf.int[2],
        Coefficient = NA,
        Coef_CI_Lower = NA,
        Coef_CI_Upper = NA  
      )
      
      # Linear regression to get the coefficient and its CI
      lm_result <- lm(rE_Value ~ Age, data = df_subpop)
      coef_x <- summary(lm_result)$coefficients[2, 1]
      ci <- confint(lm_result)[2, ]
      
      # Update the row with the coefficient and its CI
      row$Coefficient <- coef_x
      row$Coef_CI_Lower <- ci[1]
      row$Coef_CI_Upper <- ci[2]
      
      Age_correlation_rE_woOutliers <- rbind(Age_correlation_rE_woOutliers, row)
    }, error = function(e) {
      cat("Error in Pearson's Correlation for gene", genes[i], "& Subpopulation with Age", subpop, "- skipping this comparison.\n")
    })
    j <- j + 1
  }
  i <- i + 1
}
file_name <- paste("Statistics_LinReg_woOutliers/", "Pearson's Correlation_Age_rE_woOutliers.csv", sep = "")
write.csv(Age_correlation_rE_woOutliers, file_name, row.names = FALSE)

## Sex Corr. w/o Outliners-----------------------------------------
createFolder("Statistics_LinReg_woOutliers")

# Pearson's Correlation uses linear relationship to correlate the data
Age_correlation_Sex_rE_woOutliers <- data.frame(Subpopulation = character(), Gene = character(), GeneID = numeric(), Gene_Group = character(), N = numeric(), p.value = numeric(), Estimate = numeric(), Est_CI_Lower = numeric(), Est_CI_Upper = numeric(), Coefficient = numeric(), Coef_CI_Lower = numeric(), Coef_CI_Upper = numeric())

i <- 1
for (i in 2:length(gene_list)) {
  df <- get(gene_list[i])
  for (sex in unique(df$Sex)) {
    df_sex <- df %>% filter(Sex == sex)
    for (subpop in unique(df_sex$Subpopulation)) {
      df_subpop <- df_sex %>% filter(Subpopulation == subpop)
      df_subpop <- anti_join(df_subpop, outliers_Sex_rE, by = c("SampleID", "Gene", "Subpopulation", "rE_Value"))
      if (nrow(df_subpop) <= 0) {
        next
      }
      tryCatch({
        x <- cor.test(df_subpop$rE_Value, df_subpop$Age)
        row <- data.frame(
          Subpopulation = subpop,
          Sex = sex,
          Gene = df_subpop$Gene[1],
          GeneID = df_subpop$GeneID[1],
          Gene_Group = df_subpop$Gene_Group[1],
          N = nrow(df_subpop),
          p.value = x$p.value,
          Estimate = x$estimate,
          Est_CI_Lower = x$conf.int[1],
          Est_CI_Upper = x$conf.int[2],
          Coefficient = NA, 
          Coef_CI_Lower = NA, 
          Coef_CI_Upper = NA  
        )
        
        # Linear regression to get the coefficient and its CI
        lm_result <- lm(rE_Value ~ Age, data = df_subpop)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        
        # Update the row with the coefficient and its CI
        row$Coefficient <- coef_x
        row$Coef_CI_Lower <- ci[1]
        row$Coef_CI_Upper <- ci[2]
        
        Age_correlation_Sex_rE_woOutliers <- rbind(Age_correlation_Sex_rE_woOutliers, row)
      }, error = function(e) {
        cat("Error in Pearson's Correlation for gene ", genes[i], "& Subpopulation & Sex with Age", subpop, "- skipping this comparison.\n")
      })
    }
  }
}

file_name <- paste("Statistics_LinReg_woOutliers/", "Pearson's Correlation_Age_Sex_rE_woOutliers.csv", sep = "")
write.csv(Age_correlation_Sex_rE_woOutliers, file_name, row.names = FALSE)

## FACS Corr. w/o Outliers-----------------------------------------
createFolder("Plots_LinReg_FACS_woOutliers")

Age_correlation_FACS_woOutliers <- data.frame(Cells = character(), Sex = character(), N = numeric(), p.value = numeric(), Coefficient = numeric(), Down95= numeric(), Up95= numeric())

FACSdata_woOutliers <- FACSdata

i=1
for (i in 1:nrow(outliers_FACS)) {
    EventCondition <- outliers_FACS$Outlier[i]
    SampleCondition <- outliers_FACS$SampleID[i]
   if (outliers_FACS$Sex[i] == "All") {
    FACSdata_woOutliers[[EventCondition]][FACSdata_woOutliers$SampleID == SampleCondition] <- NA
   }
    i<-i+1
}

#i <- 7
for (i in 7:19) {
      x <- cor.test(FACSdata_woOutliers[,i], FACSdata_woOutliers$Age)
      row <- data.frame(Cells = colnames(FACSdata_woOutliers[i]),
                        Sex = "All",
                        N=length(FACSdata_woOutliers[,i]),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA, 
                        Coef_CI_Lower = NA, 
                        Coef_CI_Upper = NA 
                        )
        # Linear regression to get the coefficient and its CI
        lm_result <- lm(FACSdata_woOutliers[,i] ~ FACSdata_woOutliers$Age)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        
        # Update the row with the coefficient and its CI
        row$Coefficient <- coef_x
        row$Coef_CI_Lower <- ci[1]
        row$Coef_CI_Upper <- ci[2]
      Age_correlation_FACS_woOutliers <- rbind(Age_correlation_FACS_woOutliers, row)
}

Sex_correlation_FACS_woOutliers <- data.frame(Cells = character(), Sex = character(), N= numeric(), p.value = numeric(), Coefficient = numeric(), Down95= numeric(), Up95= numeric())

FACSdata_Sex_woOutliers <- FACSdata

i=1
for (i in 1:nrow(outliers_FACS)) {
    EventCondition <- outliers_FACS$Outlier[i]
    SampleCondition <- outliers_FACS$SampleID[i]
    if (outliers_FACS$Sex[i] != "All") {
    FACSdata_Sex_woOutliers[[EventCondition]][FACSdata_Sex_woOutliers$SampleID == SampleCondition] <- NA
    }
}

FACS_female_woOutliers <- FACSdata_Sex_woOutliers %>% filter(FACSdata_Sex_woOutliers$Sex == "Female")
FACS_male_woOutliers <- FACSdata_Sex_woOutliers %>% filter(FACSdata_Sex_woOutliers$Sex == "Male")

for (i in 7:19) {
      x_f <- cor.test(FACS_female_woOutliers[,i], FACS_female_woOutliers$Age)
      row_f <- data.frame(Cells = colnames(FACS_female_woOutliers[i]),  
                          Sex = "Female", 
                          N=length(FACS_female_woOutliers[,i]),
                          p.value = x_f$p.value,
                          Estimate = x_f$estimate,
                          Est_CI_Lower = x_f$conf.int[1],
                          Est_CI_Upper = x_f$conf.int[2],
                          Coefficient = NA, 
                          Coef_CI_Lower = NA, 
                          Coef_CI_Upper = NA 
                          )
        # Linear regression to get the coefficient and its CI
        lm_result_f <- lm(FACS_female_woOutliers[,i] ~ FACS_female_woOutliers$Age)
        coef_x_f <- summary(lm_result_f)$coefficients[2, 1]
        ci_f <- confint(lm_result_f)[2, ]
        
        # Update the row with the coefficient and its CI
        row_f$Coefficient <- coef_x_f
        row_f$Coef_CI_Lower <- ci_f[1]
        row_f$Coef_CI_Upper <- ci_f[2]
        
      x_m <- cor.test(FACS_male_woOutliers[,i], FACS_male_woOutliers$Age)
      row_m <- data.frame(Cells = colnames(FACS_male_woOutliers[i]),   
                          Sex = "Male", 
                          N=length(FACS_male[,i]), 
                          p.value = x_m$p.value,
                          Estimate = x_m$estimate,
                          Est_CI_Lower = x_m$conf.int[1],
                          Est_CI_Upper = x_m$conf.int[2],
                          Coefficient = NA, 
                          Coef_CI_Lower = NA, 
                          Coef_CI_Upper = NA 
                          )
      # Linear regression to get the coefficient and its CI
        lm_result_m <- lm(FACS_male_woOutliers[,i] ~ FACS_male_woOutliers$Age)
        coef_x_m <- summary(lm_result_m)$coefficients[2, 1]
        ci_m <- confint(lm_result_m)[2, ]
        
        # Update the row with the coefficient and its CI
        row_m$Coefficient <- coef_x_m
        row_m$Coef_CI_Lower <- ci_m[1]
        row_m$Coef_CI_Upper <- ci_m[2]
      Sex_correlation_FACS_woOutliers <- rbind(Sex_correlation_FACS_woOutliers, row_f, row_m)
}

Correlation_FACS_woOutliers <- rbind(Age_correlation_FACS_woOutliers, Sex_correlation_FACS_woOutliers)

file_name <- paste("Statistics_LinReg_woOutliers/", "Pearson's Correlation_Age_FACSdata.csv", sep = "")
write.csv(Age_correlation_FACS_woOutliers, file_name, row.names = FALSE)
file_name <- paste("Statistics_LinReg_woOutliers/", "Pearson's Correlation_Sex_FACSdata.csv", sep = "")
write.csv(Sex_correlation_FACS_woOutliers, file_name, row.names = FALSE)
file_name <- paste("Statistics_LinReg_woOutliers/", "Pearson's Correlation_FACSdata.csv", sep = "")
write.csv(Correlation_FACS_woOutliers, file_name, row.names = FALSE)

# Define variables and titles for each plot
plotVariables <- list(
  list(xvar = "Age", yvar = "Alive", title = "LinReg of alive PBMCs with age", filename = "Alive_Plot"),
  list(xvar = "Age", yvar = "Bcells", title = "LinReg of Bcells with age", filename = "Bcells_Plot"),
  list(xvar = "Age", yvar = "NKcells", title = "LinReg of NKcells with age", filename = "NKcells_Plot"),
  list(xvar = "Age", yvar = "Tcells", title = "LinReg of Tcells with age", filename = "Tcells_Plot"),
  list(xvar = "Age", yvar = "Monocytes", title = "LinReg of Monocytes with age", filename = "Monocytes_Plot"),
  list(xvar = "Age", yvar = "Classical_Monocytes", title = "LinReg of classical monocytes with age", filename = "Classical_Monocytes_Plot"),
  list(xvar = "Age", yvar = "Intermediate_Monocytes", title = "LinReg of intermediate monocytes with age", filename = "Intermediate_Monocytes_Plot"),
  list(xvar = "Age", yvar = "NonClassical_Monocytes", title = "LinReg of non-classical monocytes with age", filename = "NonClassical_Monocytes_Plot")
)

# Loop through plotVariables and generate/save the plots
for (plotVar in plotVariables) {
  generateLinRegPlot(FACSdata_woOutliers, plotVar$xvar, plotVar$yvar, plotVar$title, "Plots_LinReg_FACS_woOutliers/", plotVar$filename)
    generateLinRegPlot_Sex(FACSdata_Sex_woOutliers, plotVar$xvar, plotVar$yvar, plotVar$title, "Plots_LinReg_FACS_woOutliers/Sex_", plotVar$filename)
}
## Plot all linear regressions-----------------------------------------
new_folder_name <- "Plot_LinReg"
createFolder(new_folder_name)

i <- 1
pop <- unique(dataFINALmean$Subpopulation)
for (i in 1:length(gene_list)) {
  df <- get(gene_list[i])
  sig <- gene_list[i]
  j <- 1
  for (j in 1:length(pop)) {
    sub <- pop[j]
    df_filtered <- filter(df, Subpopulation == sub)
    
    # Check if there is data for the current subpopulation
    if (nrow(df_filtered) > 0) {
      titel <- paste(sig, " expression in ", sub, " monocytes (LinReg, rE, combined)", sep="")
      file_name <- paste(new_folder_name, "/", titel, ".png", sep = "")
      
      ggplot(data = df_filtered, aes(x = Age, y = rE_Value)) +
        geom_smooth(method = "glm", color = "black") +
        geom_point(aes(color = Sex), size = 2) +
        ggtitle(titel) +
        xlab("Age") +
        ylab("Gene expression relative to B2M") +
        #facet_wrap(~Subpopulation, ncol = 2) +
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
      
      ggsave(filename = file_name)
      
      titel <- paste(sig, " expression in ",sub," monocytes (Sex separated) (LinReg, rE)", sep="")
file_name <- paste(new_folder_name, "/",titel, ".png", sep = "")
#png(filename=file_name, height=750, width=750)
ggplot(data =  df_filtered, aes(x = Age, y = rE_Value, colour = Sex)) + 
    geom_smooth(method = glm) + 
    geom_point(size = 2) +
    ggtitle(titel) +
    xlab("Age") + 
    ylab("Gene expression relative to B2M") +
    #facet_wrap(~ Subpopulation, ncol = 2) +
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
    }
    j = j + 1
  }
  i = i + 1
}

## Replace the variables with data excluding the Outliers-----------------------------------------
Age_correlation_rE <- Age_correlation_rE_woOutliers
Age_correlation_Sex_rE <- Age_correlation_Sex_rE_woOutliers

## Plot all linear regressions w/o Outliers-----------------------------------------
new_folder_name <- "Plot_LinReg_woOutliers"
createFolder(new_folder_name)

i <- 1
pop <- unique(dataFINALmean$Subpopulation)
for (i in 1:length(gene_list)) {
  df <- get(gene_list[i])
  sig <- gene_list[i]
  j <- 1
  for (j in 1:length(pop)) {
    sub <- pop[j]
    df_filtered <- filter(df, Subpopulation == sub)
    #remove Age correlation outliers
    df_filtered_age <- anti_join(df_filtered, outliers_Age_rE, by = c("SampleID", "Gene", "Subpopulation", "rE_Value"))
    # Check if there is data for the current subpopulation
    if (nrow(df_filtered_age) > 0) {
      titel <- paste(sig, " expression in ", sub, " monocytes (LinReg, rE, combined)", sep="")
      file_name <- paste(new_folder_name, "/", titel, ".png", sep = "")
      
      ggplot(data = df_filtered_age, aes(x = Age, y = rE_Value)) +
        geom_smooth(method = "glm", color = "black") +
        geom_point(aes(color = Sex), size = 2) +
        ggtitle(titel) +
        xlab("Age") +
        ylab("Gene expression relative to B2M") +
        #facet_wrap(~Subpopulation, ncol = 2) +
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
      
      ggsave(filename = file_name)
      
      #remove Sex correlation outliers
      df_filtered_sex <- anti_join(df_filtered, outliers_Sex_rE, by = c("SampleID", "Gene", "Subpopulation", "rE_Value"))
      
      titel <- paste(sig, " expression in ",sub," monocytes (Sex separated) (LinReg, rE)", sep="")
file_name <- paste(new_folder_name, "/",titel, ".png", sep = "")
#png(filename=file_name, height=750, width=750)
ggplot(data =  df_filtered_sex, aes(x = Age, y = rE_Value, colour = Sex)) + 
    geom_smooth(method = glm) + 
    geom_point(size = 2) +
    ggtitle(titel) +
    xlab("Age") + 
    ylab("Gene expression relative to B2M") +
    #facet_wrap(~ Subpopulation, ncol = 2) +
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
    }
    j = j + 1
  }
  i = i + 1
}

## Matrix / Heatmap of Coefficient (not separated by Sex)-----------------------------------------
Pearsson_Matrix <- acast(Age_correlation_rE_woOutliers, Gene  ~ Subpopulation,value.var = "Coefficient" )

file_name <- paste("Statistics_LinReg/", "Pearsson's Matrix_Age_rE.csv", sep = "")
write.csv(Pearsson_Matrix, file_name)

Age_correlation_Sex_rE_woOutliers$Group <- paste(Age_correlation_Sex_rE_woOutliers$Subpopulation, Age_correlation_Sex_rE_woOutliers$Sex, sep="_")
Age_correlation_rE_woOutliers$Group <- paste(Age_correlation_rE_woOutliers$Subpopulation, Age_correlation_rE_woOutliers$Sex, sep="_")
Age_correlation_Sex_rE_woOutliers$Gene_Order <- paste(Age_correlation_Sex_rE_woOutliers$GeneID, Age_correlation_Sex_rE_woOutliers$Gene, sep="_")
Age_correlation_rE_woOutliers$Gene_Order <- paste(Age_correlation_rE_woOutliers$GeneID, Age_correlation_rE_woOutliers$Gene, sep="_")
Age_correlation_Sex_rE_woOutliers$Info <- paste(Age_correlation_Sex_rE_woOutliers$Sex, Age_correlation_Sex_rE_woOutliers$Subpopulation, sep="_")
Age_correlation_rE_woOutliers$Info <- paste(Age_correlation_rE_woOutliers$Sex, Age_correlation_rE_woOutliers$Subpopulation, sep="_")

Pearsson_Sex_Matrix_woOutliers <- acast(Age_correlation_Sex_rE_woOutliers, Gene  ~ Group,value.var = "Coefficient" )

file_name <- paste("Statistics_LinReg/", "Pearsson's Matrix_Age_Sex_rE_woOutliers.csv", sep = "")
write.csv(Pearsson_Sex_Matrix_woOutliers, file_name)

Age_correlation_rE_woOutliers$Sex <- "All"
All_Age_Pearson_rE <- rbind(Age_correlation_Sex_rE_woOutliers, Age_correlation_rE_woOutliers)
file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_Age_All_rE_woOutliers.csv", sep = "")
write.csv(All_Age_Pearson_rE, file_name)

Age_correlation_rE_Matrix <- dcast(Age_correlation_rE_woOutliers, Gene_Order ~ Subpopulation, value.var = "Coefficient")
Age_correlation_Sex_rE_Matrix <- dcast(Age_correlation_Sex_rE_woOutliers, Gene_Order ~ Info, value.var = "Coefficient")

write.table(Age_correlation_rE_Matrix, "Statistics_LinReg/Age_Pearson-Coefficient_rE_Matrix_woOutliers.txt", sep="\t", row.names=FALSE)
write.table(Age_correlation_Sex_rE_Matrix, "Statistics_LinReg/Age_Pearson-Coefficient_Sex_rE_Matrix_woOutliers.txt", sep="\t", row.names=FALSE)

Age_correlation_rE_Matrix_est <- dcast(Age_correlation_rE_woOutliers, Gene_Order ~ Subpopulation, value.var = "Estimate")
Age_correlation_Sex_rE_Matrix_est <- dcast(Age_correlation_Sex_rE_woOutliers, Gene_Order ~ Info, value.var = "Estimate")

write.table(Age_correlation_rE_Matrix_est, "Statistics_LinReg/Age_Pearson-Estimate_rE_Matrix_woOutliers.txt", sep="\t", row.names=FALSE)
write.table(Age_correlation_Sex_rE_Matrix_est, "Statistics_LinReg/Age_Pearson-Estimate_Sex_rE_Matrix_woOutliers.txt", sep="\t", row.names=FALSE)

#*Significant results* -----------------------------------------
## Extract significant Data (rE):-----------------------------------------

# Pearson Correlation:
Age_sigGenes_Pearson <- All_Age_Pearson_rE %>% filter(All_Age_Pearson_rE$p.value <0.05)
file_name <- paste("Statistics_LinReg/", "Pearsson's_All_significant_rE_woOutliners.csv", sep = "")
write.csv(Age_sigGenes_Pearson, file_name)

## Extract Data of Age significant Genes-----------------------------------------
i <- 1
j <- length(Age_sigGenes_Pearson)
for (i in 2:j) {
    sub <- Age_sigGenes_Pearson$Subpopulation[i]
    df <- get(Age_sigGenes_Pearson$Gene[i])
    df <- filter(df, Subpopulation == sub)
     }

df <- filter(df, !is.na(Age_Group), !is.na(rE_Value), !is.na(Sex))

## LinearRegression (only significant Gene)-----------------------------------------
createFolder("Significant Genes relating to Age")

i <- 1
j <- length(Age_sigGenes_Pearson$Group)
for (i in 1:j) {
    sub <- Age_sigGenes_Pearson$Subpopulation[i]
    sig <- Age_sigGenes_Pearson$Gene[i]
    df <- get(sig)
    df <- filter(df, Subpopulation == sub)
    s <- Age_sigGenes_Pearson$Sex[i]
    if (s=="All"){
        titel <- paste(sig, " expression in ", sub, " monocytes (LinReg, rE, combined)", sep="")
        #remove Age correlation outliers
    df <- anti_join(df, outliers_Age_rE, by = c("SampleID", "Gene", "Subpopulation", "rE_Value"))
file_name <- paste("Significant Genes relating to Age/",titel, ".png", sep = "")
    ggplot(data = df, aes(x = Age, y = rE_Value)) +
  geom_smooth(method = "glm", color = "black") +
  geom_point(aes(color = Sex), size = 2) +
  ggtitle(titel) +
  xlab("Age") +
  ylab("Gene expression relative to B2M") +
  facet_wrap(~Subpopulation, ncol = 2) +
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
        #remove Sex correlation outliers
      df <- anti_join(df, outliers_Sex_rE, by = c("SampleID", "Gene", "Subpopulation", "rE_Value"))
    titel <- paste(sig, " expression in ",sub," monocytes (",s,") (LinReg, rE)", sep="")
file_name <- paste("Significant Genes relating to Age/",titel, ".png", sep = "")
#png(filename=file_name, height=750, width=750)
ggplot(data = df, aes(x = Age, y = rE_Value, colour = Sex)) + 
    geom_smooth(method = glm) + 
    geom_point(size = 2) +
    ggtitle(titel) +
    xlab("Age") + 
    ylab("Gene expression relative to B2M") +
    facet_wrap(~ Subpopulation, ncol = 2) +
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

## LinReg for FACS data-----------------------------------------
createFolder("Plots_LinReg_FACS")

Age_correlation_FACS <- data.frame(Cells = character(), Sex = character(), N = numeric(), p.value = numeric(), Coefficient = numeric(), Down95= numeric(), Up95= numeric())

#i <- 7
for (i in 7:19) {
      x <- cor.test(FACSdata[,i], FACSdata$Age)
      row <- data.frame(Cells = colnames(FACSdata[i]),
                        Sex = "All",
                        N=length(FACSdata[,i]),
                        p.value = x$p.value,
                        Estimate = x$estimate,
                        Est_CI_Lower = x$conf.int[1],
                        Est_CI_Upper = x$conf.int[2],
                        Coefficient = NA, 
                        Coef_CI_Lower = NA, 
                        Coef_CI_Upper = NA 
                        )
        # Linear regression to get the coefficient and its CI
        lm_result <- lm(FACSdata[,i] ~ FACSdata$Age)
        coef_x <- summary(lm_result)$coefficients[2, 1]
        ci <- confint(lm_result)[2, ]
        
        # Update the row with the coefficient and its CI
        row$Coefficient <- coef_x
        row$Coef_CI_Lower <- ci[1]
        row$Coef_CI_Upper <- ci[2]
      Age_correlation_FACS <- rbind(Age_correlation_FACS, row)
}

Sex_correlation_FACS <- data.frame(Cells = character(), Sex = character(), N= numeric(), p.value = numeric(), Coefficient = numeric(), Down95= numeric(), Up95= numeric())

FACS_female <- FACSdata %>% filter(FACSdata$Sex == "Female")
FACS_male <- FACSdata %>% filter(FACSdata$Sex == "Male")
for (i in 7:19) {
      x_f <- cor.test(FACS_female[,i], FACS_female$Age)
      row_f <- data.frame(Cells = colnames(FACS_female[i]),  
                          Sex = "Female", 
                          N=length(FACS_female[,i]),
                          p.value = x_f$p.value,
                          Estimate = x_f$estimate,
                          Est_CI_Lower = x_f$conf.int[1],
                          Est_CI_Upper = x_f$conf.int[2],
                          Coefficient = NA, 
                          Coef_CI_Lower = NA, 
                          Coef_CI_Upper = NA 
                          )
        # Linear regression to get the coefficient and its CI
        lm_result_f <- lm(FACS_female[,i] ~ FACS_female$Age)
        coef_x_f <- summary(lm_result_f)$coefficients[2, 1]
        ci_f <- confint(lm_result_f)[2, ]
        
        # Update the row with the coefficient and its CI
        row_f$Coefficient <- coef_x_f
        row_f$Coef_CI_Lower <- ci_f[1]
        row_f$Coef_CI_Upper <- ci_f[2]
        
      x_m <- cor.test(FACS_male[,i], FACS_male$Age)
      row_m <- data.frame(Cells = colnames(FACS_male[i]),   
                          Sex = "Male", 
                          N=length(FACS_male[,i]), 
                          p.value = x_m$p.value,
                          Estimate = x_m$estimate,
                          Est_CI_Lower = x_m$conf.int[1],
                          Est_CI_Upper = x_m$conf.int[2],
                          Coefficient = NA, 
                          Coef_CI_Lower = NA, 
                          Coef_CI_Upper = NA 
                          )
      # Linear regression to get the coefficient and its CI
        lm_result_m <- lm(FACS_male[,i] ~ FACS_male$Age)
        coef_x_m <- summary(lm_result_m)$coefficients[2, 1]
        ci_m <- confint(lm_result_m)[2, ]
        
        # Update the row with the coefficient and its CI
        row_m$Coefficient <- coef_x_m
        row_m$Coef_CI_Lower <- ci_m[1]
        row_m$Coef_CI_Upper <- ci_m[2]
      Sex_correlation_FACS <- rbind(Sex_correlation_FACS, row_f, row_m)
}

Correlation_FACS <- rbind(Age_correlation_FACS, Sex_correlation_FACS)

file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_Age_FACSdata.csv", sep = "")
write.csv(Age_correlation_FACS, file_name, row.names = FALSE)
file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_Sex_FACSdata.csv", sep = "")
write.csv(Sex_correlation_FACS, file_name, row.names = FALSE)
file_name <- paste("Statistics_LinReg/", "Pearson's Correlation_FACSdata.csv", sep = "")
write.csv(Correlation_FACS, file_name, row.names = FALSE)

# Define variables and titles for each plot
plotVariables <- list(
  list(xvar = "Age", yvar = "Alive", title = "LinReg of alive PBMCs with age", filename = "Alive_Plot"),
  list(xvar = "Age", yvar = "Bcells", title = "LinReg of Bcells with age", filename = "Bcells_Plot"),
  list(xvar = "Age", yvar = "NKcells", title = "LinReg of NKcells with age", filename = "NKcells_Plot"),
  list(xvar = "Age", yvar = "Tcells", title = "LinReg of Tcells with age", filename = "Tcells_Plot"),
  list(xvar = "Age", yvar = "Monocytes", title = "LinReg of Monocytes with age", filename = "Monocytes_Plot"),
  list(xvar = "Age", yvar = "Classical_Monocytes", title = "LinReg of classical monocytes with age", filename = "Classical_Monocytes_Plot"),
  list(xvar = "Age", yvar = "Intermediate_Monocytes", title = "LinReg of intermediate monocytes with age", filename = "Intermediate_Monocytes_Plot"),
  list(xvar = "Age", yvar = "NonClassical_Monocytes", title = "LinReg of non-classical monocytes with age", filename = "NonClassical_Monocytes_Plot")
)

# Loop through plotVariables and generate/save the plots
for (plotVar in plotVariables) {
  generateLinRegPlot(FACSdata, plotVar$xvar, plotVar$yvar, plotVar$title, "Plots_LinReg_FACS/", plotVar$filename)
    generateLinRegPlot_Sex(FACSdata, plotVar$xvar, plotVar$yvar, plotVar$title, "Plots_LinReg_FACS/Sex_", plotVar$filename)
}
