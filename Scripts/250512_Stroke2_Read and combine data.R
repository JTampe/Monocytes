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
setwd("/Users/ju5263ta/Github/Monocytes/rawData_Stroke")
metadataP <- read_excel("Metadata_PatientID.xlsx")
Metadata_GeneID <- read_excel("Metadata_GeneID.xlsx")
Plate_info <- read_excel("241014_PlateID.xlsx")

Metadata_NeuroTest <- as.data.frame(read_excel("Metadata_NeuroTest.xlsx", 
                                               col_types = c("text", "text", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric")))

# Create NHISS_End column based on TP4 -------------------------
Metadata_NeuroTest$NHISS_End <- NA

# Get TP4 NHISS values
tp4_values <- Metadata_NeuroTest %>%
    filter(Timepoint == "TP4") %>%
    select(SampleID, NHISS) %>%
    rename(NHISS_TP4 = NHISS)

# Join TP4 NHISS values back to all timepoints
Metadata_NeuroTest <- Metadata_NeuroTest %>%
    left_join(tp4_values, by = "SampleID") %>%
    mutate(NHISS_End = NHISS_TP4) %>%
    select(-NHISS_TP4)

# Compute NHISS_Ratio for SampleIDs with both TP1 and TP4
nhiss_ratios <- Metadata_NeuroTest %>%
    filter(Timepoint %in% c("TP1", "TP4")) %>%
    group_by(SampleID) %>%
    summarise(
        NHISS_TP1 = NHISS[Timepoint == "TP1"][1],
        NHISS_TP4 = NHISS[Timepoint == "TP4"][1],
        NHISS_Ratio = if (!is.na(NHISS_TP1) && NHISS_TP1 != 0) {
            (NHISS_TP1 - NHISS_TP4) / NHISS_TP1
        } else {
            NA
        },
        .groups = "drop"
    ) %>%
    select(SampleID, NHISS_Ratio)

# Join NHISS_Ratio back to Metadata_NeuroTest
Metadata_NeuroTest <- Metadata_NeuroTest %>%
    left_join(nhiss_ratios, by = "SampleID")

# Compute NHISS_Ratio for SampleIDs with both TP1 and TP4
nhiss_diff <- Metadata_NeuroTest %>%
    filter(Timepoint %in% c("TP1", "TP4")) %>%
    group_by(SampleID) %>%
    summarise(
        NHISS_TP1 = NHISS[Timepoint == "TP1"][1],
        NHISS_TP4 = NHISS[Timepoint == "TP4"][1],
        NHISS_Diff = if (!is.na(NHISS_TP1) && NHISS_TP1 != 0) {
            (NHISS_TP1 - NHISS_TP4)
        } else {
            NA
        },
        .groups = "drop"
    ) %>%
    select(SampleID, NHISS_Diff)

# Join NHISS_Ratio back to Metadata_NeuroTest
Metadata_NeuroTest <- Metadata_NeuroTest %>%
    left_join(nhiss_diff, by = "SampleID")

dataAll <- merge(dataEXPAND, metadataP, by = "SampleID", all.x = TRUE)
#### match the right controls
#Matched_TP0_Gene <- c("Ctr01", "Ctr02", "Ctr08", "Ctr09", "Ctr11", "Ctr14", "Ctr16", "Ctr17", "Ctr18", "Ctr25", "Ctr27", "Ctr29", "Ctr31","Ctr41", "Ctr43")
Matched_TP0_Gene <- c("Ctr20", "Ctr07", "Ctr10", "Ctr14", "Ctr25", "Ctr24", "Ctr18", "Ctr19", 
                      "Ctr08", "Ctr26", "Ctr17", "Ctr09", "Ctr27","Ctr16", "Ctr11")

Unmatched_TP0_FACS <- c("Ctr16", "Ctr21", "Ctr23", "Ctr30", "Ctr39", "Ctr44", "Ctr33")

# FACS
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

# add the metadata of NHISS score of Patients
FACSdata <- merge(FACSdata, Metadata_NeuroTest, by = c("SampleID", "Timepoint"), all.x = TRUE)

# assign sex
FACSdata$Sex[grep("M",FACSdata$Sex)]  <- "Male"
FACSdata$Sex[grep("F",FACSdata$Sex)]  <- "Female"


# Combine the data frames based on common values in the "id" column
# --> here i loose all the test because they are not in the metadata
dataAll <- merge(dataEXPAND, metadataP, by = "SampleID", all.x = TRUE)

#relabel Male & Female
dataAll$Sex[grep("M",dataAll$Sex)]  <- "Male"
dataAll$Sex[grep("F",dataAll$Sex)]  <- "Female"

dataAll <- dataAll %>% filter(!is.na(Age), !is.na(SampleID), !is.na(Category), !is.na(Sex))

#### failed Genes --------------------------------------------------------
fails <- dataAll %>%
    group_by(Subpopulation, Gene, Timepoint) %>%
    summarize(
        N = n(),
        Count_Above_35_CT = sum(CT_Value >= 35, na.rm = TRUE) + sum(is.na(CT_Value)),  # Include NAs in the count
    )

fails$Per_Above_35_CT = fails$Count_Above_35_CT/fails$N * 100

failed_genes <- fails %>% filter(Per_Above_35_CT >= 90) %>% select(Gene)
wrong_threshold <- c("CLEC7A", "S100A8", "S100A9")
exclude_genes  <- unique(c("Xeno",failed_genes$Gene, wrong_threshold))
paste("The gene",unique(failed_genes$Gene), "was excluded, as more than 90% of the reads failed.")
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

#### Evaluate Housekeeping genes -----------
FailedHK <-table(c( 
    getMeTheFailedSamples(dataWORKING, "ACTB", 30),
    getMeTheFailedSamples(dataWORKING, "B2M", 30),
    getMeTheFailedSamples(dataWORKING, "GAPDH", 30)))

hist(as.numeric(FailedHK), 
     main = "Histogram of Failed Samples", 
     xlab = "Failed Samples Count", 
     breaks = 10)

failedSAMPLES <- names(FailedHK) [ which(FailedHK > 2) ]

# Set CT & dCT NA for the Sampels with failed HK gene
# Set data$Value to NA for the samples that are in the failedSAMPLES vector
dataWORKING$CT_Value <- ifelse(dataWORKING$Sample %in% failedSAMPLES, yes = NA, no = dataWORKING$CT_Value)
dataWORKING$dCT_Value <- ifelse(dataWORKING$Sample %in% failedSAMPLES, yes = NA, no = dataWORKING$dCT_Value)
dataWORKING$Comment <- ifelse(dataWORKING$Sample %in% failedSAMPLES, yes = "More than 2 houskeeing gene reads failed", no = dataWORKING$Comment)

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

# *Test: Pool TP1 & TP2 and TP3 & TP4* --------------------------------------------------------------------------------------------
# add TP2 to TP1
# dataWORKING$Timepoint[grep("TP2",dataWORKING$Timepoint)]  <- "TP1"
# add TP2 to TP1
# dataWORKING$Timepoint[grep("TP4",dataWORKING$Timepoint)]  <- "TP3"


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
dataWORKING$Recovery <- as.factor(dataWORKING$Recovery) # check the effect of ordered
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
    group_by(SampleID, Gene, Subpopulation,Timepoint, TechnicalReplicate, BiologicalReplicate, Sample, Gene1, Age, Sex, Category, Recovery) %>%    
    summarize(CT = mean(CT_Value, na.rm = TRUE),
              # CT_sd = sd(CT_Value, na.rm = TRUE), 
              # CT_N = sum(!is.na(CT_Value)), 
              # dCT_B2M = mean(dCT_Value, na.rm = TRUE), 
              # dCT_B2M_sd = sd(dCT_Value, na.rm = TRUE),
              .groups = "drop")

dataFINAL <- as.data.frame(data_imputed_mean)

# *Normalize the data vs the  Geometric mean of each Sample* ----------------

# Calculate the CT median or mean
dataFINAL$CTmedian <- NA
#dataFINAL$CTmean <- NA
#dataFINAL$CTgmean <- NA  # For geometric mean

dataFINAL <- dataFINAL %>%
    group_by(Sample) %>%
    mutate(
        CTmedian = median(CT[CT < 35], na.rm = TRUE),
        #CTmean = mean(CT[CT < 35], na.rm = TRUE),
        #CTgmean = exp(mean(log(CT[CT < 35]), na.rm = TRUE))  # Geometric mean
    )
# normalize all CT value to the Media CT value of each Sample
dataFINAL$dCT_median <- dataFINAL$CT - dataFINAL$CTmedian
#dataFINAL$dCT_mean <- dataFINAL$CT - dataFINAL$CTmean
#dataFINAL$dCT_gmean <- dataFINAL$CT - dataFINAL$CTgmean

# *Zscore: Normalize the data vs Gene expression ----------------
# as such: apply( data_summary , 1 , function(x) { (x["mean_Z"]- x["MeanForGene"])/ x["sd_for_gene"] } )
dataFINAL <- dataFINAL %>%
    group_by(Gene) %>%
    mutate(
        mean_dCT_Gene = mean(dCT_median, na.rm = TRUE),  # Mean of mean_capZ by Gene
        sd_dCT_Gene = sd(dCT_median, na.rm = TRUE),      # Standard deviation of mean_capZ by Gene
        z_score_dCT = (dCT_median - mean_dCT_Gene) / sd_dCT_Gene,  # Z-score calculation
        inverted_z_score = -z_score_dCT,
        capped_z_score = pmin(pmax(inverted_z_score, -3), 3) 
    ) %>%
    ungroup() 

#### Gene & NHISS Metadata  ----------------
## Gene Metadata & other edits:
#Gene_Info <- data.frame(Gene = Metadata_GeneID$Gene, GeneID = Metadata_GeneID$GeneID, Gene_Group = Metadata_GeneID$Gene_Group)
#dataFINAL <- merge(dataFINAL, Gene_Info, by = "Gene")
dataFINAL <- merge(dataFINAL, Metadata_GeneID, by = "Gene")

# add the metadata of NHISS score of Patients
dataFINAL <- merge(dataFINAL, Metadata_NeuroTest, by = c("SampleID", "Timepoint"), all.x = TRUE)

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


#### MEAN values for the ANOVA & so ----------------
dataFINALmean <- dataFINAL %>%
    group_by(SampleID, Sex, Age, Category, Timepoint, DaysPS, Subpopulation, 
             Gene, GeneID, Gene_Group, Recovery, NHISS, mRS, Barthel, MoCA,
             `HADS_Anxiety`, `HADS_Depression`, NHISS_End, NHISS_Ratio) %>%
    summarise(mean_Z = mean(inverted_z_score), sd_Z = sd(inverted_z_score),
              mean_capZ = mean(capped_z_score), sd_capZ = sd(capped_z_score)
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
write_xlsx(FACSdata, "FACSdata.xlsx")
write.csv(FACSdata, "FACSdata_clean.csv", row.names = FALSE)
summary(dataFINAL)
summary(dataFINALmean)
