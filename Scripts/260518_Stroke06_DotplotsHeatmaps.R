# *Dotplot of predictor Genes* -----------------------------------------------------------------------------
folder <- "Dotplots_Gene_Predictiors"
createFolder(folder)

# lets just focus on all!
Predictors  <- LinReg_NIHSS_End_sigGenes %>% filter(Subgroup == "All")
Predictors  <- Predictors %>% filter(Timepoint !=  "TP0")
Predictors  <- Predictors %>% filter(Timepoint !=  "TP4")
Predictors  <- Predictors %>% filter(Timepoint !=  "TP3")

# Relabel the timepoints: PS = post-stroke
Predictors$Timepoint[grep("TP1",Predictors$Timepoint)]  <- "24 hours PS"
Predictors$Timepoint[grep("TP2",Predictors$Timepoint)]  <- "3-5 days PS"

# Create a new column for Timepoint and Subpopulation combination
Predictors <- Predictors %>%
    mutate(Sample_Combo2 = paste(Timepoint, Subpopulation, sep = "_"))

Predictors <- Predictors %>%
    mutate(Sample_Combo = paste(Subpopulation, Timepoint, sep = "_"))

# Create a combined column for Gene and GeneID for better sorting
Predictors <- Predictors %>%
    mutate(Gene_Combined = paste(GeneID, Gene, sep = "_"))  # Combine GeneID and Gene

# Create a linear gradient function for the original plots
linear_gradient <- function() {
    # Define the colors for the gradient
    colors <- c("blue", "white", "red")
    
    # Set the breakpoints for the gradient (-1.2 to 1.2)
    scale_color_gradientn(colors = colors,
                          limits = c(-0.79, 0.79),  # Set the limits from -1.2 to 1.2
                          guide = "colorbar",
                          na.value = "grey50")  # Color for NA values
}

# Ensure Sample_Combo is ordered alphabetically
Predictors <- Predictors %>%
    mutate(Sample_Combo = factor(Sample_Combo, levels = sort(unique(Sample_Combo))))

ggplot(Predictors, aes(x = Gene_Combined, y = Sample_Combo)) +
    geom_point(aes(size = -log10(p.value), color = Estimate)) +  # Size uses -log10(p.value) for better scaling
    scale_size(range = c(1, 7), name = "p-value\n(-log10)") +  # Customize dot size range
    linear_gradient() +  # Apply the custom color gradient
    theme_minimal() +
    labs(
        title = "Dot Plot of Predictors Data",
        color = "Estimate",
        size = "Significance"
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "lightgrey", size = 0.25),
        panel.grid.minor = element_line(color = "white", size = 0.25)
    ) +
    scale_y_discrete(limits = rev(levels(Predictors$Sample_Combo)))  # Reverse y-axis order

# Save the plot
ggsave(filename = "Dotplots_Gene_Predictiors/Dotplot_Predictors.png", width = 10, height = 8)

# *Heatmaps* ----------------------------------
new_folder <- "Pearsson_Heatmaps"
createFolder(new_folder)

# Create a linear gradient function for the original plots
linear_gradient <- function() {
    # Define the colors for the gradient
    colors <- c("blue", "white", "red")
    
    # Set the breakpoints for the gradient (-1.2 to 1.2)
    scale_color_gradientn(colors = colors,
                          limits = c(-0.79, 0.79),  # Set the limits from -1.2 to 1.2
                          guide = "colorbar",
                          na.value = "grey50")  # Color for NA values
}

#### NIHSS_Ratio (Estimate) -----------------------
Predictors_Ratio <- LinReg_NIHSS_Ratio_capZ 
# Relabel the timepoints: 
Predictors_Ratio$Timepoint[grep("TP1",Predictors_Ratio$Timepoint)]  <- "24 hours"
Predictors_Ratio$Timepoint[grep("TP2",Predictors_Ratio$Timepoint)]  <- "3-5 days"

#Create a new column for Timepoint and Subpopulation combination
Predictors_Ratio  <- Predictors_Ratio  %>%
    mutate(Sample_Combo2 = paste(Timepoint, Subpopulation, sep = "_")) %>%
    mutate(Sample_Combo = paste(Subpopulation, Timepoint, sep = "_"))  %>%
    mutate(Gene_Combined = paste(GeneID, Gene, sep = "_"))
# Ensure Sample_Combo is ordered alphabetically
Predictors_Ratio  <- Predictors_Ratio  %>%
    mutate(Sample_Combo = factor(Sample_Combo, levels = sort(unique(Sample_Combo)))) %>%
    mutate(Sample_Combo2 = factor(Sample_Combo2, levels = sort(unique(Sample_Combo2)))) %>%
    mutate(Gene_Combined = paste(GeneID, Gene, sep = "_"))

Predictors_Ratio$Significance <- ifelse(Predictors_Ratio $p.value < 0.05, "*", "")

#### Heatmap - all genes
name_plot <- "Correlation of acute (24 hours & 3-5 days PS) gene expr. with relative change in NIHSS (24h vs 3 months PS)"
ggplot(Predictors_Ratio , aes(Gene_Combined, Sample_Combo2, fill= Estimate)) + 
    geom_tile(color = "white") +  # Add white grid lines
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-0.7, 0.8)) +
    geom_text(aes(label = Significance), color = "white", size = 10, 
              hjust = 0.5, vjust = 0.5) +  # Center the text in the squares
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, y = "Timepoint & Suptype", x = "Gene", fill = "Estimate")+
    scale_y_discrete(limits = rev(levels(Predictors_Ratio $Sample_Combo2))) 
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 15, height = 4.5, dpi = 300)

#### NIHSS_Ratio - sig. Genes: -----------------
Predictors_Ratio  <- merge(Predictors_Ratio , Metadata_GeneID, by = "Gene")

Predictors_Ratio  <- Predictors_Ratio  %>%
    mutate(Gene_Combined2 = paste(GeneID_Ratio, Gene, sep = "_"))  # Combine GeneID and Gene

# Filter rows where p-value is <= 0.5
significant_genes <- Predictors_Ratio  %>% filter(Predictors_Ratio$p.value <0.05)
significant_genes <- unique(significant_genes$Gene)

# adjust data frame accordingly
Predictors_Ratio  <- Predictors_Ratio  %>%
    filter(Gene %in% significant_genes)

#### Heatmap - sig. Genes
name_plot <- "Correlation of acute (24 hours & 3-5 days PS) gene expr. with relative change in NIHSS (24h vs 3 months PS) - only sig. Genes"
ggplot(Predictors_Ratio , aes(Gene_Combined2, Sample_Combo2, fill= Estimate)) + 
    geom_tile(color = "white") +  # Add white grid lines
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-0.7, 0.8)) +
    geom_text(aes(label = Significance), color = "white", size = 12, 
              hjust = 0.5, vjust = 0.5) +  # Center the text in the squares
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, y = "Timepoint & Suptype", x = "Gene", fill = "Estimate")+
    scale_y_discrete(limits = rev(levels(Predictors_Ratio $Sample_Combo2))) 
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 9, height = 5, dpi = 300)

#### NIHSS_Ratio - TP1/TP2 sig. Genes: -----------------
Predictors_Ratio_TP1  <- Predictors_Ratio %>% filter(Timepoint ==  "24 hours")
Predictors_Ratio_TP2  <- Predictors_Ratio %>% filter(Timepoint ==  "3-5 days")
# Filter rows where p-value is <= 0.5
significant_genes_TP1 <- Predictors_Ratio_TP1  %>% filter(Predictors_Ratio_TP1$p.value <0.05)
significant_genes_TP1 <- unique(significant_genes_TP1$Gene)
significant_genes_TP2 <- Predictors_Ratio_TP2  %>% filter(Predictors_Ratio_TP2$p.value <0.05)
significant_genes_TP2 <- unique(significant_genes_TP2$Gene)

# adjust data frame accordingly
Predictors_Ratio_TP1 <- Predictors_Ratio_TP1 %>% filter(Gene %in% significant_genes_TP1)
Predictors_Ratio_TP2 <- Predictors_Ratio_TP2 %>% filter(Gene %in% significant_genes_TP2)

#### Heatmap - TP1 sig. Genes
name_plot <- "Correlation of acute (24 hours) gene expr. with relative change in NIHSS (24h vs 3 months PS) - only sig. Genes"
ggplot(Predictors_Ratio_TP1, aes(Gene_Combined2, Sample_Combo2, fill= Estimate)) + 
    geom_tile(color = "white") +  # Add white grid lines
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-0.65, 0.71)) +
    geom_text(aes(label = Significance), color = "white", size = 12, 
              hjust = 0.5, vjust = 0.8) +  # Center the text in the squares
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, y = "Suptype", x = "Gene", fill = "Estimate")+
    scale_y_discrete(limits = rev(levels(Predictors_Ratio_TP1$Sample_Combo2))) 
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 6, height = 5, dpi = 300)

#### Heatmap - TP2 sig. Genes
name_plot <- "Correlation of sub-acute (3-5 days) gene expr. with relative change in NIHSS (24h vs 3 months PS) - only sig. Genes"
ggplot(Predictors_Ratio_TP2, aes(Gene_Combined2, Sample_Combo2, fill= Estimate)) + 
    geom_tile(color = "white") +  # Add white grid lines
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-0.68, 0.8)) +
    geom_text(aes(label = Significance), color = "white", size = 12, 
              hjust = 0.5, vjust = 0.8) +  # Center the text in the squares
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, y = "Suptype", x = "Gene", fill = "Estimate")+
    scale_y_discrete(limits = rev(levels(Predictors_Ratio_TP2$Sample_Combo2))) 
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 6, height = 5, dpi = 300)
#### NIHSS (Estimate) ---------------
LinReg_NIHSS_capZ <- LinReg_NIHSS_capZ %>%
    mutate(Timepoint = factor(Timepoint, levels = sort(unique(Timepoint))))

LinReg_NIHSS_capZ <- LinReg_NIHSS_capZ %>%
    mutate(Gene_Ordered = paste(GeneID, Gene, sep = "_"))
LinReg_NIHSS_capZ <- LinReg_NIHSS_capZ %>%
    mutate(Groups = paste(Subpopulation, Timepoint, sep = "_"))

# Estimate: - 0.99 to 0.92
name_plot <- paste("Heatmap of NIHSS Correlation Estimates")
# Add significance levels
LinReg_NIHSS_capZ$Significance <- ifelse(LinReg_NIHSS_capZ$p.value < 0.05, "*", "")

# Create the heatmap with significance levels
ggplot(LinReg_NIHSS_capZ, aes(x = Gene_Ordered, y = Groups, fill = Estimate)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "purple", 
                         midpoint = 0, limits = c(-0.72, 0.78)) +
    geom_text(aes(label = Significance), color = "black", size = 3) + # Add significance levels with white text
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, x = "Gene", y = "Timepoint & Subpopulation", fill = "Estimate")
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 10, height = 5, dpi = 300)

name_plot <- paste(name_plot," (transposed)")
ggplot(LinReg_NIHSS_capZ, aes(Groups, Gene_Ordered, fill= Estimate)) + 
    geom_tile() +
    scale_fill_gradient2(low = "yellow", mid = "white", high = "purple", 
                         midpoint = 0, limits = c(-0.72, 0.78)) +
    geom_text(aes(label = Significance), color = "black", size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = name_plot, x = "Timepoint & Subpopulation", y = "Timepoint", fill = "Estimate")
ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
       plot = last_plot(), width = 5, height = 10, dpi = 300)

#### NIHSS Split by Subpopulation (Estimate) ---------------
timepoint_labels <- c("TP0" = "Healthy", "TP1" = "24 hours post-stroke", "TP2" = "3-5 days post-stroke", 
                      "TP3" = "1 month post-stroke", "TP4" = "3 months post-stroke")

for (sup in unique(LinReg_NIHSS_capZ$Subpopulation)) {
    df <- filter(LinReg_NIHSS_capZ, Subpopulation == sup)
    sigGenes <- filter(LinReg_NIHSS_capZ, p.value <= 0.05) 
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of NIHSS Correlation Estimates - ", sup, " Monocytes") 
    ggplot(df, aes(Gene_Ordered, Timepoint, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "yellow", mid = "white", high = "purple", 
                             midpoint = 0, limits = c(-0.72, 0.78)) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Timepoint", fill = "Coefficient")+
        scale_y_discrete(labels = timepoint_labels)+# Apply custom y-axis labels
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}
for (sup in unique(LinReg_NIHSS_capZ$Subpopulation)) {
    df <- filter(LinReg_NIHSS_capZ, Subpopulation == sup)
    # sigGenes <- filter(LinReg_NIHSS_capZ, p.value <= 0.05) if you want to compare all of them
    sigGenes <- filter(df, p.value <= 0.05)
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of NIHSS Correlation Estimates - ", sup, " Monocytes (only sigGenes)") 
    ggplot(df, aes(Gene_Ordered, Timepoint, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "yellow", mid = "white", high = "purple", 
                             midpoint = 0, limits = c(-0.72, 0.78)) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Timepoint", fill = "Coefficient")+
        scale_y_discrete(labels = timepoint_labels)+# Apply custom y-axis labels
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}

#### NIHSS Split by Timepoint (Estimate) ---------------
for (tp in unique(LinReg_NIHSS_capZ$Timepoint)) {
    df <- filter(LinReg_NIHSS_capZ, Timepoint == tp)
    sigGenes <- filter(LinReg_NIHSS_capZ, p.value <= 0.05) 
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of NIHSS Correlation Estimates - ", tp) 
    ggplot(df, aes(Gene_Ordered, Subpopulation, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "yellow", mid = "white", high = "purple", 
                             midpoint = 0, limits = c(-0.72, 0.78)) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Monocyte Suptype", fill = "Coefficient")+
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}
for (tp in unique(LinReg_NIHSS_capZ$Timepoint)) {
    df <- filter(LinReg_NIHSS_capZ, Timepoint == tp)
    # sigGenes <- filter(LinReg_NIHSS_capZ, p.value <= 0.05) if you want to compare all of them
    sigGenes <- filter(df, p.value <= 0.05)
    sigGenes <- unique(sigGenes$Gene)
    df <- df %>%  filter(Gene %in% sigGenes)
    df <- df %>%  filter(!Gene %in% c("ACTB", "B2M", "GAPDH"))
    
    name_plot <- paste("Heatmap of NIHSS Correlation Estimates - ", tp, " (only sigGenes)") 
    ggplot(df, aes(Gene_Ordered, Subpopulation, fill= Estimate)) + 
        geom_tile() +
        scale_fill_gradient2(low = "yellow", mid = "white", high = "purple", 
                             midpoint = 0, limits = c(-0.72, 0.78)) +
        geom_text(aes(label = Significance), color = "black", size = 6) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = name_plot, x = "Gene", y = "Monocyte Suptype", fill = "Coefficient")+
        scale_x_discrete(labels = function(x) str_sub(x, start = 4)) + # Remove first three characters
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14), # Adjust font size
            axis.text.y = element_text(size = 14),                        # Adjust font size
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold")
        )
    ggsave(filename = paste(new_folder, "/", name_plot, ".png", sep = ""),
           plot = last_plot(), width = 10, height = 5, dpi = 300)
}

sigGenes <- filter(LinReg_NIHSS_capZ, p.value <= 0.05) 
unique(sigGenes$Gene)
