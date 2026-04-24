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
    colors <- c("blue", "white", "red")
    
    # Set the breakpoints for the gradient (-1.2 to 1.2)
    scale_color_gradientn(colors = colors,
                          limits = c(-1.2, 1.2),  # Set the limits from -1.2 to 1.2
                          guide = "colorbar",
                          na.value = "grey50")  # Color for NA values
}

ggplot(data_summary, aes(x = Sample_Combo, y = Gene_Combined)) +
    geom_point(aes(size = Dot_Size, color = mean_capZ)) +  # Use Dot_Size for scaling
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
    geom_point(aes(size = Dot_Size, color = mean_capZ)) +  # Use Dot_Size for scaling
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
    geom_point(aes(size = Dot_Size, color = mean_capZ)) +  # Use Dot_Size for scaling
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
    geom_point(aes(size = Dot_Size, color = mean_capZ)) +  # Use Dot_Size for scaling
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

#### PCA -------
dataFINALmean <- dataFINALmean %>%
    mutate(PCA_Sample = paste(Timepoint, Subpopulation, SampleID, sep = "_"))

# Reverse melt to create matrix for PCA
data_wide <- dcast(dataFINALmean, PCA_Sample ~ Gene, value.var = "mean_capZ")
rownames(data_wide) <- data_wide$PCA_Sample
data_wide <- data_wide[, -1]  # Remove the Sample column to leave only the matrix

data.meta <- t( data.frame(stringr::str_split( rownames( data_wide), "_")))
colnames(data.meta) <- c("Timepoint", "Subpopulation", "SampleID")
data.meta <- merge(data.meta , metadataP, by = "SampleID", all.x = TRUE)
#data.meta <- cbind( data.meta,  dataFINALmean[rownames(data.meta), "Age"])

pca_result <- prcomp(data_wide, scale. = TRUE)

plot_my_pca(pca_result, data.meta, "Timepoint", my_colors, legend_pos = "topleft")

subtype_col <- c("grey", "red", "orange", "green")
plot_my_pca(pca_result, data.meta, "Subpopulation", subtype_col, "topleft")
Category_col <- c("red", "orange", "green")
plot_my_pca(pca_result, data.meta, "Category", Category_col, "topleft")
plot_my_pca(pca_result, data.meta, "Recovery", c("red", "darkgreen","grey"), "topleft")
plot_my_pca(pca_result, data.meta, "Sex", c("#A67C00", "#1D04C2"), "topleft")

data.meta$Age <- factor(data.meta$Age, levels = sort(unique(data.meta$Age)), ordered = TRUE)
# Define Gradient_colour based on the levels of Age
Gradient_colour <- setNames(
    colorRampPalette(c("blue", "red"))(length(levels(data.meta$Age))), 
    levels(data.meta$Age)
)

# Run the function to plot PCA, mapping Age levels to colors
plot_my_pca(pca_result, data.meta, "Age", Gradient_colour[data.meta$Age])