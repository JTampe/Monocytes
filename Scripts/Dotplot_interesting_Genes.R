# *Dotplot of interesting genes genes Genes* -----------------------------------------------------------------------------
# rdefine the genes from control paper
interesting_Genes <- c("CD36", "TNFA", "ANXA1", "TLR8", "TGFB1", "CD91" )

dataFINALmean_sigGENES<- dataFINALmean %>%
    filter(Gene %in% interesting_Genes & !Gene %in% c("ACTB", "B2M", "GAPDH"))

unique(dataFINALmean_sigGENES$Gene_Combined)

# Using the function for data_summary and data_summary2
data_summary <- calculate_summary(dataFINALmean_sigGENES, "Sample_Combo", "Gene_Combined")
data_summary2 <- calculate_summary(dataFINALmean_sigGENES, "Sample_Combo2", "Gene_Combined")

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
ggsave(filename = "Dotplots_interesting_Genes_summary/Exp_by_TP.png")

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
ggsave(filename = "Dotplots_interesting_Genes_summary/Exp_by_Suptype.png")

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
ggsave(filename = "Dotplots_interesting_Genes_summary/Exp_by_TP_grey.png")

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
ggsave(filename = "Dotplots_interesting_Genes_summary/Exp_by_Suptype_grey.png")