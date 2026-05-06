library(ggplot2)
library(dplyr)
library(ggrepel)

volcano_plot_linreg <- function(results_df,
                                timepoint = NULL,
                                subpopulation = NULL,
                                gene_group = NULL,
                                p_cutoff = 0.05,
                                coef_cutoff = 0,
                                label_top_n = 10,
                                effect_col = "Coefficient",
                                plot_title = "Volcano plot",
                                output_file = NULL) {
    
    df <- results_df %>%
        filter(!is.na(p.value), !is.na(.data[[effect_col]]))
    
    # Optional filters
    if (!is.null(timepoint)) {
        df <- df %>% filter(Timepoint %in% timepoint)
    }
    if (!is.null(subpopulation)) {
        df <- df %>% filter(Subpopulation %in% subpopulation)
    }
    if (!is.null(gene_group)) {
        df <- df %>% filter(Gene_Group %in% gene_group)
    }
    
    # Add volcano metrics
    df <- df %>%
        mutate(
            negLog10P = -log10(p.value),
            Significance = case_when(
                p.value < p_cutoff & .data[[effect_col]] >  coef_cutoff ~ "Up",
                p.value < p_cutoff & .data[[effect_col]] < -coef_cutoff ~ "Down",
                TRUE ~ "NS"
            )
        )
    
    # Select genes to label
    label_df <- df %>%
        filter(Significance != "NS") %>%
        arrange(p.value) %>%
        slice_head(n = label_top_n)
    
    p <- ggplot(df, aes(x = .data[[effect_col]], y = negLog10P)) +
        geom_point(aes(color = Significance), alpha = 0.8, size = 2) +
        scale_color_manual(
            values = c(
                "Up" = "#B2182B",
                "Down" = "#2166AC",
                "NS" = "grey70"
            )
        ) +
        geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
        geom_vline(xintercept = c(-coef_cutoff, coef_cutoff),
                   linetype = "dashed", color = "black") +
        ggrepel::geom_text_repel(
            data = label_df,
            aes(label = Gene),
            size = 4,
            max.overlaps = Inf
        ) +
        labs(
            title = plot_title,
            x = "Regression coefficient",
            y = expression(-logp.value),
            color = "Effect"
        ) +
        theme_bw(base_size = 14) +
        theme(
            panel.grid.minor = element_blank(),
            legend.position = "right"
        )
    
    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 8, height = 6, dpi = 300)
    }
    
    return(p)
}

subpops <- sort(unique(LinReg_NHISS_Ratio_capZ$Subpopulation))
timepoints <- sort(unique(LinReg_NHISS_Ratio_capZ$Timepoint))


for (sp in subpops) {
    for (tp in timepoints) {
        
        df_check <- LinReg_NHISS_Ratio_capZ %>%
            filter(
                Subpopulation == sp,
                Timepoint == tp,
                Subgroup == "All",
                !is.na(p.value),
                !is.na(Coefficient)
            )
        
        # Skip combinations with no data
        if (nrow(df_check) == 0) {
            message("Skipping: ", sp, " | ", tp, " (no data)")
            next
        }
        
        message("Plotting: ", sp, " | ", tp)
        
        volcano_plot_linreg(
            results_df    = LinReg_NHISS_Ratio_capZ,
            subpopulation = sp,
            timepoint     = tp,
            p_cutoff      = 0.05,
            coef_cutoff   = 0,
            label_top_n   = 10,
            plot_title    = paste(
                tools::toTitleCase(sp),
                "Monocytes (", tp, "): NHISS Ratio ~ Gene expression",
                sep = ""
            ),
            output_file   = paste0(
                "LinReg_Results_capZ/Volcano_",
                sp, "_", tp, "_NHISS_Ratio.png"
            )
        )
    }
}
