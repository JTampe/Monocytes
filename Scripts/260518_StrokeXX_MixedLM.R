library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggpubr)

automate_mixedlm_extraction_Category <- function(
        results_folder,
        plots_folder,
        results_name,
        plot_save,
        dataset,
        loop_vars,
        timepoint_col = "Timepoint",
        subject_col = "SampleID",
        color_var
) {
    
    if (!dir.exists(results_folder)) dir.create(results_folder)
    if (!dir.exists(plots_folder)) dir.create(plots_folder)
    
    results <- data.frame(
        Category = character(),
        Variable = character(),
        Global_P_Value = numeric(),
        TP1_P_Value = numeric(),
        TP2_P_Value = numeric(),
        TP3_P_Value = numeric(),
        TP4_P_Value = numeric(),
        stringsAsFactors = FALSE
    )
    
    for (var_name in loop_vars) {
        
        df_test <- dataset %>%
            select(
                value = all_of(var_name),
                Timepoint = all_of(timepoint_col),
                Subject = all_of(subject_col),
                Category = all_of(color_var)
            ) %>%
            filter(!is.na(value), !is.na(Timepoint), !is.na(Subject), !is.na(Category))
        
        df_test$Timepoint <- factor(
            df_test$Timepoint,
            levels = c("TP0", "TP1", "TP2", "TP3", "TP4")
        )
        
        if (n_distinct(df_test$Timepoint) < 2) next
        
        for (cat in unique(df_test$Category)) {
            
            df_cat <- df_test %>% filter(Category == cat)
            
            # Need both controls and at least one longitudinal timepoint
            if (!("TP0" %in% df_cat$Timepoint) || nrow(df_cat) < 5) next
            
            # Mixed‑effects model
            model <- tryCatch(
                lmer(value ~ Timepoint + (1 | Subject), data = df_cat),
                error = function(e) NULL
            )
            
            if (is.null(model)) next
            
            # Global p‑value
            anova_tab <- anova(model)
            global_p <- anova_tab["Timepoint", "Pr(>F)"]
            
            # Post‑hoc vs control
            em <- emmeans(model, ~ Timepoint)
            ctr <- contrast(
                em,
                method = "trt.vs.ctrl",
                ref = "TP0",
                adjust = "holm"
            )
            
            ctr_tab <- as.data.frame(ctr)
            
            get_p <- function(tp) {
                row <- ctr_tab %>% filter(grepl(tp, contrast))
                if (nrow(row) == 0) NA else row$p.value
            }
            
            results <- rbind(
                results,
                data.frame(
                    Category = cat,
                    Variable = var_name,
                    Global_P_Value = global_p,
                    TP1_P_Value = get_p("TP1"),
                    TP2_P_Value = get_p("TP2"),
                    TP3_P_Value = get_p("TP3"),
                    TP4_P_Value = get_p("TP4")
                )
            )
            
            # ======================
            # Optional plotting
            # ======================
            if (plot_save == "y") {
                
                my_colors <- c(
                    "Female" = "#1D04C2",
                    "Male" = "#A67C00",
                    "MINOR" = "#A67C00",
                    "MODERATE" = "#700606",
                    "Bad" = "#700606",
                    "Good" = "#A67C00"
                )
                
                used_colors <- my_colors[names(my_colors) %in% df_cat$Category]
                
                box_plot <- ggboxplot(
                    df_cat,
                    x = "Timepoint",
                    y = "value",
                    color = "Category",
                    add = "jitter"
                ) +
                    scale_x_discrete(labels = c(
                        "TP0" = "Control",
                        "TP1" = "24 hours",
                        "TP2" = "3–5 days",
                        "TP3" = "1 month",
                        "TP4" = "3 months"
                    )) +
                    scale_color_manual(values = used_colors) +
                    ylab(paste(var_name, "percentage")) +
                    ggtitle(paste(var_name, "-", cat))
                
                ggsave(
                    filename = file.path(
                        plots_folder,
                        paste0(var_name, "_", cat, "_MixedLM.png")
                    ),
                    plot = box_plot,
                    width = 7.5,
                    height = 5
                )
            }
        }
    }
    
    write.csv(
        results,
        file.path(results_folder, paste0(results_name, "_results.csv")),
        row.names = FALSE
    )
    
    return(results)
}


MixedLM_FACS <- automate_mixedlm_extraction_Category(
    results_folder = output_location,
    plots_folder   = "MixedLM_FACS_Plots",
    results_name   = "MixedLM_FACS",
    plot_save      = "y",
    dataset        = ANOVA_FACSdata,
    loop_vars      = colnames(ANOVA_FACSdata)[9:12],
    timepoint_col  = "Timepoint",
    subject_col    = "SampleID",
    color_var      = "Sex"   # or severity / group
)
