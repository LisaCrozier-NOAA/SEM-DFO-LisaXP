var_lookup_NCC_AK<-read.csv("LisaXP/outputs_4/var_lookup_NCC_AK.SAR.csv",row.names=1);names(var_lookup_NCC_AK)

indicator_summary_lisa<-read.csv("LisaXP/outputs_4/DAG1_summarytable_dAIC3.csv",row.names = 1)
indicator_summary_lisa


importance_topvar<-read.csv("LisaXP/outputs_4/importance_topvar_daic3.csv",row.names=1);importance_topvar
head(importance_topvar)

# Define colors that actually match the methodology workflow visuals
cartoon_colors <- c(
  "Foundations"    = "#A3C1AD", # Sage Green (The background of Panel 1 & 3)
  "Data/Assembly"  = "#D9E4EC", # Pale Blue (The paper/background in Panel 2)
  "Functional"     = "#5D8AA8", # Steel Blue (The "Bottom-Up" bin in Panel 3)
  "Simulation"     = "#E67E22", # Burnt Orange (The "Salmon" bin and SEM gears)
  "Dimension/DFA"  = "#7F8C8D", # Slate Gray (The DFA cylinders in Panel 5)
  "Selection/AIC"  = "#F39C12", # Goldenrod (The AIC medal/ribbon in Panel 6)
  "Latent/AK"      = "#2E4053"  # Deep Navy (The deep ocean in Panel 7)
)


#Barchart=================
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggh4x) 

# 1. DAG1C only---------
      # Reshape and Refine Data
      df_long <- importance_topvar %>%
        mutate(
          Trophic = fct_reorder(Trophic, order),
          Region = fct_reorder(Region, order),
          plot_names = fct_rev(fct_reorder(Lisaname, order))
        ) %>%
        pivot_longer(
          cols = c(DAG1C_long, DAG1C_short),
      #    cols = c(DAG1A_long, DAG1B_long, DAG1A_short, DAG1B_short),
          names_to = c("variable", "type"),
          names_sep = "_",
          values_to = "importance"
        )

      # add Missing data labels
      missing_labels <- df_long %>%
        group_by(plot_names, type, Trophic, Region) %>%
        summarize(all_missing = all(is.na(importance)), .groups = 'drop') %>%
        filter(all_missing == TRUE)
      
      # 2. Create the Plot
      # scale_fill_gradient2(low = "#2E4053", 
      #                      mid = "white", 
      #                      high = "#D35400", 
                           
      importance_topvar_barchart<-
        ggplot(df_long, aes(x = plot_names, y = importance, fill = variable)) +
        geom_col(position = "dodge") +
      #  scale_fill_manual(values = c("DAG1A" = "#2E4053", "DAG1B" = "#D35400")) +
        scale_fill_manual(values = c("DAG1C" = "#2E4053")) +
        geom_text(data = missing_labels, 
                  aes(x = plot_names, y = 0, label = "Missing data"), 
                  inherit.aes = FALSE,
                  hjust = -0.1, 
                  color = "grey40", 
                  fontface = "italic",
                  size = 3) +
        coord_flip() + 
        facet_nested(Region + Trophic ~ type, 
                     scales = "free_y", 
                     space = "free_y",
                     nest_line = element_line(color = "black"),
                     resect = unit(2, "mm")) + 
        # Fix: Defined specific breaks for the vertical grid lines
        scale_y_continuous(
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          limits = c(0, 1.1), # Ensures space for labels, adjust to 1 if preferred
          expand = expansion(mult = c(0, 0.05))
        ) +
        theme_minimal() +
        labs(
          x = "Indicator", 
          y = "Importance",
          fill = NULL 
        ) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          # Keeps only the vertical lines at our specific breaks
          panel.grid.minor.x = element_blank(), 
          strip.text.y = element_text(angle = 0),
          panel.spacing.x = unit(2, "lines"), 
          panel.spacing.y = unit(1.5, "lines"),
          strip.background = element_blank()
        )
      
      importance_topvar_barchart
      
      ggsave( "LisaXP/outputs_4/DAG1C_importance_topvar_barchart.png",importance_topvar_barchart)

      
      
# 2. DAG1A and DAG1B -------------
       # Reshape and Refine Data
      df_long <- importance_topvar %>%
        mutate(
          Trophic = fct_reorder(Trophic, order),
          Region = fct_reorder(Region, order),
          plot_names = fct_rev(fct_reorder(Lisaname, order))
        ) %>%
        pivot_longer(
          #cols = c(DAG1C_long, DAG1C_short),
              cols = c(DAG1A_long, DAG1B_long, DAG1A_short, DAG1B_short),
          names_to = c("variable", "type"),
          names_sep = "_",
          values_to = "importance"
        )
      
      # add Missing data labels
      missing_labels <- df_long %>%
        group_by(plot_names, type, Trophic, Region) %>%
        summarize(all_missing = all(is.na(importance)), .groups = 'drop') %>%
        filter(all_missing == TRUE)
      
      # 2. Create the Plot
      # scale_fill_gradient2(low = "#2E4053", 
      #                      mid = "white", 
      #                      high = "#D35400", 
      
      importance_topvar_barchart<-
        ggplot(df_long, aes(x = plot_names, y = importance, fill = variable)) +
        geom_col(position = "dodge") +
          scale_fill_manual(values = c("DAG1A" = "#2E4053", "DAG1B" = "#D35400")) +
        #scale_fill_manual(values = c("DAG1C" = "#2E4053")) +
        geom_text(data = missing_labels, 
                  aes(x = plot_names, y = 0, label = "Missing data"), 
                  inherit.aes = FALSE,
                  hjust = -0.1, 
                  color = "grey40", 
                  fontface = "italic",
                  size = 3) +
        coord_flip() + 
        facet_nested(Region + Trophic ~ type, 
                     scales = "free_y", 
                     space = "free_y",
                     nest_line = element_line(color = "black"),
                     resect = unit(2, "mm")) + 
        # Fix: Defined specific breaks for the vertical grid lines
        scale_y_continuous(
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          limits = c(0, 1.1), # Ensures space for labels, adjust to 1 if preferred
          expand = expansion(mult = c(0, 0.05))
        ) +
        theme_minimal() +
        labs(
          x = "Indicator", 
          y = "Importance",
          fill = NULL 
        ) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          # Keeps only the vertical lines at our specific breaks
          panel.grid.minor.x = element_blank(), 
          strip.text.y = element_text(angle = 0),
          panel.spacing.x = unit(2, "lines"), 
          panel.spacing.y = unit(1.5, "lines"),
          strip.background = element_blank()
        )
      
      importance_topvar_barchart
      
      ggsave( "LisaXP/outputs_4/DAG1A_DAG1B_importance_topvar_barchart.png",importance_topvar_barchart)

      
      
# 3. All 3 DAGS on one plot--------      
      # Reshape and Refine Data
      df_long <- importance_topvar %>%
        mutate(
          Trophic = fct_reorder(Trophic, order),
          Region = fct_reorder(Region, order),
          plot_names = fct_rev(fct_reorder(Lisaname, order))
        ) %>%
        pivot_longer(
          cols = c(DAG1A_long, DAG1B_long, DAG1A_short, DAG1B_short,DAG1C_long, DAG1C_short),
          names_to = c("variable", "type"),
          names_sep = "_",
          values_to = "importance"
        )
      
      # add Missing data labels
      missing_labels <- df_long %>%
        group_by(plot_names, type, Trophic, Region) %>%
        summarize(all_missing = all(is.na(importance)), .groups = 'drop') %>%
        filter(all_missing == TRUE)
      
      # 2. Create the Plot
      # scale_fill_gradient2(low = "#2E4053", 
      #                      mid = "white", 
      #                      high = "#D35400", 
      
      importance_topvar_barchart<-
        ggplot(df_long, aes(x = plot_names, y = importance, fill = variable)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c("DAG1A" = "#2E4053", "DAG1B" = "#D35400", "DAG1C" = "#A3C1AD")) +
        geom_text(data = missing_labels, 
                  aes(x = plot_names, y = 0, label = "Missing data"), 
                  inherit.aes = FALSE,
                  hjust = -0.1, 
                  color = "grey40", 
                  fontface = "italic",
                  size = 3) +
        coord_flip() + 
        facet_nested(Region + Trophic ~ type, 
                     scales = "free_y", 
                     space = "free_y",
                     nest_line = element_line(color = "black"),
                     resect = unit(2, "mm")) + 
        # Fix: Defined specific breaks for the vertical grid lines
        scale_y_continuous(
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          limits = c(0, 1.1), # Ensures space for labels, adjust to 1 if preferred
          expand = expansion(mult = c(0, 0.05))
        ) +
        theme_minimal() +
        labs(
          x = "Indicator", 
          y = "Importance",
          fill = NULL 
        ) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          # Keeps only the vertical lines at our specific breaks
          panel.grid.minor.x = element_blank(), 
          strip.text.y = element_text(angle = 0),
          panel.spacing.x = unit(2, "lines"), 
          panel.spacing.y = unit(1.5, "lines"),
          strip.background = element_blank()
        )
      
      importance_topvar_barchart
      
      ggsave( "LisaXP/outputs_4/DAG1A_DAG1B_DAG1C_importance_topvar_barchart.png",importance_topvar_barchart)
      
      
      