library(ggplot2)
library(dplyr)

# Handle NA values in NeoExpressionLevel by converting to "Low"
annotation_with_tmb_neo_clean <- annotation_with_tmb_neo %>%
    mutate(NeoExpressionLevel = ifelse(is.na(NeoExpressionLevel), "Low", as.character(NeoExpressionLevel))) %>%
    # Reorder NeoExpressionLevel to High, Medium, Low
    mutate(NeoExpressionLevel = factor(NeoExpressionLevel, levels = c("High", "Medium", "Low"))) %>%
    # Reorder Type factor to put Wild first
    mutate(Type = factor(Type, levels = c("Wild", "Mut")))

# Calculate sample sizes for each group and neoexpression level
sample_sizes <- annotation_with_tmb_neo_clean %>%
    group_by(Type, NeoExpressionLevel) %>%
    summarise(n = n(), .groups = 'drop')

# Print sample sizes to verify
print("Sample sizes by group:")
print(sample_sizes)

# Create the compact plot
p_compact <- annotation_with_tmb_neo_clean %>%
    ggplot(aes(x = Type, y = TMB, color = NeoExpressionLevel)) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.8, size = 1.2) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.2, color = "black", linewidth = 0.3) +
    facet_wrap(~ NeoExpressionLevel, nrow = 1, strip.position = "top") +
    scale_color_manual(values = c("High" = "#E31A1C", "Medium" = "#FF7F00", "Low" = "#1F78B4")) +
    scale_y_continuous(trans = "log10", 
                       breaks = c(1, 10, 100, 1000),
                       labels = c("1", "10", "100", "1000")) +
    scale_x_discrete(labels = function(x) {
        sapply(x, function(type) {
            # Get sample size for this type in current facet
            n_val <- sample_sizes$n[sample_sizes$Type == type & 
                                        sample_sizes$NeoExpressionLevel == "High"][1]  # Will be replaced per facet
            base_label <- ifelse(type == "Wild", "pole_wild_type", "pole_mutated")
            return(paste0(base_label, "\n(n=", n_val, ")"))
        })
    }) +
    labs(
        x = "",
        y = "TMB (log10)"
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 8, face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.2),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.3),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        panel.spacing.x = unit(0.05, "cm"),
        panel.spacing.y = unit(0.05, "cm"),
        plot.margin = margin(3, 3, 3, 3, "pt")
    )

# Get the overall y-axis range for consistent scaling
y_min <- min(annotation_with_tmb_neo_clean$TMB, na.rm = TRUE)
y_max <- max(annotation_with_tmb_neo_clean$TMB, na.rm = TRUE)

# Create separate plots for each facet with correct sample sizes and consistent scaling
create_facet_plot <- function(neo_level) {
    facet_data <- annotation_with_tmb_neo_clean %>%
        filter(NeoExpressionLevel == neo_level)
    
    facet_sample_sizes <- sample_sizes %>%
        filter(NeoExpressionLevel == neo_level)
    
    ggplot(facet_data, aes(x = Type, y = TMB, color = NeoExpressionLevel)) +
        geom_jitter(width = 0.1, height = 0, alpha = 0.8, size = 1.2) +
        geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.2, color = "black", linewidth = 0.3) +
        scale_color_manual(values = c("High" = "#E31A1C", "Medium" = "#FF7F00", "Low" = "#1F78B4")) +
        scale_y_continuous(limits = c(y_min, y_max), 
                           breaks = seq(0, ceiling(y_max/50)*50, by = 50)) +
        scale_x_discrete(labels = function(x) {
            sapply(x, function(type) {
                n_val <- facet_sample_sizes$n[facet_sample_sizes$Type == type]
                if(length(n_val) == 0) n_val <- 0
                base_label <- ifelse(type == "Wild", "pole_wild_type", "pole_mutated")
                return(paste0(base_label, "\n(n=", n_val, ")"))
            })
        }) +
        labs(title = neo_level, x = "", y = if(neo_level == "High") "TMB" else "") +
        theme_classic() +
        theme(
            axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7, color = "black"),
            axis.title.y = element_text(size = 8, color = "black"),
            plot.title = element_text(size = 8, face = "bold", color = "black", hjust = 0.5),
            legend.position = "none",
            panel.grid.major.y = element_line(color = "grey95", linewidth = 0.2),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.3),
            axis.ticks = element_line(color = "black", linewidth = 0.2),
            plot.margin = margin(3, 3, 3, 3, "pt")
        )
}

# Create individual plots
p_high <- create_facet_plot("High")
p_medium <- create_facet_plot("Medium")
p_low <- create_facet_plot("Low")

# Combine plots using patchwork or gridExtra
library(gridExtra)
p_final <- grid.arrange(p_high, p_medium, p_low, ncol = 3)

# Display the plot
print(p_final)

# Save the plot
today_date <- format(Sys.Date(), "%Y%m%d")

# For gridExtra output, save differently
png(filename = paste0("TMB_neoantigen_compact_", today_date, ".png"),
    width = 4.2, height = 2.4, units = "in", res = 300)
grid.arrange(p_high, p_medium, p_low, ncol = 3)
dev.off()

pdf(file = paste0("TMB_neoantigen_compact_", today_date, ".pdf"),
    width = 4.2, height = 2.4)
grid.arrange(p_high, p_medium, p_low, ncol = 3)
dev.off()
