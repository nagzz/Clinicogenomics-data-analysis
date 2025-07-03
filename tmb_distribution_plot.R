# Load required libraries
library(ggplot2)
library(dplyr)

# Reorder the Type factor to put Wild first
tmb$Type <- factor(tmb$Type, levels = c("Wild", "Mut"))

# Calculate sample sizes for each group
sample_sizes <- tmb %>%
    group_by(Type) %>%
    summarise(n = n(), .groups = 'drop')

# Create labels with sample sizes (Wild will be first due to factor ordering)
labels_with_n <- c(
    "Wild" = paste0("pole_wild_type\n(n=", sample_sizes$n[sample_sizes$Type == "Wild"], ")"),
    "Mut" = paste0("pole_mutated\n(n=", sample_sizes$n[sample_sizes$Type == "Mut"], ")")
)

# Create scatter plot with box plot overlay
ggplot(tmb, aes(x = Type, y = TMB)) +
    # Add jittered points
    geom_jitter(width = 0.2, alpha = 0.7, color = "#3498DB", size = 2) +
    # Add box plot overlay
    geom_boxplot(alpha = 0.3, width = 0.5, outlier.shape = NA, 
                 fill = "transparent", color = "black") +
    # Customize labels
    labs(title = "TMB Distribution by Type",
         x = "Group",
         y = "TumorMutationalBurden") +
    # Apply minimal theme
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    # Adjust x-axis labels to match your plot
    scale_x_discrete(labels = labels_with_n)

# Alternative version with more customization to match your exact style
ggplot(tmb, aes(x = Type, y = TMB)) +
    geom_jitter(width = 0.1, alpha = 0.8, color = "#2E86AB", size = 1.2) +
    geom_boxplot(alpha = 0, width = 0.25, outlier.shape = NA, 
                 color = "black", linewidth = 0.6) +
    labs(x = "Group", y = "TMB") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor.y = element_line(color = "grey95", linewidth = 0.2)
    ) +
    scale_x_discrete(labels = labels_with_n) +
    scale_y_continuous(breaks = seq(0, 500, 50))

# Save the plot with today's date
today_date <- format(Sys.Date(), "%Y%m%d")
ggsave(paste0("tmb_distribution_plot_", today_date, ".svg"), width = 4, height = 4, dpi = 300)

# Alternative save options with date:
# ggsave(paste0("tmb_distribution_plot_", today_date, ".pdf"), width = 6, height = 5)
# ggsave(paste0("tmb_distribution_plot_", today_date, ".jpg"), width = 6, height = 5, dpi = 300)
# ggsave(paste0("tmb_distribution_plot_", today_date, ".tiff"), width = 6, height = 5, dpi = 300)
ggsave(paste0("tmb_distribution_plot_", today_date, ".png"), width = 5, height = 4, dpi = 300)
