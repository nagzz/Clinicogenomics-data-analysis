# neoantigen count dotplot
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)

df <- read_delim("SLIDlist_neoantigenCount.txt", delim = "\t", col_names = FALSE)

# Extract relevant columns and rename for clarity
plot_data <- df %>%
    select(neoantigen_count = X2, group = X4, type = X6) %>%
    # Filter out any Unknown/Not Applicable values
    filter(group %in% c("I", "II", "III", "IV"), 
           type %in% c("mut", "wild")) %>%
    mutate(type = factor(type, levels = c("wild", "mut")))

# Get total sample size
n_total <- nrow(plot_data)

p <- ggplot(plot_data, aes(x = type, y = neoantigen_count)) +
    geom_dotplot(
        binaxis = 'y', 
        stackdir = 'center', 
        stackratio = 1, 
        dotsize = 0.5,
        binwidth = 15,
        aes(fill = type),
        colour = "black"
    ) +
    scale_fill_manual(values = c("wild" = "#4DAFAC", "mut" = "#E8756C")) +
    scale_y_continuous(breaks = seq(0, 400, by = 50)) +
    facet_wrap(~ group, ncol = 4) +
    labs(
        title = "Neoantigen Count by PolE Type and Group",
        x = "Pt - PolE Type",
        y = "Neoantigen Count",
        fill = "Type",
        caption = paste("n =", n_total)
    ) +
    theme_bw() +
    theme(
        panel.grid.major = element_line(color = "#EEEEEE"),
        panel.grid.minor = element_line(color = "#EEEEEE"),
        strip.background = element_rect(fill = "#DDDDDD"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        # Position caption under the legend
        plot.caption = element_text(
            hjust = 1.06,  # Right align
            vjust = 50,
            size = 16,  # Match other text size
            margin = margin(t = 0, b = 10)
        )
    )
ggsave("neoantigen_count_dotplot.png", plot = p, width = 18, height = 6, dpi = 300)
