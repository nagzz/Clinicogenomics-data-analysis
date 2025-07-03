#Interaction among immunecells classification, Neoantigen expression levels and POLe mutations
library(dplyr)

# Create summary data with NA filtering
summary_data <- annotations_data1 %>%
    filter(!is.na(Macrophages_M1_Classification) & 
               !is.na(T_cells_CD8_Classification) & 
               !is.na(NeoExpressionLevel) &
               !is.na(Protein_change)) %>%
    count(Protein_change, Macrophages_M1_Classification, 
          T_cells_CD8_Classification, NeoExpressionLevel) %>%
    group_by(Protein_change) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    unite("Classification_Combo", T_cells_CD8_Classification, NeoExpressionLevel, sep = "_")

# Get top protein changes for readability
top_proteins <- summary_data %>%
    group_by(Protein_change) %>%
    summarise(total = sum(n)) %>%
    arrange(desc(total)) %>%
    slice_head(n = 15) %>%
    pull(Protein_change)

# Create the plot with font size 12 and black circle outlines
p <- ggplot(summary_data %>% filter(Protein_change %in% top_proteins), 
       aes(x = Classification_Combo, y = Protein_change, 
           fill = prop, size = n)) +
    geom_point(shape = 21, color = "black") +  # Changed from "white" to "black"
    facet_wrap(~paste("Macrophages M1:", Macrophages_M1_Classification), ncol = 2) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    scale_size_continuous(range = c(1, 8), 
                          breaks = c(1, 2, 3),
                          name = "Count") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    labs(title = "Protein Changes vs Classifications",
         x = "T_cells_CD8 & NeoExpression Level",
         y = "Protein Change",
         fill = "Proportion",
         size = "Count")

# Display the plot with specific dimensions
print(p)

# To save with specific dimensions:
# ggsave("protein_classifications_plot.png", plot = p, 
#        width = 14, height = 10, dpi = 300, units = "in")
