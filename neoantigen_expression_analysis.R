# Neoantigen Expression Analysis
# 
# This script performs two main tasks:
# 1. Categorizing samples into high, medium, and low neoantigen expression levels
# 2. Generating a heatmap visualization of the neoantigen expression data

# Load required libraries
library(readr)      # For reading delimited files
library(dplyr)      # For data manipulation
library(pheatmap)   # For heatmap generation
library(ggplot2)    # For saving plots
library(grid)       # For grid graphics
library(gridExtra)  # For grid operations

#==============================================================================
# PART 1: Categorize samples by neoantigen expression levels
#==============================================================================

# Read annotation file
annotations <- read_delim("SLID_type.txt", delim = "\t", col_names = FALSE)

# Convert stages to standard format
annotations$X2 <- ifelse(annotations$X2 %in% c("I", "II", "III", "IV"), 
                          annotations$X2, 
                          "Uk")

# Create annotation dataframe
annotation_col <- data.frame(
    Type = factor(annotations$X3),    # Type from third column
    Stage = factor(annotations$X2, levels = c("I", "II", "III", "IV", "Uk")),   # Stage with ordered levels
    row.names = annotations$X1
)

# Calculate mean TPM for each sample
sample_means <- colMeans(final_mergedTPM[, 3:ncol(final_mergedTPM)], na.rm = TRUE)
names(sample_means) <- colnames(final_mergedTPM)[3:ncol(final_mergedTPM)]

# Determine thresholds for expression categories using quantiles
quantiles <- quantile(sample_means, probs = c(0.33, 0.66))
low_threshold <- quantiles[1]
high_threshold <- quantiles[2]

# Classify samples based on mean expression
expression_levels <- cut(sample_means, 
                         breaks = c(-Inf, low_threshold, high_threshold, Inf),
                         labels = c("Low", "Medium", "High"))

# Create dataframe with expression levels
sample_expression_df <- data.frame(
  Sample_ID = names(sample_means),
  NeoMeanExpression = sample_means,
  ExpressionLevel = expression_levels
)

# Add expression levels to annotation dataframe
annotation_with_neoexpression <- annotation_col
annotation_with_neoexpression$NeoExpressionLevel <- NA

for (sample_id in rownames(annotation_col)) {
  if (sample_id %in% names(sample_means)) {
    expression <- as.character(expression_levels[names(sample_means) == sample_id])
    annotation_with_neoexpression[sample_id, "NeoExpressionLevel"] <- expression
  }
}

# Remove unknown stage samples
annotation_with_neoexpression <- annotation_with_neoexpression[annotation_with_neoexpression$Stage != "Uk", ]
# Remove rows with NA in NeoExpressionLevel
annotation_clean <- annotation_with_neoexpression[!is.na(annotation_with_neoexpression$NeoExpressionLevel), ]

#==============================================================================
# PART 2: Generate neoantigen expression heatmap
#==============================================================================

# Select data from the third column onwards 
data_subset <- final_mergedTPM[, 3:ncol(final_mergedTPM)]

# Find column IDs with at least 2 gene expression values >= 2
valid_cols <- colSums(data_subset >= 2) >= 2
filtered_data <- data_subset[, valid_cols]

# Log transform filtered data
log2_data <- log2(filtered_data + 1)

# Create matrix and assign row names
matrix_data <- as.matrix(log2_data)
rownames(matrix_data) <- final_mergedTPM$Gene_ID

# z-score normalization on rows of matrix_data
scaled_data <- t(scale(t(matrix_data)))

# Calculate variance for selecting top genes
gene_vars <- apply(scaled_data, 1, var)
top <- names(sort(gene_vars, decreasing = TRUE))

# Calculate column consistency
col_cors <- cor(scaled_data[top,])
col_consistency <- colMeans(abs(col_cors))
col_order <- order(col_consistency, decreasing = TRUE)

# Create final data matrix with proper ordering
final_data <- scaled_data[top, col_order]

# Create breaks ensuring proper handling of zeros
breaks <- c(-Inf, seq(min(final_data[final_data > 0]), max(final_data), length.out = 99))

# Extract matched annotations based on column names of final_data
matched_cols <- colnames(final_data)

# Create annotation column data frame with all three annotations
annotation_col <- data.frame(
    Type = factor(
        annotation_with_neoexpression$Type[match(matched_cols, rownames(annotation_with_neoexpression))],
        levels = c("wild", "mut")
    ),
    Stage = factor(
        annotation_with_neoexpression$Stage[match(matched_cols, rownames(annotation_with_neoexpression))],
        levels = c("I", "II", "III", "IV")
    ),
    NeoExpressionLevel = factor(
        annotation_with_neoexpression$NeoExpressionLevel[match(matched_cols, rownames(annotation_with_neoexpression))],
        levels = c("High", "Medium", "Low")
    ),
    row.names = matched_cols
)

# Define colors for annotations
ann_colors <- list(
    Type = c(wild = "#1B9E77", mut = "#D95F02"),  # Colors for Type
    Stage = c(  # Colors for Stage with sequential blues
        "I" = "#4575B4",    # Light blue
        "II" = "#74ADD1",   # Medium blue
        "III" = "#ABD9E9",  # Dark blue
        "IV" = "#E0F3F8"    # Very dark blue
    ),
    NeoExpressionLevel = c(  # Colors for NeoExpressionLevel - orange-red family with more contrast
        "High" = "#FF0000",     # Bright red for High
        "Medium" = "#FF8C00",   # Dark orange for Medium
        "Low" = "#FFCC80"       # Light orange for Low
    )
)

# Create ordered indices based on Type, Stage, and NeoExpressionLevel
# First get wild samples
wild_indices <- which(annotation_col$Type == "wild")
wild_ordered <- wild_indices[order(
    annotation_col$Stage[wild_indices], 
    annotation_col$NeoExpressionLevel[wild_indices]
)]

# Then get mut samples
mut_indices <- which(annotation_col$Type == "mut")
mut_ordered <- mut_indices[order(
    annotation_col$Stage[mut_indices], 
    annotation_col$NeoExpressionLevel[mut_indices]
)]

# Combine the ordered indices
new_order <- c(wild_ordered, mut_ordered)

# Reorder the data and annotations
final_data_ordered <- final_data[, new_order]
annotation_col_ordered <- annotation_col[new_order, , drop=FALSE]

# Create the pheatmap object without plotting
pheatmap_obj <- pheatmap(
    final_data_ordered, 
    scale = "none", 
    cluster_cols = FALSE, 
    clustering_distance_rows = "correlation", 
    show_rownames = FALSE, 
    show_colnames = TRUE, 
    angle_col = 90, 
    fontsize_row = 6, 
    fontsize_col = 4, 
    color = c("navy", colorRampPalette(c("#2F2", "red", "orange"))(99)), 
    breaks = breaks, 
    margins = c(12, 12), 
    annotation_col = annotation_col_ordered, 
    annotation_colors = ann_colors, 
    annotation_legend_side = "bottom", 
    annotation_names_side = "right", 
    legend = TRUE, 
    legend_position = "bottom",
    silent = TRUE  # Don't display the plot yet
)

# Convert pheatmap object to a grob for better rendering
pheatmap_grob <- grid.grabExpr(grid.draw(pheatmap_obj$gtable))

# Create a ggplot object with the pheatmap grob
p <- ggplot() + 
     annotation_custom(pheatmap_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
     theme_void()

# Save as PNG with ggsave
ggsave(
    filename = "NeoantigenExpression.png", 
    plot = p, 
    width = 12, 
    height = 6, 
    units = "in", 
    dpi = 300
)

# Save as SVG (vector format for publication)
ggsave(
    filename = "NeoantigenExpression.svg", 
    plot = p, 
    width = 12, 
    height = 6, 
    units = "in",
    device = "svg"
)
