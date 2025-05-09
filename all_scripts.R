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

ls /fs/scratch/PAS1348/pvacseqoutput_annotated_somatic/ | grep -f - SLIDlist_pole_mutated_wild.txt | while IFS=$'\t' read -r _ i j _; do (echo -e "Gene_ID\tEnsembl_ID\t$(echo $j | tr '-' '_')"; cut -f 10 /fs/scratch/PAS1348/pvacseqoutput_annotated_somatic/$j/MHC_Class_I/$j*.filtered.tsv | tail -n +2 | grep -f - /fs/scratch/PAS1348/orien/2024_01_10/RNAseq/gene_and_transcript_expression_results/*$i*.genes.results | grep -v '_PAR_' | awk 'BEGIN{FS="\t"}{print $1"\t"$3"\t"$8}' | awk 'BEGIN{FS=OFS="\t"}{sub(/.[0-9][0-9]$/, "", $2); sub(/.[0-9]$/, "", $2)}1') > /fs/scratch/PAS1348/pvacseqoutput_annotated_somatic/${j%}.filtered.genecount.TPM.results; done

setwd("/fs/scratch/PAS1348/pvacseqoutput_annotated_somatic")
library(tidyverse)
files <- list.files(path = "/fs/scratch/PAS1348/pvacseqoutput_annotated_somatic/", pattern = "*.filtered.genecount.TPM.results$", full.names = TRUE)
merged_data <- files %>% map(~ {df <- read.table(., header = TRUE, sep = "\t") 
    df$Gene_ID <- as.character(df$Gene_ID)
    df$Ensembl_ID <- as.character(df$Ensembl_ID)
    return(df)
    }) %>% 
    reduce(full_join, by = c('Gene_ID', 'Ensembl_ID'))
merged_data1 <- merged_data %>% mutate(across(where(is.logical), as.numeric)) %>% mutate(across(everything(), .fns = ~replace_na(.,0)))
merged_data1 <- merged_data1 %>% mutate_if(is.numeric, round)
write.table(merged_data1, file = "merged.filtered.neogenecount.TPM.results", quote = FALSE, sep = "\t", row.names = FALSE)

# Read neoantigen representing gene TPM file
library(readr)
mergedTPM <- read_delim("merged.filtered.neogenecount.TPM.results.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
View(mergedTPM)
counts_ge_2 <- rowSums(mergedTPM[, 3:ncol(mergedTPM)] >= 2)
counts_le_2 <- rowSums(mergedTPM[, 3:ncol(mergedTPM)] <= 2)
filtered_mergedTPM <- mergedTPM[counts_ge_2 > 1 & counts_le_2 < ncol(mergedTPM[, 3:ncol(mergedTPM)]), ]
print(dim(filtered_mergedTPM))
[1] 1549  172
View(filtered_mergedTPM)
counts_ge_2 <- colSums(filtered_mergedTPM[, 3:ncol(filtered_mergedTPM)] >= 2)
counts_lt_2 <- colSums(filtered_mergedTPM[, 3:ncol(filtered_mergedTPM)] < 2)
keep_cols <- counts_ge_2 > 1 & counts_lt_2 < nrow(filtered_mergedTPM)
final_keep_cols <- c(rep(TRUE, 2), keep_cols)
final_mergedTPM <- filtered_mergedTPM[, final_keep_cols]
print(dim(final_mergedTPM))
[1] 1549  105
 

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
# Replace NA values with "Low" in the NeoExpressionLevel column
# annotation_with_expression$NeoExpressionLevel[is.na(annotation_with_expression$NeoExpressionLevel)] <- "Low"

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



#process gene results files
#files location /fs/scratch/PAS1348/orien/2024_01_10/RNAseq/gene_and_transcript_expression_results/

library(tidyverse)

# Function to read and prepare SLID mapping
read_slid_mapping <- function(slid_file) {
  read_delim(slid_file, 
             delim = "\t", 
             col_names = c("file_id", "sample_id"),
             col_types = cols(.default = "c"),
             trim_ws = TRUE) %>%
    mutate(
      file_id = str_trim(file_id),
      sample_id = str_trim(sample_id)
    )
}

# Function to check if file has the correct format
is_valid_gene_file <- function(file_path) {
  # Check if file name matches the expected pattern
  if (!str_detect(basename(file_path), "\\.genes\\.results$")) {
    return(FALSE)
  }
  
  # Try to read the first few lines and check columns
  tryCatch({
    df <- read_delim(file_path, delim = "\t", n_max = 1, 
                     escape_double = FALSE, trim_ws = TRUE, 
                     show_col_types = FALSE)
    return(all(c("gene_symbol", "gene_id", "expected_count") %in% colnames(df)))
  }, error = function(e) {
    return(FALSE)
  })
}

# Function to process a single gene results file
process_gene_file <- function(file_path, sample_id) {
  # Read the file
  df <- read_delim(file_path, 
                   delim = "\t", 
                   escape_double = FALSE, 
                   trim_ws = TRUE,
                   show_col_types = FALSE)
  
  # Select required columns, clean gene_id, and convert expected_count to integer
  df <- df %>%
    select(gene_symbol, gene_id, expected_count) %>%
    mutate(
      gene_id = str_replace(gene_id, "\\.[0-9]+$", ""),  # Remove the version number after the dot
      expected_count = as.integer(round(expected_count))  # Round and convert to integer
    ) %>%
    rename_with(~sample_id, .cols = "expected_count")  # Rename expected_count column with the mapped sample ID
  
  return(df)
}

# Process all files in the directory
process_all_files <- function(directory_path, slid_file, pattern = "\\.genes\\.results$") {
  # Read SLID mapping
  slid_mapping <- read_slid_mapping(slid_file)
  message("Loaded SLID mapping with ", nrow(slid_mapping), " entries")
  
  # Get list of all gene result files
  files <- list.files(
    path = directory_path,
    pattern = pattern,
    full.names = TRUE
  )
  
  # Filter for valid files
  valid_files <- files[sapply(files, is_valid_gene_file)]
  
  if(length(valid_files) == 0) {
    stop("No valid gene results files found in the directory")
  }
  
  # Create a named vector for file mapping
  file_mapping <- setNames(
    valid_files,
    basename(valid_files) %>%
      str_extract("^[^.]+") %>%
      str_remove("R$")
  )
  
  # Filter files based on SLID mapping and create processing pairs
  processing_pairs <- slid_mapping %>%
    filter(file_id %in% names(file_mapping)) %>%
    mutate(
      full_path = file_mapping[file_id],
      sample_id = str_replace_all(sample_id, "-", "_")  # Replace hyphens with underscores in sample IDs
    ) %>%
    select(full_path, sample_id)
  
  # Report on filtered files
  unmapped_files <- setdiff(names(file_mapping), slid_mapping$file_id)
  if(length(unmapped_files) > 0) {
    message("\nSkipping ", length(unmapped_files), " files not found in SLID mapping:")
    message(paste("- ", unmapped_files, collapse = "\n"))
  }
  
  if(nrow(processing_pairs) == 0) {
    stop("No files remaining after filtering for SLID mapping")
  }
  
  # Print the files being processed with their mapped IDs
  message("\nProcessing ", nrow(processing_pairs), " files with the following mappings:")
  for(i in 1:nrow(processing_pairs)) {
    message("- ", basename(processing_pairs$full_path[i]), " -> ", processing_pairs$sample_id[i])
  }
  
  # Process each file and combine results
  result <- map2(
    processing_pairs$full_path,
    processing_pairs$sample_id,
    ~process_gene_file(.x, .y)
  ) %>%
    reduce(full_join, by = c("gene_symbol", "gene_id"))
  
  return(result)
}

# Usage example:
directory_path <- "/fs/scratch/PAS1348/orien/2024_01_10/RNAseq/gene_and_transcript_expression_results/"
slid_file <- "/fs/ess/PAS1348/nagesh/orien/scripts/SLIDlist.txt"
combined_data <- process_all_files(directory_path, slid_file)

# Immune_infiltration_heatmap
# Load required libraries
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Step 1: Data preparation
# Extract sample IDs from annotation dataframe
annotation_with_expression$SampleID <- rownames(annotation_with_expression)
annotation_ids <- annotation_with_expression$SampleID

# Filter cibersortx data for matching IDs
# Create a clean ID column by removing any potential whitespace
cibersortx <- cibersortx %>% 
    mutate(ID = trimws(Mixture))

# Filter for matching IDs
matched_cibersortx <- cibersortx %>%
    filter(ID %in% annotation_ids)

# Extract immune cell columns (specifically columns 2:23)
immune_columns <- colnames(cibersortx)[2:23]

# Step 2: Restructure the data to swap rows and columns
# Create a matrix with immune cells as rows and samples as columns
immune_data_transposed <- matched_cibersortx %>%
    select(Mixture, all_of(immune_columns)) %>%
    column_to_rownames("Mixture") %>%
    t() %>%
    as.matrix()

# Use absolute values directly without log transformation
immune_matrix_abs <- immune_data_transposed

# Step 3: Prepare annotation dataframe for columns (formerly rows)
# Filter annotation data for matching IDs and create a new dataframe
matched_annotation <- annotation_with_expression %>%
    filter(SampleID %in% colnames(immune_matrix_abs))

# Convert NeoExpressionLevel and Type to factors with specific orders
matched_annotation <- matched_annotation %>%
    mutate(
        NeoExpressionLevel = factor(NeoExpressionLevel, levels = c("High", "Medium", "Low")),
        Type = factor(Type, levels = c("wild", "mut")),
        Stage = factor(Stage, levels = c("I", "II", "III", "IV"))
    )

# Sort annotation data in a specific order (first by Type, then by Stage, then by NeoExpressionLevel)
matched_annotation <- matched_annotation %>%
    arrange(Type, Stage, NeoExpressionLevel)

# Create annotation dataframe
annotation_df <- matched_annotation %>%
    select(Type, Stage, NeoExpressionLevel)

# Set row names manually
rownames(annotation_df) <- matched_annotation$SampleID

# Reorder the columns of the matrix to match the sorted annotation order
immune_matrix_abs <- immune_matrix_abs[, rownames(annotation_df)]

# Step 4: Define colors for annotations
# Colors for NeoExpressionLevel
neo_colors <- c("High" = "#FF0000", "Medium" = "#FF8C00", "Low" = "#FFCC80")

# Colors for Stage
stage_colors <- c("I" = "#4575B4", "II" = "#74ADD1", "III" = "#ABD9E9", "IV" = "#E0F3F8")

# Colors for Type
type_colors <- c("wild" = "#1B9E77", "mut" = "#D95F02")

# Combine all annotation colors
ann_colors <- list(
    NeoExpressionLevel = neo_colors,
    Stage = stage_colors,
    Type = type_colors
)

# Step 5: Create the heatmap with absolute values
# Calculate breaks for the color scale - ensure they are unique
min_val <- min(immune_matrix_abs)
max_val <- max(immune_matrix_abs)

# Create a color palette that transitions from navy (low) to white (medium) to red (high)
# Fix the breaks to ensure they are unique
if(min_val == 0) {
    # If min value is 0, create breaks starting slightly above 0
    breaks <- c(
        0,
        seq(0.00001, max_val, length.out = 99)
    )
} else {
    # Calculate midpoint
    mid_val <- (min_val + max_val) / 2
    # Create breaks with slight offset to ensure uniqueness
    breaks <- unique(c(
        seq(min_val, mid_val, length.out = 50),
        seq(mid_val + 0.00001, max_val, length.out = 50)
    ))
}

pheatmap(
    mat = immune_matrix_abs,
    color = colorRampPalette(c("white", "red"))(100),
    breaks = breaks,
    annotation_col = annotation_df,  # Now using column annotations instead of row
    annotation_colors = ann_colors,
    cluster_cols = FALSE,  # Don't cluster columns (keep samples in sorted order)
    cluster_rows = FALSE,  # Don't cluster rows (immune cells)
    show_rownames = TRUE,  # Show immune cell types (now rows)
    show_colnames = TRUE,  # Show sample IDs (now columns)
    fontsize_row = 8,      # Adjust for immune cell names
    fontsize_col = 4,      # Make column font smaller to fit all sample IDs
    main = paste0("Immune Cell Infiltration (n=", ncol(immune_matrix_abs), ")"),
    legend = TRUE,
    annotation_legend = TRUE, 
    filename = "immune_infiltration_heatmap.png",
    border_color = NA,  # Remove cell borders
    cellwidth = NA,     # Auto-calculate cell width
    cellheight = 12,    # Fixed cell height
    gaps_col = c(),     # No column gaps
    gaps_row = c(),     # No row gaps
    width = 20,
    height = 10,
    res = 300
)

#Violin plots split by stage, type, and assign color annotation to mean expression levels (high, medium, and low)
library(readr)      # For reading delimited files
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

# Using cibersortx and annotation_with_neoexpression data frames
# Prepare the data by reshaping from wide to long format and joining with annotations
cibersortx <- read_delim("CIBERSORTx_Job1_Results.csv", delim = ",", col_names = T)
cell_columns <- 2:23  # From "B cells naive" to "Neutrophils"
cibersortx_long <- melt(cibersortx[, c(1, cell_columns)], 
                        id.vars = "Mixture", 
                        variable.name = "Cell_Type", 
                        value.name = "Fraction")

# Join with annotation_with_expression data
cibersortx_with_anno <- right_join(
  cibersortx_long,
  cbind(Mixture = rownames(annotation_with_neoexpression), as.data.frame(annotation_with_neoexpression)),
  by = "Mixture"
)

# Check if join worked
if(any(is.na(cibersortx_with_anno$Type)) || any(is.na(cibersortx_with_anno$Stage)) || any(is.na(cibersortx_with_anno$NeoExpressionLevel))) {
  warning("Some samples couldn't be matched with annotation data")
  print(table(is.na(cibersortx_with_anno$Type)))
  print(table(is.na(cibersortx_with_anno$Stage)))
  print(table(is.na(cibersortx_with_anno$NeoExpressionLevel)))
}

# Remove any rows with NA values
cibersortx_with_anno <- cibersortx_with_anno %>% 
  filter(!is.na(Fraction) & !is.na(Type) & !is.na(Stage) & !is.na(NeoExpressionLevel))

# Define colors for NeoExpressionLevel
expression_colors <- c("Low" = "#3498db", "Medium" = "#f39c12", "High" = "#e74c3c")

# Create directory for the new plots
dir.create("type_stage_expression_plots1", showWarnings = FALSE)

# Function to create violin plots split by Type, Stage, and colored by NeoExpressionLevel
create_type_stage_expression_violin_plots <- function() {
  # Get unique cell types
  cell_types <- unique(cibersortx_with_anno$Cell_Type)
  
  # For each cell type, create a plot that splits by Type, then by Stage
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    
    # Subset data for the current cell type
    subset_data <- cibersortx_with_anno %>% filter(Cell_Type == cell_type)
    
    # Get unique types and stages for consistent ordering
    types <- c("wild", "mut")
    stages <- sort(unique(subset_data$Stage))
    
    # Create separate plots for wild and mut types
    plot_list <- list()
    
    for (type_val in types) {
      # Filter data for this type
      type_data <- subset_data %>% filter(Type == type_val)
      
      if (nrow(type_data) > 0) {
        # Create plot for this type with stages, using uniform color for violins
        # and highlighting only the points by NeoExpressionLevel
        p <- ggplot(type_data, aes(x = Stage, y = Fraction)) +
          geom_violin(alpha = 0.7, fill = "#CCCCCC") +  # Uniform light gray for all violins
          geom_boxplot(width = 0.1, alpha = 0.5, fill = "white", outlier.shape = NA) +
          geom_jitter(aes(color = NeoExpressionLevel), width = 0.1, alpha = 0.9, size = 2.5) +
          scale_x_discrete(limits = stages) +
          scale_color_manual(values = expression_colors) +
          labs(title = paste(cell_type, "-", type_val, "Type"),
               y = "Cell Fraction",
               x = "Stage") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 11),
                legend.position = "bottom",
                legend.title = element_text(size = 9),
                legend.text = element_text(size = 8))
        
        plot_list[[type_val]] <- p
      }
    }
    
    # Arrange the plots side by side
    if (length(plot_list) == 2) {
      combined_plot <- gridExtra::grid.arrange(
        plot_list[["wild"]], plot_list[["mut"]], 
        ncol = 2,
        top = grid::textGrob(paste("Distribution of", cell_type, "by Stage and NeoExpressionLevel"), 
                            gp = grid::gpar(fontsize = 14))
      )
    } else if (length(plot_list) == 1) {
      # If we only have one type, just use that plot
      combined_plot <- plot_list[[1]] +
        labs(title = paste("Distribution of", cell_type, "by Stage and NeoExpressionLevel"))
    } else {
      next  # Skip if no plots were created
    }
    
    # Save the combined plot
    filename <- paste0("type_stage_expression_plots1/", gsub(" ", "_", cell_type), "_type_stage_expression.png")
    ggsave(filename, combined_plot, width = 12, height = 6, dpi = 300)
    
    # Print status message
    cat("Saved type-stage-expression plot for", cell_type, "\n")
  }
}

# Create the type-stage-expression split violin plots
create_type_stage_expression_violin_plots()

#process gene results files expected_count matrix
#files location /fs/scratch/PAS1348/orien/2024_01_10/RNAseq/gene_and_transcript_expression_results/

library(tidyverse)

# Function to read and prepare SLID mapping
read_slid_mapping <- function(slid_file) {
  read_delim(slid_file, 
             delim = "\t", 
             col_names = c("file_id", "sample_id"),
             col_types = cols(.default = "c"),
             trim_ws = TRUE) %>%
    mutate(
      file_id = str_trim(file_id),
      sample_id = str_trim(sample_id)
    )
}

# Function to check if file has the correct format
is_valid_gene_file <- function(file_path) {
  # Check if file name matches the expected pattern
  if (!str_detect(basename(file_path), "\\.genes\\.results$")) {
    return(FALSE)
  }
  
  # Try to read the first few lines and check columns
  tryCatch({
    df <- read_delim(file_path, delim = "\t", n_max = 1, 
                     escape_double = FALSE, trim_ws = TRUE, 
                     show_col_types = FALSE)
    return(all(c("gene_symbol", "gene_id", "expected_count") %in% colnames(df)))
  }, error = function(e) {
    return(FALSE)
  })
}

# Function to process a single gene results file
process_gene_file <- function(file_path, sample_id) {
  # Read the file
  df <- read_delim(file_path, 
                   delim = "\t", 
                   escape_double = FALSE, 
                   trim_ws = TRUE,
                   show_col_types = FALSE)
  
  # Select required columns, clean gene_id, and convert expected_count to integer
  df <- df %>%
    select(gene_symbol, gene_id, expected_count) %>%
    mutate(
      gene_id = str_replace(gene_id, "\\.[0-9]+$", ""),  # Remove the version number after the dot
      expected_count = as.integer(round(expected_count))  # Round and convert to integer
    ) %>%
    rename_with(~sample_id, .cols = "expected_count")  # Rename expected_count column with the mapped sample ID
  
  return(df)
}

# Process all files in the directory
process_all_files <- function(directory_path, slid_file, pattern = "\\.genes\\.results$") {
  # Read SLID mapping
  slid_mapping <- read_slid_mapping(slid_file)
  message("Loaded SLID mapping with ", nrow(slid_mapping), " entries")
  
  # Get list of all gene result files
  files <- list.files(
    path = directory_path,
    pattern = pattern,
    full.names = TRUE
  )
  
  # Filter for valid files
  valid_files <- files[sapply(files, is_valid_gene_file)]
  
  if(length(valid_files) == 0) {
    stop("No valid gene results files found in the directory")
  }
  
  # Create a named vector for file mapping
  file_mapping <- setNames(
    valid_files,
    basename(valid_files) %>%
      str_extract("^[^.]+") %>%
      str_remove("R$")
  )
  
  # Filter files based on SLID mapping and create processing pairs
  processing_pairs <- slid_mapping %>%
    filter(file_id %in% names(file_mapping)) %>%
    mutate(
      full_path = file_mapping[file_id],
      sample_id = str_replace_all(sample_id, "-", "_")  # Replace hyphens with underscores in sample IDs
    ) %>%
    select(full_path, sample_id)
  
  # Report on filtered files
  unmapped_files <- setdiff(names(file_mapping), slid_mapping$file_id)
  if(length(unmapped_files) > 0) {
    message("\nSkipping ", length(unmapped_files), " files not found in SLID mapping:")
    message(paste("- ", unmapped_files, collapse = "\n"))
  }
  
  if(nrow(processing_pairs) == 0) {
    stop("No files remaining after filtering for SLID mapping")
  }
  
  # Print the files being processed with their mapped IDs
  message("\nProcessing ", nrow(processing_pairs), " files with the following mappings:")
  for(i in 1:nrow(processing_pairs)) {
    message("- ", basename(processing_pairs$full_path[i]), " -> ", processing_pairs$sample_id[i])
  }
  
  # Process each file and combine results
  result <- map2(
    processing_pairs$full_path,
    processing_pairs$sample_id,
    ~process_gene_file(.x, .y)
  ) %>%
    reduce(full_join, by = c("gene_symbol", "gene_id"))
  
  return(result)
}

# Usage example:
directory_path <- "/fs/scratch/PAS1348/orien/2024_01_10/RNAseq/gene_and_transcript_expression_results/"
slid_file <- "/fs/ess/PAS1348/nagesh/orien/scripts/SLIDlist.txt"
combined_data <- process_all_files(directory_path, slid_file)

#Differential expression analysis wildtype Vs Mut
library(limma)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

# Create output directory
dir.create("results", showWarnings = FALSE)

# Prepare count data from combined_data
# Convert tibble to data frame and remove duplicates
counts <- as.data.frame(combined_data[!duplicated(combined_data$gene_symbol), ])

# Set gene symbols as row names
rownames(counts) <- counts$gene_symbol

# Remove gene_symbol and gene_id columns
counts <- counts[, !(colnames(counts) %in% c("gene_symbol", "gene_id"))]

# Convert to matrix
counts <- as.matrix(counts)

# Verify sample names match between count data and metadata
if(!all(colnames(counts) %in% sample_info$sample_id)) {
    stop("Sample names in count data don't match metadata")
}

# Order sample info to match count data column order
sample_info <- sample_info[match(colnames(counts), sample_info$sample_id), ]

# Create simple design matrix for type comparison
design_simple <- model.matrix(~type, data = sample_info)
colnames(design_simple) <- c("Intercept", "typemut")

# Create DGEList object
dge <- DGEList(counts = counts)

# Filter low expression genes
keep <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[keep, ]
cat("Number of genes after filtering:", nrow(dge), "\n")

# Calculate normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# Run voom transformation
v <- voom(dge, design_simple, plot = FALSE)

# Fit linear model
fit <- lmFit(v, design_simple)

# Create contrast matrix
cont.matrix <- makeContrasts(
    mutvswild = typemut,
    levels = design_simple
)

# Fit contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Get results
results <- topTable(fit2, coef = "mutvswild", n = Inf)
cat("\nNumber of total results:", nrow(results), "\n")

# Filter significant genes
significant <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
cat("Number of significant genes:", nrow(significant), "\n")

# Write results to files
write.table(results, 
           file = "results/mutvswild_all_results.txt", 
           sep = "\t", 
           quote = FALSE)

write.table(significant, 
           file = "results/mutvswild_significant_genes.txt", 
           sep = "\t", 
           quote = FALSE)

# Generate MA plot
pdf("results/mutvswild_MA_plot.pdf", width = 10, height = 8)
limma::plotMA(fit2, coef = "mutvswild")
abline(h = c(-1, 1), col = "red", lty = 2)
dev.off()

# Generate Volcano plot
pdf("results/mutvswild_volcano_plot.pdf", width = 12, height = 10)

# Get top genes for labeling
top_genes <- results[with(results, 
    adj.P.Val < 0.05 & abs(logFC) > 1), ]
top_genes <- top_genes[order(top_genes$adj.P.Val), ]
top_genes <- head(top_genes, 30)  # Label only top 30 most significant genes

print(EnhancedVolcano(results,
    lab = rownames(results),
    x = 'logFC',
    y = 'adj.P.Val',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 1.0,
    labSize = 3.0,
    title = 'Mutant vs Wild Type',
    subtitle = paste('Total significant genes:', nrow(significant)),
    legendPosition = 'right',
    legendLabSize = 14,  # Increased legend font size
    col = c('grey', 'blue', 'red', 'green'),
    colAlpha = 0.4,
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    max.overlaps = 40,
    selectLab = rownames(top_genes),
    boxedLabels = TRUE,
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'partial',
    borderWidth = 1.5,
    borderColour = 'black'))
dev.off()

# Generate summary statistics
summary_stats <- data.frame(
    Category = c("Total_Genes", "Significant_Genes", "Upregulated", "Downregulated"),
    Count = c(
        nrow(results),
        nrow(significant),
        sum(significant$logFC > 0),
        sum(significant$logFC < 0)
    )
)

# Write summary statistics
write.table(summary_stats, 
           "results/mutvswild_summary.txt", 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE)

# Generate heatmap for top 100 genes
top_genes <- results[order(results$adj.P.Val), ]
top_genes <- head(top_genes, 100)
top_genes_names <- rownames(top_genes)

cat("Number of selected genes:", length(top_genes_names), "\n")

# Extract expression data for top genes
expr_matrix <- v$E[top_genes_names, ]

# Order samples by Type and Stage
sample_order <- order(sample_info$type, sample_info$stage)

# Reorder expression matrix columns
expr_matrix <- expr_matrix[, sample_order]

# Create annotation dataframe with ordered samples
annotation_col <- data.frame(
    Stage = sample_info$stage[sample_order],
    Type = sample_info$type[sample_order],
    row.names = sample_info$sample_id[sample_order]
)

# Create color schemes
stage_colors <- brewer.pal(n = min(5, 9), name = "Set1")
type_colors <- brewer.pal(n = 3, name = "Set2")[1:2]

ann_colors <- list(
    Stage = setNames(stage_colors, levels(sample_info$stage)),
    Type = setNames(type_colors, levels(sample_info$type))
)

# Calculate dimensions for square cells
n_rows <- nrow(expr_matrix)
n_cols <- ncol(expr_matrix)
cell_size <- 20 
width <- n_cols * cell_size + 300 
height <- n_rows * cell_size + 300 

# Generate heatmap with PNG format
png("results/top_100_genes_heatmap.png", 
    width = width,     
    height = height,    
    res = 300)     

pheatmap(expr_matrix,
         scale = "row",  
         show_rownames = TRUE,  
         show_colnames = TRUE,  
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = "Top 100 Most Significant Genes\n(Mutant vs Wild Type)",
         fontsize_row = 4,      
         fontsize_col = 4,
         fontsize_main = 16,     
         fontsize = 12,         
         cluster_cols = FALSE,   
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         border_color = NA,
         angle_col = 90,      
         annotation_names_col = TRUE,
         annotation_names_row = TRUE,
         annotation_legend = TRUE,
         cellwidth = cell_size/5,    
         cellheight = cell_size/5,   
         treeheight_row = 0,        
         treeheight_col = 0) 

dev.off()
 

#Heatmap for specific gene set
library(limma)
library(edgeR)
library(pheatmap)
library(RColorBrewer)

# Create output directory
dir.create("results_geneset", showWarnings = FALSE)

# Define the gene set of interest
target_genes <- c("CGAS", "TMEM173", "MAVS", "IFI27", 
                 "IFI44L", "IFIT1", "ISG15", "RSAD2", "SIGLEC1")

# Prepare count data from combined_data
# Convert tibble to data frame and remove duplicates
counts <- as.data.frame(combined_data[!duplicated(combined_data$gene_symbol), ])

# Set gene symbols as row names
rownames(counts) <- counts$gene_symbol

# Remove gene_symbol and gene_id columns
counts <- counts[, !(colnames(counts) %in% c("gene_symbol", "gene_id"))]

# Convert to matrix
counts <- as.matrix(counts)

# Filter for target genes
counts_subset <- counts[rownames(counts) %in% target_genes, ]

# Verify all genes were found
found_genes <- rownames(counts_subset)
missing_genes <- setdiff(target_genes, found_genes)
if(length(missing_genes) > 0) {
    warning("Missing genes: ", paste(missing_genes, collapse = ", "))
}

# Create DGEList object
dge <- DGEList(counts = counts_subset)

# Calculate normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# Create simple design matrix for type comparison
design_simple <- model.matrix(~type, data = sample_info)
colnames(design_simple) <- c("Intercept", "typemut")

# Run voom transformation
v <- voom(dge, design_simple, plot = FALSE)

# Extract expression data
expr_matrix <- v$E

# Order samples by Type and Stage
sample_order <- order(sample_info$type, sample_info$stage)

# Reorder expression matrix columns
expr_matrix <- expr_matrix[, sample_order]

# Create annotation dataframe with ordered samples
annotation_col <- data.frame(
    Stage = sample_info$stage[sample_order],
    Type = sample_info$type[sample_order],
    row.names = sample_info$sample_id[sample_order]
)

# Create color schemes
stage_colors <- brewer.pal(n = min(5, 9), name = "Set1")
type_colors <- brewer.pal(n = 3, name = "Set2")[1:2]

ann_colors <- list(
    Stage = setNames(stage_colors, levels(sample_info$stage)),
    Type = setNames(type_colors, levels(sample_info$type))
)

# Calculate dimensions for square cells
n_rows <- nrow(expr_matrix)
n_cols <- ncol(expr_matrix)
cell_size <- 25  # increased pixels per cell
width <- n_cols * cell_size + 200  # adjusted padding
height <- n_rows * cell_size + 200  # adjusted padding

# Generate heatmap
png(file.path("results_geneset", "geneset_heatmap.png"),
    width = width,
    height = 900,
    res = 300)

pheatmap(expr_matrix,
         scale = "row",  # Scale by row to show relative expression
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = "Expression Heatmap of Selected Genes",
         fontsize_row = 4,  
         fontsize_col = 4,   
         fontsize_main = 14,
         cluster_cols = FALSE, 
         cluster_rows = TRUE,  
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         border_color = NA,
         angle_col = 90,
         annotation_names_col = TRUE,
         annotation_names_row = TRUE,
         annotation_legend = TRUE,
         cellwidth = cell_size/5,
         cellheight = cell_size/5)

dev.off()
 

#Whisker plots - Part 2
# This script is a continuation of the heatmap analysis and generates whisker plots
# for the same set of genes from the previous analysis

# Load required libraries
library(limma)
library(edgeR)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(cowplot)

# Create output directory
dir.create("results_geneset", showWarnings = FALSE)

# Define the gene set of interest (using the same set as in heatmap script)
target_genes <- c("CGAS", "TMEM173", "MAVS", "IFI27", 
                 "IFI44L", "IFIT1", "ISG15", "RSAD2", "SIGLEC1")

# Reorder genes to match the original order
gene_order <- target_genes

# Prepare count data from combined_data
# Convert tibble to data frame and remove duplicates
counts <- as.data.frame(combined_data[!duplicated(combined_data$gene_symbol), ])

# Set gene symbols as row names
rownames(counts) <- counts$gene_symbol

# Remove gene_symbol and gene_id columns
counts <- counts[, !(colnames(counts) %in% c("gene_symbol", "gene_id"))]

# Convert to matrix
counts <- as.matrix(counts)

# Filter for target genes and maintain the original order
counts_subset <- counts[rownames(counts) %in% target_genes, ]
counts_subset <- counts_subset[gene_order[gene_order %in% rownames(counts_subset)], ]

# Verify all genes were found
found_genes <- rownames(counts_subset)
missing_genes <- setdiff(target_genes, found_genes)
if(length(missing_genes) > 0) {
    warning("Missing genes: ", paste(missing_genes, collapse = ", "))
}

# Create DGEList object
dge <- DGEList(counts = counts_subset)

# Calculate normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# Create simple design matrix for type comparison
design_simple <- model.matrix(~type, data = sample_info)
colnames(design_simple) <- c("Intercept", "typemut")

# Run voom transformation
v <- voom(dge, design_simple, plot = FALSE)

# Extract expression data
expr_matrix <- v$E

# Prepare data for plotting
plot_data <- as.data.frame(expr_matrix) %>%
    mutate(Gene = rownames(expr_matrix)) %>%
    gather(key = "Sample", value = "Expression", -Gene)

# Add sample information and filter out Uk stage
plot_data <- merge(plot_data,
                   data.frame(Sample = sample_info$sample_id,
                              Type = sample_info$type,
                              Stage = sample_info$stage),
                   by = "Sample") %>%
    filter(Stage != "Uk")  # Remove Uk stage

# If annotation_with_expression is available, add NeoantigenExpression information
if(exists("annotation_with_expression")) {
    # Convert IDs to character to ensure proper matching
    annotation_with_expression$WES <- as.character(annotation_with_expression$WES)
    annotation_with_expression$NeoantigenExpresion <- as.character(annotation_with_expression$NeoantigenExpresion)
    
    # Create mapping from sample ID to NeoantigenExpression
    neoantigen_map <- setNames(
        annotation_with_expression$NeoantigenExpresion,
        annotation_with_expression$WES
    )
    
    # Add NeoantigenExpression column
    plot_data$NeoantigenExpression <- neoantigen_map[plot_data$Sample]
    
    # For mut samples with NA NeoantigenExpression, set to "Low"
    is_mut_na <- plot_data$Type == "mut" & is.na(plot_data$NeoantigenExpression)
    if(any(is_mut_na)) {
        plot_data$NeoantigenExpression[is_mut_na] <- "Low"
    }
    
    # Create Type_NeoAntigen column for visualization
    plot_data$Type_NeoAntigen <- as.character(plot_data$Type)
    is_mut <- plot_data$Type == "mut"
    plot_data$Type_NeoAntigen[is_mut] <- paste0("mut_", plot_data$NeoantigenExpression[is_mut])
    
    # Make sure Type_NeoAntigen is a factor with the correct order
    plot_data$Type_NeoAntigen <- factor(
        plot_data$Type_NeoAntigen,
        levels = c("wild", "mut_High", "mut_Medium", "mut_Low")
    )
    
    # Print summary of Type_NeoAntigen categories
    cat("Type_NeoAntigen categories:\n")
    print(table(plot_data$Type_NeoAntigen))
}

# Set factor levels for other variables
plot_data$Gene <- factor(plot_data$Gene, levels = target_genes)
plot_data$Stage <- factor(plot_data$Stage, levels = c("I", "II", "III", "IV"))
plot_data$Type <- factor(plot_data$Type)

# Create color scheme
type_colors <- brewer.pal(n = 3, name = "Set2")[1:2]

# Function to format p-value
format_pvalue <- function(p) {
    sapply(p, function(x) {
        if (x < 0.001) return("p < 0.001")
        if (x < 0.01) return(paste0("p = ", sprintf("%.3f", x)))
        return(paste0("p = ", sprintf("%.2f", x)))
    })
}

# Calculate Kruskal-Wallis test results
kw_results <- plot_data %>%
    group_by(Gene) %>%
    summarise(
        KW_pvalue = kruskal.test(Expression ~ interaction(Type, Stage))$p.value,
        .groups = "drop"
    ) %>%
    mutate(
        KW_significant = KW_pvalue < 0.05,
        KW_label = format_pvalue(KW_pvalue),
        Significance = case_when(
            KW_pvalue < 0.001 ~ "***",
            KW_pvalue < 0.01 ~ "**",
            KW_pvalue < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    )

# Function to create whisker plot for a single gene
create_gene_plot <- function(data, gene, kw_results, show_legend = FALSE) {
    # Get KW test results for this gene
    kw_result <- kw_results %>% filter(Gene == gene)
    
    # Filter data for this gene
    gene_data <- data %>% filter(Gene == gene)
    
    # Create the base plot - using position_dodge for side-by-side boxplots
    if(exists("annotation_with_expression") && "Type_NeoAntigen" %in% colnames(data)) {
        # Use Type_NeoAntigen if available
        p <- ggplot(gene_data, 
               aes(x = Stage, y = Expression, fill = Type_NeoAntigen))
        
        # Define color palette for Type_NeoAntigen
        neo_colors <- c(
            "wild" = "#66C2A5",       # green for wild
            "mut_High" = "#B10026",   # dark red for High
            "mut_Medium" = "#E31A1C", # medium red for Medium
            "mut_Low" = "#FB9A99"     # light red for Low
        )
        
        p <- p + scale_fill_manual(values = neo_colors,
                              name = "Type NeoAntigen")
    } else {
        # Use Type if Type_NeoAntigen is not available
        p <- ggplot(gene_data, 
               aes(x = Stage, y = Expression, fill = Type))
               
        p <- p + scale_fill_manual(values = type_colors,
                              name = "Type")
    }
    
    # Complete the plot with side-by-side boxplots
    p <- p + geom_boxplot(position = position_dodge(width = 0.8),
                   outlier.shape = NA) +
        geom_point(position = position_jitterdodge(dodge.width = 0.8),
                   alpha = 0.4,
                   size = 1,
                   shape = 16) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 10),
            legend.position = if(show_legend) "right" else "none",
            panel.grid.minor = element_blank()
        ) +
        labs(
            title = gene,
            subtitle = paste0(
                "Kruskal-Wallis test: ",
                kw_result$KW_label,
                " ",
                kw_result$Significance
            ),
            x = "Stage",
            y = "Expression Level"
        ) +
        stat_summary(
            fun = median,
            geom = "line",
            aes(group = Type_NeoAntigen),
            position = position_dodge(width = 0.8)
        )
    
    return(p)
}

# Generate individual plots for the grid
gene_plots <- lapply(target_genes, function(gene) {
    create_gene_plot(plot_data, gene, kw_results, show_legend = FALSE)
})

# Get the legend from a temporary plot
legend_plot <- create_gene_plot(plot_data, target_genes[1], kw_results, show_legend = TRUE)
legend <- cowplot::get_legend(legend_plot)

# Create significance text
sig_text <- textGrob(
    "Significance levels: *** p<0.001, ** p<0.01, * p<0.05, ns: not significant",
    gp = gpar(fontsize = 10)
)

# Arrange plots in a grid
plots_arranged <- do.call(gridExtra::arrangeGrob, c(gene_plots, ncol = 3))

# Create bottom panel with legend and significance
bottom <- gridExtra::arrangeGrob(
    legend,
    sig_text,
    ncol = 1,
    heights = c(1, 1)
)

# Combine everything
combined_plot <- gridExtra::grid.arrange(
    plots_arranged,
    bottom,
    heights = c(10, 2),
    ncol = 1
)

# Save combined plot in both PDF and PNG formats
ggsave(
    filename = file.path("results_geneset", "all_genes_whisker_plots.pdf"),
    plot = combined_plot,
    width = 12,
    height = 15,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

ggsave(
    filename = file.path("results_geneset", "all_genes_whisker_plots.png"),
    plot = combined_plot,
    width = 12,
    height = 15,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

# Also save SVG format for better scalability
ggsave(
    filename = file.path("results_geneset", "all_genes_whisker_plots.svg"),
    plot = combined_plot,
    width = 12,
    height = 15,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

image.png

#Whisker plots - Part 2 (Modified)
# This script is a continuation of the heatmap analysis and generates whisker plots
# for the same set of genes, excluding stage I and low neoantigenexpression IDs,
# and combining high and medium neoantigenexpression IDs

# Load required libraries
library(limma)
library(edgeR)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(cowplot)

# Create output directory
dir.create("results_geneset_modified", showWarnings = FALSE)

# Define the gene set of interest (using the same set as in heatmap script)
target_genes <- c("CGAS", "TMEM173", "MAVS", "IFI27", 
                 "IFI44L", "IFIT1", "ISG15", "RSAD2", "SIGLEC1")

# Reorder genes to match the original order
gene_order <- target_genes

# Prepare count data from combined_data
# Convert tibble to data frame and remove duplicates
counts <- as.data.frame(combined_data[!duplicated(combined_data$gene_symbol), ])

# Set gene symbols as row names
rownames(counts) <- counts$gene_symbol

# Remove gene_symbol and gene_id columns
counts <- counts[, !(colnames(counts) %in% c("gene_symbol", "gene_id"))]

# Convert to matrix
counts <- as.matrix(counts)

# Filter for target genes and maintain the original order
counts_subset <- counts[rownames(counts) %in% target_genes, ]
counts_subset <- counts_subset[gene_order[gene_order %in% rownames(counts_subset)], ]

# Verify all genes were found
found_genes <- rownames(counts_subset)
missing_genes <- setdiff(target_genes, found_genes)
if(length(missing_genes) > 0) {
    warning("Missing genes: ", paste(missing_genes, collapse = ", "))
}

# Create DGEList object
dge <- DGEList(counts = counts_subset)

# Calculate normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# Create simple design matrix for type comparison
design_simple <- model.matrix(~type, data = sample_info)
colnames(design_simple) <- c("Intercept", "typemut")

# Run voom transformation
v <- voom(dge, design_simple, plot = FALSE)

# Extract expression data
expr_matrix <- v$E

# Prepare data for plotting
plot_data <- as.data.frame(expr_matrix) %>%
    mutate(Gene = rownames(expr_matrix)) %>%
    gather(key = "Sample", value = "Expression", -Gene)

# Add sample information and filter out Uk stage and stage I
plot_data <- merge(plot_data,
                   data.frame(Sample = sample_info$sample_id,
                              Type = sample_info$type,
                              Stage = sample_info$stage),
                   by = "Sample") %>%
    filter(Stage != "Uk" & Stage != "I")  # Remove Uk stage and Stage I

# If annotation_with_expression is available, add NeoantigenExpression information
if(exists("annotation_with_expression")) {
    # Convert IDs to character to ensure proper matching
    annotation_with_expression$WES <- as.character(annotation_with_expression$WES)
    annotation_with_expression$NeoantigenExpresion <- as.character(annotation_with_expression$NeoantigenExpresion)
    
    # Create mapping from sample ID to NeoantigenExpression
    # Note: WES in annotation_with_expression must match sample_id in sample_info
    neoantigen_map <- setNames(
        annotation_with_expression$NeoantigenExpresion,
        annotation_with_expression$WES
    )
    
    # Add NeoantigenExpression column
    # Using the column name "NeoantigenExpression" with proper spelling for consistency
    plot_data$NeoantigenExpression <- neoantigen_map[plot_data$Sample]
    
    # For mut samples with NA NeoantigenExpression, set to "Low"
    is_mut_na <- plot_data$Type == "mut" & is.na(plot_data$NeoantigenExpression)
    if(any(is_mut_na)) {
        plot_data$NeoantigenExpression[is_mut_na] <- "Low"
    }
    
    # Filter out Low neoantigenexpression IDs
    plot_data <- plot_data %>% 
        filter(!(Type == "mut" & NeoantigenExpression == "Low"))
    
    # Combine High and Medium into a single group called "High_Medium"
    # Make sure we capture all variations regardless of case
    plot_data$NeoantigenExpression[plot_data$NeoantigenExpression %in% c("High", "Medium", "high", "medium")] <- "High_Medium"
    
    # Print summary to verify the consolidation
    cat("NeoantigenExpression categories after combining High/Medium:\n")
    print(table(plot_data$NeoantigenExpression, useNA = "always"))
    
    # Create Type_NeoAntigen column for visualization
    plot_data$Type_NeoAntigen <- as.character(plot_data$Type)
    is_mut <- plot_data$Type == "mut"
    plot_data$Type_NeoAntigen[is_mut] <- paste0("mut_", plot_data$NeoantigenExpression[is_mut])
    
    # Make sure Type_NeoAntigen is a factor with the correct order
    plot_data$Type_NeoAntigen <- factor(
        plot_data$Type_NeoAntigen,
        levels = c("wild", "mut_High_Medium")
    )
    
    # Print summary of Type_NeoAntigen categories
    cat("Type_NeoAntigen categories:\n")
    print(table(plot_data$Type_NeoAntigen))
}

# Set factor levels for other variables
plot_data$Gene <- factor(plot_data$Gene, levels = target_genes)
plot_data$Stage <- factor(plot_data$Stage, levels = c("II", "III", "IV"))
plot_data$Type <- factor(plot_data$Type)

# Create color scheme
type_colors <- brewer.pal(n = 3, name = "Set2")[1:2]

# Function to format p-value
format_pvalue <- function(p) {
    sapply(p, function(x) {
        if (x < 0.001) return("p < 0.001")
        if (x < 0.01) return(paste0("p = ", sprintf("%.3f", x)))
        return(paste0("p = ", sprintf("%.2f", x)))
    })
}

# Calculate Kruskal-Wallis test results
kw_results <- plot_data %>%
    group_by(Gene) %>%
    summarise(
        KW_pvalue = kruskal.test(Expression ~ interaction(Type, Stage))$p.value,
        .groups = "drop"
    ) %>%
    mutate(
        KW_significant = KW_pvalue < 0.05,
        KW_label = format_pvalue(KW_pvalue),
        Significance = case_when(
            KW_pvalue < 0.001 ~ "***",
            KW_pvalue < 0.01 ~ "**",
            KW_pvalue < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    )

# Function to create whisker plot for a single gene
create_gene_plot <- function(data, gene, kw_results, show_legend = FALSE) {
    # Get KW test results for this gene
    kw_result <- kw_results %>% filter(Gene == gene)
    
    # Filter data for this gene
    gene_data <- data %>% filter(Gene == gene)
    
    # Create the base plot - using position_dodge for side-by-side boxplots
    if(exists("annotation_with_expression") && "Type_NeoAntigen" %in% colnames(data)) {
        # Use Type_NeoAntigen if available
        p <- ggplot(gene_data, 
               aes(x = Stage, y = Expression, fill = Type_NeoAntigen))
        
        # Define color palette for Type_NeoAntigen (only two colors needed now)
        # Changed color for mut_High_Medium from #E31A1C to #FF9999 (light red)
        neo_colors <- c(
            "wild" = "#66C2A5",        # green for wild
            "mut_High_Medium" = "#FF9999" # light red for combined High_Medium
        )
        
        # Add debugging message to verify the levels
        cat("Type_NeoAntigen levels in plot function:\n")
        print(levels(gene_data$Type_NeoAntigen))
        
        p <- p + scale_fill_manual(values = neo_colors,
                              name = "Type NeoAntigen")
    } else {
        # Use Type if Type_NeoAntigen is not available
        p <- ggplot(gene_data, 
               aes(x = Stage, y = Expression, fill = Type))
               
        p <- p + scale_fill_manual(values = type_colors,
                              name = "Type")
    }
    
    # Complete the plot with side-by-side boxplots
    p <- p + geom_boxplot(position = position_dodge(width = 0.8),
                   outlier.shape = NA) +
        geom_point(position = position_jitterdodge(dodge.width = 0.8),
                   alpha = 0.4,
                   size = 1,
                   shape = 16) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 10),
            legend.position = if(show_legend) "right" else "none",
            panel.grid.minor = element_blank()
        ) +
        labs(
            title = gene,
            subtitle = paste0(
                "Kruskal-Wallis test: ",
                kw_result$KW_label,
                " ",
                kw_result$Significance
            ),
            x = "Stage",
            y = "Expression Level"
        )
    
    return(p)
}

# Check if we have any data after all the filtering
if(nrow(plot_data) == 0) {
    stop("No data remains after filtering. Please check your filtering criteria.")
}

# Print data summary before creating plots
cat("Summary of filtered data:\n")
cat("Total samples:", length(unique(plot_data$Sample)), "\n")
cat("Samples by Type and Stage:\n")
print(table(Type = plot_data$Type, Stage = plot_data$Stage))

# Generate individual plots for the grid
gene_plots <- lapply(target_genes, function(gene) {
    gene_data <- plot_data %>% filter(Gene == gene)
    if(nrow(gene_data) == 0) {
        warning("No data available for gene: ", gene)
        # Return empty plot with a message
        return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = paste("No data for", gene)) +
               theme_void())
    }
    create_gene_plot(plot_data, gene, kw_results, show_legend = FALSE)
})

# Get the legend from a temporary plot
legend_plot <- create_gene_plot(plot_data, target_genes[1], kw_results, show_legend = TRUE)
legend <- cowplot::get_legend(legend_plot)

# Create significance text
sig_text <- textGrob(
    "Significance levels: *** p<0.001, ** p<0.01, * p<0.05, ns: not significant",
    gp = gpar(fontsize = 10)
)

# Arrange plots in a grid
plots_arranged <- do.call(gridExtra::arrangeGrob, c(gene_plots, ncol = 3))

# Create bottom panel with legend and significance
bottom <- gridExtra::arrangeGrob(
    legend,
    sig_text,
    ncol = 1,
    heights = c(1, 1)
)

# Combine everything
combined_plot <- gridExtra::grid.arrange(
    plots_arranged,
    bottom,
    heights = c(10, 2),
    ncol = 1
)

# Save combined plot in both PDF and PNG formats
ggsave(
    filename = file.path("results_geneset_modified", "all_genes_whisker_plots_modified.pdf"),
    plot = combined_plot,
    width = 12,
    height = 15,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

ggsave(
    filename = file.path("results_geneset_modified", "all_genes_whisker_plots_modified.png"),
    plot = combined_plot,
    width = 12,
    height = 15,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

# Also save SVG format for better scalability
ggsave(
    filename = file.path("results_geneset_modified", "all_genes_whisker_plots_modified.svg"),
    plot = combined_plot,
    width = 12,
    height = 15,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

# Load required packages
library(survival)
library(survminer)
library(ggplot2)

# Function to read the specialized format
read_survival_data <- function(file_path) {
    # Read all lines from the file
    lines <- readLines(file_path)
    
    # Initialize empty data frame for results
    all_data <- data.frame()
    
    # Variables to track current group
    current_group <- NULL
    header_line <- NULL
    
    # Process each line
    for (i in 1:length(lines)) {
        line <- trimws(lines[i])
        
        # Skip empty lines
        if (line == "") next
        
        # Check if this is a group header
        if (grepl("^\\([A-Z]\\)", line)) {
            # Extract group name
            current_group <- gsub("^\\([A-Z]\\)\\s+", "", line)
        } 
        # Check if this is a header line
        else if (grepl("^Case ID", line) && !is.null(current_group)) {
            header_line <- unlist(strsplit(line, "\t"))
        } 
        # This is a data line
        else if (!is.null(current_group) && !is.null(header_line) && grepl("\t", line)) {
            # Parse data line
            values <- unlist(strsplit(line, "\t"))
            
            # Create row data frame
            if (length(values) == length(header_line)) {
                row_data <- as.data.frame(t(values), stringsAsFactors = FALSE)
                names(row_data) <- header_line
                
                # Add group information
                row_data$Group <- current_group
                
                # Append to result data frame
                all_data <- rbind(all_data, row_data)
            }
        }
    }
    
    # Convert numeric columns
    numeric_cols <- c("Number at Risk", "Survival Rate", "Time (months)")
    for (col in numeric_cols) {
        if (col %in% names(all_data)) {
            all_data[[col]] <- as.numeric(all_data[[col]])
        }
    }
    
    # Convert status to binary event indicator (1 = event, 0 = censored)
    all_data$Event <- ifelse(all_data$Status == "deceased", 1, 0)
    
    return(all_data)
}

# Read the data
data <- read_survival_data("overallSurvival data/pole_mut_stage2and3_cd8_overall.txt")

# Get the max time value
max_time <- max(data$`Time (months)`)
print(paste("Maximum time value:", max_time))

# Create custom breaks - regular intervals up to 240, then the max value
custom_breaks <- c(seq(0, 240, by = 24), max_time)

# Create survival object
surv_obj <- Surv(time = data$`Time (months)`, event = data$Event)

# Fit survival curve by group
fit <- survfit(surv_obj ~ Group, data = data)

# Create the survival plot with custom x-axis
survplot <- ggsurvplot(
    fit,
    data = data,
    risk.table = FALSE,
    pval = TRUE,
    conf.int = FALSE,
    xlim = c(0, max_time),
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    title = "Overall Survival",
    palette = c("#4285F4", "#EA4335"),
    legend.title = "Group",
    legend.labs = c("pole_mut_stage2and3_cd8_high", "pole_mut_stage2and3_cd8_low"),
    ggtheme = theme_bw(),
    font.main = 16,
    font.x = 14,
    font.y = 14,
    font.tickslab = 12,
    font.legend = 12,
    break.time.by = NULL  # Remove default breaks
)

# Modify the x-axis breaks using the survplot$plot object
survplot$plot <- survplot$plot + 
    scale_x_continuous(breaks = custom_breaks)

# Save as PNG (high resolution)
ggsave(
    filename = "overallSurvival data/pole_mut_stage2and3_cd8_overall.png", 
    plot = print(survplot), 
    width = 10, 
    height = 8, 
    dpi = 300, 
    bg = "white"
)

# Save as SVG
ggsave(
    filename = "overallSurvival data/pole_mut_stage2and3_cd8_overall.svg", 
    plot = print(survplot), 
    width = 10, 
    height = 8, 
    device = "svg", 
    bg = "white"
)

# Gene Expression Comparison (Mutation CD8 High vs Low)
# Script to generate boxplots showing gene expression by stage and CD8 stratification

# Load required libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Define the gene groups
gene_groups <- list(
  original = c("CGAS", "TMEM173", "MAVS", "IFI27", "IFI44L", "IFIT1", "ISG15", "RSAD2", "SIGLEC1"),
  t_cell_markers = c("CD8A", "CD3E", "CD28"),
  naive_memory = c("CCR7", "IL7R", "TCF7", "SELL", "SATB1", "GPR183", "LTB", "LEF1", "S100A10"),
  cytotoxic = c("PRF1", "GZMA", "GZMK", "NKG7"),
  exhausted = c("CXCL13", "HSPB1", "IRF4", "LAYN", "GIMAP6", "HSPH1", "CXCR6", "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT")
)

# Group names for plot titles
group_names <- list(
  original = "ISG Genes",
  t_cell_markers = "T-cell Marker Genes",
  naive_memory = "Naive/Memory CD8+ T Cells",
  cytotoxic = "Cytotoxic CD8+ T Cells",
  exhausted = "Exhausted CD8+ T Cells"
)

# Create output directory
dir.create("gene_plots", showWarnings = FALSE)

# Calculate the maximum number of genes across all groups
max_genes <- max(sapply(gene_groups, length))
# Set standard dimensions for all plots
std_num_cols <- min(3, max_genes)  # Maximum 3 columns
std_num_rows <- ceiling(max_genes / std_num_cols)
std_plot_height <- 8 + (std_num_rows * 3)
std_plot_width <- 12

cat("Using standardized plot dimensions for all gene groups:\n")
cat("Maximum genes in any group:", max_genes, "\n")
cat("Standard plot width:", std_plot_width, "inches\n")
cat("Standard plot height:", std_plot_height, "inches\n")

# Function to prepare data and generate plots
generate_expression_plots <- function(count, annotation_with_neoexpression, cibersortx, gene_groups, group_names) {
  # Ensure annotation data is in the right format
  if(!is.data.frame(annotation_with_neoexpression)) {
    ann_df <- as.data.frame(annotation_with_neoexpression)
  } else {
    ann_df <- annotation_with_neoexpression
  }
  
  # Ensure sample IDs are available
  if(!is.null(rownames(ann_df)) && all(rownames(ann_df) != seq_len(nrow(ann_df)))) {
    ann_df$Mixture <- rownames(ann_df)
  }
  
  # Filter to include only stages II and III
  stages_annotation <- ann_df %>%
    filter(Stage %in% c("II", "III"))
  
  # Prepare the count matrix
  expr_matrix <- count
  
  # Handle rownames in count matrix
  if(!identical(rownames(expr_matrix), as.character(1:nrow(expr_matrix)))) {
    # Matrix already has rownames
    cat("Count matrix already has rownames\n")
  } else {
    # Need to set rownames from first column
    row_names_col <- names(expr_matrix)[1]
    expr_matrix <- as.data.frame(expr_matrix)
    rownames(expr_matrix) <- expr_matrix[[row_names_col]]
    expr_matrix <- expr_matrix[, -1, drop = FALSE]
  }
  
  # Process CIBERSORTx data for CD8 T cell stratification
  cibersortx <- as.data.frame(cibersortx)
  
  # Extract CD8 T cell values
  cd8_data <- cibersortx %>%
    select(Mixture, `T cells CD8`) %>%
    rename(CD8_level = `T cells CD8`)
  
  # Merge CD8 data with annotation
  combined_annotation <- stages_annotation %>%
    left_join(cd8_data, by = "Mixture")
  
  # Calculate maximum CD8 in wild type samples AFTER removing outliers
  wild_cd8 <- combined_annotation %>%
    filter(Type == "wild") %>%
    pull(CD8_level)
  
  # Remove outliers using IQR method
  Q1 <- quantile(wild_cd8, 0.25, na.rm = TRUE)
  Q3 <- quantile(wild_cd8, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # Filter out outliers and get max
  wild_cd8_no_outliers <- wild_cd8[wild_cd8 >= lower_bound & wild_cd8 <= upper_bound]
  max_cd8_wild <- max(wild_cd8_no_outliers, na.rm = TRUE)
  
  cat("Wild type CD8 statistics:\n")
  cat("  Original range:", min(wild_cd8, na.rm = TRUE), "to", max(wild_cd8, na.rm = TRUE), "\n")
  cat("  After outlier removal:", min(wild_cd8_no_outliers, na.rm = TRUE), "to", max_cd8_wild, "\n")
  cat("  Cutoff value for mutation stratification:", max_cd8_wild, "\n")
  
  # Create new type classification with CD8 stratification
  combined_annotation <- combined_annotation %>%
    mutate(Type_CD8 = case_when(
      Type == "mut" & CD8_level > max_cd8_wild ~ "mut_CD8_high",
      Type == "mut" & CD8_level <= max_cd8_wild ~ "mut_CD8_low",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Type_CD8))  # This will filter out wild type samples
  
  # Find matching samples
  matching_samples <- colnames(expr_matrix)[colnames(expr_matrix) %in% combined_annotation$Mixture]
  
  # Filter expression matrix to keep only selected stages and samples with CD8 data
  expr_matrix_filtered <- expr_matrix[, matching_samples]
  
  # Create sample info with new type classification and stage
  sample_info <- combined_annotation %>%
    filter(Mixture %in% matching_samples) %>%
    rename(sample_id = Mixture) %>%
    arrange(match(sample_id, colnames(expr_matrix_filtered)))
  
  # Ensure sample order matches
  expr_matrix_filtered <- expr_matrix_filtered[, sample_info$sample_id]
  
  # Calculate sample counts for each group
  sample_counts <- sample_info %>%
    group_by(Stage, Type_CD8) %>%
    summarise(count = n(), .groups = "drop")
  
  cat("Samples by stage and CD8-stratified type (mut only):\n")
  print(sample_counts)
  
  # Create DGEList and normalize
  dge <- DGEList(counts = expr_matrix_filtered)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Convert to log-CPM
  log_cpm <- cpm(dge, log = TRUE)
  
  # Process each gene group separately
  results_list <- list()
  
  for(group_id in names(gene_groups)) {
    target_genes <- gene_groups[[group_id]]
    group_title <- group_names[[group_id]]
    
    cat("Processing gene group:", group_title, "\n")
    
    # Check if all target genes are present in the data
    missing_genes <- setdiff(target_genes, rownames(log_cpm))
    if(length(missing_genes) > 0) {
      cat("Warning: The following genes were not found in the dataset:", paste(missing_genes, collapse = ", "), "\n")
      # Filter to only include genes that are present in the dataset
      target_genes <- intersect(target_genes, rownames(log_cpm))
    }
    
    if(length(target_genes) == 0) {
      cat("Error: No genes from this group were found in the dataset. Skipping group.\n")
      next
    }
    
    # Extract expression data for target genes
    target_expr <- log_cpm[target_genes, ]
    
    # Prepare data for plotting
    plot_data <- data.frame()
    
    for (gene in target_genes) {
      gene_expr <- target_expr[gene, ]
      gene_df <- data.frame(
        Gene = gene,
        Expression = gene_expr,
        Sample = colnames(target_expr),
        Type_CD8 = sample_info$Type_CD8[match(colnames(target_expr), sample_info$sample_id)],
        Stage = sample_info$Stage[match(colnames(target_expr), sample_info$sample_id)]
      )
      plot_data <- rbind(plot_data, gene_df)
    }
    
    # Factor the variables for consistent ordering
    plot_data$Type_CD8 <- factor(plot_data$Type_CD8, levels = c("mut_CD8_low", "mut_CD8_high"))
    plot_data$Stage <- factor(plot_data$Stage, levels = c("II", "III"))
    
    # Calculate fold changes between CD8 high and low for each gene and stage
    fold_change_results <- plot_data %>%
      group_by(Gene, Stage, Type_CD8) %>%
      summarise(Mean_Expr = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(
        names_from = Type_CD8,
        values_from = Mean_Expr
      ) %>%
      mutate(
        # Calculate fold change (since we're in log space, subtract rather than divide)
        FoldChange = mut_CD8_high - mut_CD8_low,
        # Calculate fold change in non-log scale
        FoldChange_raw = 2^FoldChange
      )
    
    # Perform statistical tests - Wilcoxon rank-sum test for each gene and stage
    stat_results <- data.frame()
    
    for (gene_name in unique(plot_data$Gene)) {
      for (stage_name in unique(plot_data$Stage)) {
        gene_stage_data <- plot_data %>% 
          filter(Gene == gene_name, Stage == stage_name)
        
        if (all(c("mut_CD8_low", "mut_CD8_high") %in% gene_stage_data$Type_CD8)) {
          wilcox_res <- wilcox.test(
            Expression ~ Type_CD8, 
            data = gene_stage_data,
            exact = FALSE
          )
          
          stat_results <- rbind(
            stat_results,
            data.frame(
              Gene = gene_name,
              Stage = stage_name,
              p_value = wilcox_res$p.value
            )
          )
        }
      }
    }
    
    # Add significance symbols
    stat_results <- stat_results %>%
      mutate(
        significance = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**", 
          p_value < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        p_label = sapply(p_value, function(p) {
          if (p < 0.001) return("p < 0.001")
          else if (p < 0.01) return(paste0("p = ", format(round(p, 3), nsmall = 3)))
          else return(paste0("p = ", format(round(p, 2), nsmall = 2)))
        })
      )
    
    # Create individual plots for each gene
    plot_list <- list()
    
    # Custom color palette for the two groups
    color_palette <- c("mut_CD8_low" = "#E41A1C", "mut_CD8_high" = "#377EB8")
    
    # Using standard grid layout for consistency across all groups
    num_cols <- std_num_cols
    num_rows <- std_num_rows
    
    # Calculate number of empty plots needed to complete the grid
    num_genes <- length(target_genes)
    num_empty_plots <- (num_cols * num_rows) - num_genes
    
    for (i in 1:num_genes) {
      gene <- target_genes[i]
      gene_data <- plot_data %>% filter(Gene == gene)
      
      # Extract stats for this gene
      gene_stats <- stat_results %>% filter(Gene == gene)
      
      # Add fold change info
      gene_fc <- fold_change_results %>% filter(Gene == gene) %>%
        mutate(FC_label = paste0("FC = ", format(round(FoldChange_raw, 2), nsmall = 2)))
      
      # Find max y value for appropriate annotation placement
      max_y <- max(gene_data$Expression, na.rm = TRUE)
      
      # Create plot
      p <- ggplot(gene_data, aes(x = Stage, y = Expression, fill = Type_CD8)) +
        geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
        geom_point(aes(color = Type_CD8), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.7, size = 0.8) +
        scale_fill_manual(values = color_palette) +
        scale_color_manual(values = color_palette) +
        labs(
          title = gene,
          x = "Stage",
          y = "Expression Level (log-CPM)"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",  # No individual legends
          panel.grid.minor = element_blank()
        )
      
      # Add p-value and fold change annotations
      for(stage_val in c("II", "III")) {
        stage_stat <- gene_stats %>% filter(Stage == stage_val)
        stage_fc <- gene_fc %>% filter(Stage == stage_val)
        
        if(nrow(stage_stat) > 0) {
          # Add p-value
          p <- p + annotate(
            "text",
            x = which(levels(plot_data$Stage) == stage_val),
            y = max_y + 1.2,
            label = paste0(stage_stat$p_label, " ", stage_stat$significance),
            size = 2.5
          )
          
          # Add fold change
          if(nrow(stage_fc) > 0) {
            p <- p + annotate(
              "text",
              x = which(levels(plot_data$Stage) == stage_val),
              y = max_y + 0.8,
              label = stage_fc$FC_label,
              size = 2.5
            )
          }
        }
      }
      
      plot_list[[gene]] <- p
    }
    
    # Add empty plots if needed to complete the grid
    if (num_empty_plots > 0) {
      for (i in 1:num_empty_plots) {
        # Create an empty plot with the same theme
        empty_plot <- ggplot() + 
          theme_void() +
          theme(
            plot.background = element_rect(fill = "white", color = NA)
          )
        plot_list[[paste0("empty_", i)]] <- empty_plot
      }
    }
    
    # Arrange plots in a grid
    combined_plot <- plot_grid(
      plotlist = plot_list,
      ncol = num_cols,
      nrow = num_rows
    )
    
    # Create a manual color legend panel
    color_legend_panel <- ggdraw() +
      draw_line(
        x = c(0.35, 0.40), 
        y = c(0.5, 0.5), 
        color = color_palette["mut_CD8_low"], 
        size = 3
      ) +
      draw_line(
        x = c(0.60, 0.65), 
        y = c(0.5, 0.5), 
        color = color_palette["mut_CD8_high"], 
        size = 3
      ) +
      draw_label(
        "CD8 Low", 
        x = 0.30, 
        y = 0.5, 
        hjust = 1, 
        size = 12
      ) +
      draw_label(
        "CD8 High", 
        x = 0.55, 
        y = 0.5, 
        hjust = 1, 
        size = 12
      )
    
    # Create sample size text for all groups
    sample_size_text <- sample_counts %>%
      mutate(
        type_label = case_when(
          Type_CD8 == "mut_CD8_low" ~ "CD8 Low",
          Type_CD8 == "mut_CD8_high" ~ "CD8 High"
        ),
        group_label = paste0("Stage ", Stage, " ", type_label, ": n=", count)
      ) %>%
      pull(group_label) %>%
      paste(collapse = ", ")
    
    # Add significance legend
    significance_text <- ggdraw() + 
      draw_label(
        "Significance: *** p<0.001, ** p<0.01, * p<0.05, ns: not significant",
        x = 0.5, y = 0.7, hjust = 0.5, vjust = 0.5, size = 10
      )
    
    # Add sample size information under significance legend
    sample_text <- ggdraw() + 
      draw_label(
        paste("Sample sizes:", sample_size_text),
        x = 0.5, y = 0.3, hjust = 0.5, vjust = 0.5, size = 10
      )
    
    # Combine significance and sample size into one text area
    info_text <- plot_grid(
      significance_text,
      sample_text,
      ncol = 1,
      rel_heights = c(1, 1)
    )
    
    # Add group title
    group_title_text <- ggdraw() + 
      draw_label(
        group_title,
        x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5, size = 16, fontface = "bold"
      )
    
    # Add CD8 cutoff and fold change reference information
    cutoff_text <- ggdraw() + 
      draw_label(
        paste("CD8 T cell cutoff value (wild type max after outlier removal):", round(max_cd8_wild, 4),
              "| Fold change (FC) = CD8 High / CD8 Low"),
        x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5, size = 10
      )
    
    # Combine everything
    final_plot <- plot_grid(
      group_title_text,  # Group title
      combined_plot,     # Gene plots
      color_legend_panel, # Color legend
      plot_grid(info_text, cutoff_text, ncol = 1, rel_heights = c(2, 1)), # Info text
      ncol = 1,
      rel_heights = c(1, 20, 1, 3)
    )
    
    # Use standardized dimensions for all plots
    # Save the plot in SVG format
    svglite::svglite(
      file.path("gene_plots", paste0("gene_expression_", group_id, ".svg")),
      width = std_plot_width,
      height = std_plot_height
    )
    print(final_plot)
    dev.off()
    
    # Also save PNG for easy viewing
    ggsave(
      file.path("gene_plots", paste0("gene_expression_", group_id, ".png")),
      final_plot,
      width = std_plot_width,
      height = std_plot_height,
      dpi = 300
    )
    
    cat("Plot for group '", group_title, "' saved to gene_plots directory (", std_plot_width, "x", std_plot_height, " inches).\n", sep="")
    
    # Create fold change summary table for this gene group
    fc_summary <- fold_change_results %>%
      left_join(stat_results, by = c("Gene", "Stage")) %>%
      mutate(
        FoldChange_raw = round(FoldChange_raw, 2),
        p_value = format(p_value, digits = 3),
        significant = p_value < 0.05
      ) %>%
      select(Gene, Stage, FoldChange_raw, p_value, significance)
    
    # Save the fold change summary for this group
    write.csv(
      fc_summary, 
      file.path("gene_plots", paste0("fold_change_", group_id, ".csv")), 
      row.names = FALSE
    )
    
    cat("Fold change summary for '", group_title, "' saved.\n", sep="")
    
    # Store results for this group
    results_list[[group_id]] <- list(
      plot_data = plot_data,
      fold_change = fold_change_results,
      stat_results = stat_results
    )
  }
  
  # Return all results
  return(list(
    group_results = results_list,
    cd8_cutoff = max_cd8_wild,
    sample_counts = sample_counts
  ))
}

# Run the function to generate the plots
# Make sure svglite package is installed first
if(!require(svglite)) {
  install.packages("svglite")
  library(svglite)
}

results <- generate_expression_plots(
  count, 
  annotation_with_neoexpression, 
  cibersortx, 
  gene_groups, 
  group_names
)

# Display summary information
cat("\nCD8 T cell cutoff value (wild type max after outlier removal):", round(results$cd8_cutoff, 4), "\n")
cat("\nSample counts by group:\n")
print(results$sample_counts)

# Display fold change and p-value summaries for each gene group
for(group_id in names(gene_groups)) {
  if(group_id %in% names(results$group_results)) {
    group_title <- group_names[[group_id]]
    cat("\n========================================\n")
    cat("Summary for gene group:", group_title, "\n")
    cat("========================================\n")
    
    group_results <- results$group_results[[group_id]]
    
    cat("\nFold change summary (CD8 High vs CD8 Low):\n")
    print(group_results$fold_change %>% select(Gene, Stage, FoldChange_raw))
    cat("\nStatistical test results:\n")
    print(group_results$stat_results %>% select(Gene, Stage, p_value, significance))
  }
}


# R script for cell population analysis with SVG Venn diagrams
# Load required libraries
library(tidyverse)
library(VennDiagram)  # For proper Venn diagrams
library(grid)         # Required for grid.draw function

# Define cell types and stages to analyze
cell_types <- c("Macrophages M1", "T cells CD8")
stages <- c("II", "III")

# Define consistent colors for cell types
m1_color <- "steelblue"  # Color for Macrophages M1
cd8_color <- "darkred"   # Color for T cells CD8

# Step 1: Merge the datasets and fix the Type column
merged_data <- annotation_clean %>%
    mutate(Type = as.character(Type)) %>%
    inner_join(cibersortx, by = c("SampleID" = "Mixture"))

# Step 2: Process each stage
for(stage in stages) {
    cat("\n\n========= ANALYSIS FOR STAGE", stage, "=========\n")
    
    # Filter for current stage
    stage_data <- merged_data %>% filter(Stage == stage)
    
    # Process each cell type  
    results_list <- list()
    
    for(cell_type in cell_types) {
        cat("\n--- Analysis for", cell_type, "---\n")
        
        # Step 3: Separate wild type and mutant data
        wt_data <- stage_data %>% filter(Type == "wild")
        mut_data <- stage_data %>% filter(Type == "mut")
        
        cat("Number of wild type samples:", nrow(wt_data), "\n")
        cat("Number of mutant samples:", nrow(mut_data), "\n")
        
        # Step 4: Remove outliers from wild type data
        if(nrow(wt_data) > 0) {
            q1_wt <- quantile(wt_data[[cell_type]], 0.25, na.rm = TRUE)
            q3_wt <- quantile(wt_data[[cell_type]], 0.75, na.rm = TRUE)
            iqr_wt <- q3_wt - q1_wt
            
            lower_bound_wt <- q1_wt - 1.5 * iqr_wt
            upper_bound_wt <- q3_wt + 1.5 * iqr_wt
            
            wt_filtered <- wt_data %>% 
                filter(!!sym(cell_type) >= lower_bound_wt & !!sym(cell_type) <= upper_bound_wt)
            
            cat("Wild type samples after outlier removal:", nrow(wt_filtered), "\n")
            
            # Step 5: Use highest value from wild type as cutoff
            cutoff <- max(wt_filtered[[cell_type]])
            cat("Cutoff value (highest wild type value after outlier removal):", cutoff, "\n")
        } else {
            # If no wild type data, use highest value of mutant data after outlier removal
            q1_mut <- quantile(mut_data[[cell_type]], 0.25, na.rm = TRUE)
            q3_mut <- quantile(mut_data[[cell_type]], 0.75, na.rm = TRUE)
            iqr_mut <- q3_mut - q1_mut
            
            lower_bound_mut <- q1_mut - 1.5 * iqr_mut
            upper_bound_mut <- q3_mut + 1.5 * iqr_mut
            
            mut_filtered_for_cutoff <- mut_data %>% 
                filter(!!sym(cell_type) >= lower_bound_mut & !!sym(cell_type) <= upper_bound_mut)
            
            cutoff <- max(mut_filtered_for_cutoff[[cell_type]])
            cat("Cutoff value (highest mutant value after outlier removal):", cutoff, "\n")
        }
        
        # Step 6: Remove outliers from mutant data
        if(nrow(mut_data) > 0) {
            q1_mut <- quantile(mut_data[[cell_type]], 0.25, na.rm = TRUE)
            q3_mut <- quantile(mut_data[[cell_type]], 0.75, na.rm = TRUE)
            iqr_mut <- q3_mut - q1_mut
            
            lower_bound_mut <- q1_mut - 1.5 * iqr_mut
            upper_bound_mut <- q3_mut + 1.5 * iqr_mut
            
            mut_filtered <- mut_data %>% 
                filter(!!sym(cell_type) >= lower_bound_mut & !!sym(cell_type) <= upper_bound_mut)
            
            cat("Mutant samples after outlier removal:", nrow(mut_filtered), "\n")
            
            # Step 7: Classify mutant patients as high or low based on cutoff
            mut_classified <- mut_filtered %>%
                mutate(Classification = ifelse(!!sym(cell_type) > cutoff, "High", "Low"))
            
            # Step 8: Extract high and low patient lists
            high_patients <- mut_classified %>%
                filter(Classification == "High") %>%
                pull(SampleID)
            
            low_patients <- mut_classified %>%
                filter(Classification == "Low") %>%
                pull(SampleID)
            
            # Step 9: Generate and print summary statistics
            summary_stats <- mut_classified %>%
                group_by(Classification) %>%
                summarize(
                    Count = n(),
                    Mean = mean(!!sym(cell_type)),
                    Median = median(!!sym(cell_type)),
                    Min = min(!!sym(cell_type)),
                    Max = max(!!sym(cell_type))
                )
            
            print(summary_stats)
            
            # Step 10: Store results for later comparison
            results_list[[cell_type]] <- list(
                filtered_data = mut_classified,
                high_patients = high_patients,
                low_patients = low_patients,
                cutoff = cutoff,
                summary_stats = summary_stats
            )
            
            # Create classification dataframe with proper column naming
            classification_df <- mut_classified %>%
                select(SampleID, !!sym(cell_type))
            
            # Add classification column with a conventional name
            classification_col_name <- paste0(gsub(" ", "_", cell_type), "_Classification")
            classification_df[[classification_col_name]] <- mut_classified$Classification
            
            # Create a new path for the output file
            output_file <- paste0("Stage_", stage, "_", gsub(" ", "_", cell_type), "_classification.csv")
            
            # Check and attempt to remove the file if it exists (to avoid permission errors)
            if (file.exists(output_file)) {
                try(file.remove(output_file), silent = TRUE)
            }
            
            # Join with merged_df to get Patient IDs
            classification_with_patient_id <- classification_df %>%
                inner_join(merged_df, by = c("SampleID" = "WES_id")) %>%
                select(`Patient ID`, `Sample ID`, SampleID, everything())
            
            # Write results to CSV with error handling
            tryCatch({
                write.csv(
                    classification_with_patient_id,
                    output_file,
                    row.names = FALSE
                )
                cat("Successfully saved:", output_file, "\n")
            }, error = function(e) {
                cat("Error saving file:", output_file, "\n")
                cat("Error message:", e$message, "\n")
            })
        }
    }
    
    # Step 11: Analyze overlaps between cell types (if both analyzed)
    if(length(results_list) == 2) {
        # Get high patients from both cell types
        high_m1 <- results_list[[cell_types[1]]]$high_patients
        high_cd8 <- results_list[[cell_types[2]]]$high_patients
        
        # Get low patients from both cell types
        low_m1 <- results_list[[cell_types[1]]]$low_patients
        low_cd8 <- results_list[[cell_types[2]]]$low_patients
        
        # Find overlaps for HIGH group
        high_overlapping <- intersect(high_m1, high_cd8)
        
        # Find overlaps for LOW group
        low_overlapping <- intersect(low_m1, low_cd8)
        
        # Print overlap analysis
        cat("\n===== Overlap Analysis for Stage", stage, "=====\n")
        cat("High", cell_types[1], "patients:", length(high_m1), "\n")
        cat("High", cell_types[2], "patients:", length(high_cd8), "\n")
        cat("Overlapping patients (high in both):", length(high_overlapping), "\n")
        cat("Low", cell_types[1], "patients:", length(low_m1), "\n")
        cat("Low", cell_types[2], "patients:", length(low_cd8), "\n")
        cat("Overlapping patients (low in both):", length(low_overlapping), "\n")
        
        # Step 12: Create comprehensive classification dataframe
        m1_class_df <- results_list[[cell_types[1]]]$filtered_data %>% 
            select(SampleID, !!sym(cell_types[1]))
        m1_class_df[[paste0(gsub(" ", "_", cell_types[1]), "_Classification")]] <- 
            results_list[[cell_types[1]]]$filtered_data$Classification
        
        cd8_class_df <- results_list[[cell_types[2]]]$filtered_data %>% 
            select(SampleID, !!sym(cell_types[2]))
        cd8_class_df[[paste0(gsub(" ", "_", cell_types[2]), "_Classification")]] <- 
            results_list[[cell_types[2]]]$filtered_data$Classification
        
        # Join the classification dataframes
        comprehensive_classification <- full_join(
            m1_class_df,
            cd8_class_df,
            by = "SampleID"
        )
        
        # Output file for comprehensive classification
        comprehensive_file <- paste0("Stage_", stage, "_Comprehensive_Classification.csv")
        
        # Check and attempt to remove the file if it exists
        #if (file.exists(comprehensive_file)) {
        #    try(file.remove(comprehensive_file), silent = TRUE)
        #}
        
        # Join with merged_df to get Patient IDs
        comprehensive_with_patient_id <- comprehensive_classification %>%
            inner_join(merged_df, by = c("SampleID" = "WES_id")) %>%
            select(`Patient ID`, `Sample ID`, SampleID, everything())
        
        # Save comprehensive classification with error handling
        tryCatch({
            write.csv(
                comprehensive_with_patient_id,
                comprehensive_file,
                row.names = FALSE
            )
            cat("Successfully saved:", comprehensive_file, "\n")
        }, error = function(e) {
            cat("Error saving file:", comprehensive_file, "\n")
            cat("Error message:", e$message, "\n")
        })
        
        # Step 13: Create proper Venn diagrams with consistent colors and better sizing
        
        # Prepare SVG files for HIGH and LOW Venn diagrams
        high_venn_file <- paste0("Stage_", stage, "_HIGH_Venn_Diagram.svg")
        low_venn_file <- paste0("Stage_", stage, "_LOW_Venn_Diagram.svg")
        
        # Check and attempt to remove the files if they exist
        if (file.exists(high_venn_file)) try(file.remove(high_venn_file), silent = TRUE)
        if (file.exists(low_venn_file)) try(file.remove(low_venn_file), silent = TRUE)
        
        # Calculate percentages for labels
        high_m1_only <- setdiff(high_m1, high_cd8)
        high_cd8_only <- setdiff(high_cd8, high_m1)
        high_both <- intersect(high_m1, high_cd8)
        
        low_m1_only <- setdiff(low_m1, low_cd8)
        low_cd8_only <- setdiff(low_cd8, low_m1)
        low_both <- intersect(low_m1, low_cd8)
        
        # HIGH group Venn diagram with fixed dimensions and better labels
        # When directly saving to file, venn.diagram() returns NULL or a numeric
        # Set draw=FALSE to get a grid object we can manipulate
        high_venn_obj <- venn.diagram(
            x = list(
                "Macrophages M1 High" = high_m1,
                "T cells CD8 High" = high_cd8
            ),
            filename = NULL,  # Don't save to file yet
            fill = c(m1_color, cd8_color),
            alpha = 0.5,
            
            # Improve text scaling
            cex = 1.2,                  # Size of count/percentage text
            cat.cex = 1.0,              # Size of category names
            cat.pos = c(-20, 20),       # Position categories clearly away from circles
            cat.dist = c(0.05, 0.05),   # Distance of category names from diagrams
            
            # Improve diagram title
            main = paste("Stage", stage, "- HIGH Group Overlap"),
            main.cex = 1.2,
            
            # Add counts and percentages to labels
            cat.just = list(c(0.5, 1), c(0.5, 0)),  # Text justification
            
            # Customize label style
            cat.default.pos = "outer",
            cat.fontfamily = "sans",
            fontfamily = "sans",
            
            # Diagram styling
            lwd = 1.5,
            euler.d = TRUE,
            scaled = TRUE,
            
            # Custom label text
            cat.prompts = FALSE        # Remove default "n=" text
        )
        
        # LOW group Venn diagram - using the SAME COLORS for consistency
        low_venn_obj <- venn.diagram(
            x = list(
                "Macrophages M1 Low" = low_m1,
                "T cells CD8 Low" = low_cd8
            ),
            filename = NULL,  # Don't save to file yet
            fill = c(m1_color, cd8_color),
            alpha = 0.5,
            
            # Improve text scaling
            cex = 1.2,                  # Size of count/percentage text
            cat.cex = 1.0,              # Size of category names
            cat.pos = c(-20, 20),       # Position categories clearly away from circles
            cat.dist = c(0.05, 0.05),   # Distance of category names from diagrams
            
            # Improve diagram title
            main = paste("Stage", stage, "- LOW Group Overlap"),
            main.cex = 1.2,
            
            # Add counts and percentages to labels
            cat.just = list(c(0.5, 1), c(0.5, 0)),  # Text justification
            
            # Customize label style
            cat.default.pos = "outer",
            cat.fontfamily = "sans",
            fontfamily = "sans",
            
            # Diagram styling
            lwd = 1.5,
            euler.d = TRUE,
            scaled = TRUE,
            
            # Custom label text
            cat.prompts = FALSE        # Remove default "n=" text
        )
        
        # Save the grid objects to SVG files
        svg(high_venn_file, width = 7, height = 7)
        grid.draw(high_venn_obj)
        dev.off()
        
        svg(low_venn_file, width = 7, height = 7)
        grid.draw(low_venn_obj)
        dev.off()
        
        # Preview diagrams in R console
        cat("Previewing HIGH group Venn diagram in R console...\n")
        grid.newpage()
        grid.draw(high_venn_obj)
        
        cat("Previewing LOW group Venn diagram in R console...\n")
        grid.newpage()
        grid.draw(low_venn_obj)
        
        cat("Venn diagrams saved as SVG files with improved sizing and readability:\n")
        cat("  -", high_venn_file, "\n")
        cat("  -", low_venn_file, "\n")
        
        # Optionally: Save additional PNG versions for easier viewing
        high_venn_png <- paste0("Stage_", stage, "_HIGH_Venn_Diagram.png")
        low_venn_png <- paste0("Stage_", stage, "_LOW_Venn_Diagram.png")
        
        # Save as PNG with explicit size control
        png(high_venn_png, width = 800, height = 800, res = 120)
        grid.draw(high_venn_obj)
        dev.off()
        
        png(low_venn_png, width = 800, height = 800, res = 120)
        grid.draw(low_venn_obj)
        dev.off()
        
        cat("Also saved as PNG files for easier viewing:\n")
        cat("  -", high_venn_png, "\n")
        cat("  -", low_venn_png, "\n")
    }
} 



#Analysis of I, II and III stages cell population

# R script for cell population analysis
# Load required libraries
library(tidyverse)

# Define cell types and stages to analyze
cell_types <- c("Macrophages M1", "T cells CD8")
stages <- c("I", "II", "III")  # Added stage I as requested

# Define consistent colors for cell types (kept for potential future use)
m1_color <- "steelblue"  # Color for Macrophages M1
cd8_color <- "darkred"   # Color for T cells CD8

# Step 1: Merge the datasets and fix the Type column
merged_data <- annotation_clean %>%
    mutate(Type = as.character(Type)) %>%
    inner_join(cibersortx, by = c("SampleID" = "Mixture"))

# Step 2: Process each stage
for(stage in stages) {
    cat("\n\n========= ANALYSIS FOR STAGE", stage, "=========\n")
    
    # Filter for current stage
    stage_data <- merged_data %>% filter(Stage == stage)
    
    # Process each cell type  
    results_list <- list()
    
    for(cell_type in cell_types) {
        cat("\n--- Analysis for", cell_type, "---\n")
        
        # Step 3: Separate wild type and mutant data
        wt_data <- stage_data %>% filter(Type == "wild")
        mut_data <- stage_data %>% filter(Type == "mut")
        
        cat("Number of wild type samples:", nrow(wt_data), "\n")
        cat("Number of mutant samples:", nrow(mut_data), "\n")
        
        # Step 4: Remove outliers from wild type data
        if(nrow(wt_data) > 0) {
            q1_wt <- quantile(wt_data[[cell_type]], 0.25, na.rm = TRUE)
            q3_wt <- quantile(wt_data[[cell_type]], 0.75, na.rm = TRUE)
            iqr_wt <- q3_wt - q1_wt
            
            lower_bound_wt <- q1_wt - 1.5 * iqr_wt
            upper_bound_wt <- q3_wt + 1.5 * iqr_wt
            
            wt_filtered <- wt_data %>% 
                filter(!!sym(cell_type) >= lower_bound_wt & !!sym(cell_type) <= upper_bound_wt)
            
            cat("Wild type samples after outlier removal:", nrow(wt_filtered), "\n")
            
            # Step 5: Use highest value from wild type as cutoff
            cutoff <- max(wt_filtered[[cell_type]])
            cat("Cutoff value (highest wild type value after outlier removal):", cutoff, "\n")
        } else {
            # If no wild type data, use highest value of mutant data after outlier removal
            q1_mut <- quantile(mut_data[[cell_type]], 0.25, na.rm = TRUE)
            q3_mut <- quantile(mut_data[[cell_type]], 0.75, na.rm = TRUE)
            iqr_mut <- q3_mut - q1_mut
            
            lower_bound_mut <- q1_mut - 1.5 * iqr_mut
            upper_bound_mut <- q3_mut + 1.5 * iqr_mut
            
            mut_filtered_for_cutoff <- mut_data %>% 
                filter(!!sym(cell_type) >= lower_bound_mut & !!sym(cell_type) <= upper_bound_mut)
            
            cutoff <- max(mut_filtered_for_cutoff[[cell_type]])
            cat("Cutoff value (highest mutant value after outlier removal):", cutoff, "\n")
        }
        
        # Step 6: Remove outliers from mutant data
        if(nrow(mut_data) > 0) {
            q1_mut <- quantile(mut_data[[cell_type]], 0.25, na.rm = TRUE)
            q3_mut <- quantile(mut_data[[cell_type]], 0.75, na.rm = TRUE)
            iqr_mut <- q3_mut - q1_mut
            
            lower_bound_mut <- q1_mut - 1.5 * iqr_mut
            upper_bound_mut <- q3_mut + 1.5 * iqr_mut
            
            mut_filtered <- mut_data %>% 
                filter(!!sym(cell_type) >= lower_bound_mut & !!sym(cell_type) <= upper_bound_mut)
            
            cat("Mutant samples after outlier removal:", nrow(mut_filtered), "\n")
            
            # Step 7: Classify mutant patients as high or low based on cutoff
            mut_classified <- mut_filtered %>%
                mutate(Classification = ifelse(!!sym(cell_type) > cutoff, "High", "Low"))
            
            # Step 7a: All wild type samples (after outlier removal) are classified as "Low"
            wt_classified <- NULL
            if(nrow(wt_filtered) > 0) {
                wt_classified <- wt_filtered %>%
                    mutate(Classification = "Low")  # All wild type samples are "Low"
            }
            
            # Step 7b: Combine mutant and wild type classifications
            all_classified <- bind_rows(mut_classified, wt_classified)
            
            # Step 8: Extract high and low patient lists
            high_patients <- all_classified %>%
                filter(Classification == "High") %>%
                pull(SampleID)
            
            low_patients <- all_classified %>%
                filter(Classification == "Low") %>%
                pull(SampleID)
            
            # Step 9: Generate and print summary statistics
            summary_stats <- all_classified %>%
                group_by(Classification) %>%
                summarize(
                    Count = n(),
                    Mean = mean(!!sym(cell_type)),
                    Median = median(!!sym(cell_type)),
                    Min = min(!!sym(cell_type)),
                    Max = max(!!sym(cell_type))
                )
            
            print(summary_stats)
            
            # Step 10: Store results for later comparison
            results_list[[cell_type]] <- list(
                filtered_data = all_classified,
                high_patients = high_patients,
                low_patients = low_patients,
                cutoff = cutoff,
                summary_stats = summary_stats
            )
            
            # Create classification dataframe with proper column naming
            classification_df <- all_classified %>%
                select(SampleID, !!sym(cell_type))
            
            # Add classification column with a conventional name
            classification_col_name <- paste0(gsub(" ", "_", cell_type), "_Classification")
            classification_df[[classification_col_name]] <- all_classified$Classification
            
            # Add Type column to identify wild vs mutant samples
            classification_df$Type <- all_classified$Type
            
            # Create a new path for the output file
            output_file <- paste0("Stage_", stage, "_", gsub(" ", "_", cell_type), "_classification.csv")
            
            # Check and attempt to remove the file if it exists (to avoid permission errors)
            if (file.exists(output_file)) {
                try(file.remove(output_file), silent = TRUE)
            }
            
            # Join with merged_df to get Patient IDs
            classification_with_patient_id <- classification_df %>%
                inner_join(merged_df, by = c("SampleID" = "WES_id")) %>%
                select(`Patient ID`, `Sample ID`, SampleID, everything())
            
            # Write results to CSV with error handling
            tryCatch({
                write.csv(
                    classification_with_patient_id,
                    output_file,
                    row.names = FALSE
                )
                cat("Successfully saved:", output_file, "\n")
            }, error = function(e) {
                cat("Error saving file:", output_file, "\n")
                cat("Error message:", e$message, "\n")
            })
        }
    }
    
    # Step 11: Analyze overlaps between cell types (if both analyzed)
    if(length(results_list) == 2) {
        # Get high patients from both cell types
        high_m1 <- results_list[[cell_types[1]]]$high_patients
        high_cd8 <- results_list[[cell_types[2]]]$high_patients
        
        # Get low patients from both cell types
        low_m1 <- results_list[[cell_types[1]]]$low_patients
        low_cd8 <- results_list[[cell_types[2]]]$low_patients
        
        # Find overlaps for HIGH group
        high_overlapping <- intersect(high_m1, high_cd8)
        
        # Find overlaps for LOW group
        low_overlapping <- intersect(low_m1, low_cd8)
        
        # Print overlap analysis
        cat("\n===== Overlap Analysis for Stage", stage, "=====\n")
        cat("High", cell_types[1], "patients:", length(high_m1), "\n")
        cat("High", cell_types[2], "patients:", length(high_cd8), "\n")
        cat("Overlapping patients (high in both):", length(high_overlapping), "\n")
        cat("Low", cell_types[1], "patients:", length(low_m1), "\n")
        cat("Low", cell_types[2], "patients:", length(low_cd8), "\n")
        cat("Overlapping patients (low in both):", length(low_overlapping), "\n")
        
        # Step 12: Create comprehensive classification dataframe
        m1_class_df <- results_list[[cell_types[1]]]$filtered_data %>% 
            select(SampleID, !!sym(cell_types[1]))
        m1_class_df[[paste0(gsub(" ", "_", cell_types[1]), "_Classification")]] <- 
            results_list[[cell_types[1]]]$filtered_data$Classification
        
        cd8_class_df <- results_list[[cell_types[2]]]$filtered_data %>% 
            select(SampleID, !!sym(cell_types[2]))
        cd8_class_df[[paste0(gsub(" ", "_", cell_types[2]), "_Classification")]] <- 
            results_list[[cell_types[2]]]$filtered_data$Classification
        
        # Join the classification dataframes
        comprehensive_classification <- full_join(
            m1_class_df,
            cd8_class_df,
            by = "SampleID"
        )
        
        # Output file for comprehensive classification
        comprehensive_file <- paste0("Stage_", stage, "_Comprehensive_Classification.csv")
        
        # Join with merged_df to get Patient IDs
        comprehensive_with_patient_id <- comprehensive_classification %>%
            inner_join(merged_df, by = c("SampleID" = "WES_id")) %>%
            select(`Patient ID`, `Sample ID`, SampleID, everything())
        
        # Save comprehensive classification with error handling
        tryCatch({
            write.csv(
                comprehensive_with_patient_id,
                comprehensive_file,
                row.names = FALSE
            )
            cat("Successfully saved:", comprehensive_file, "\n")
        }, error = function(e) {
            cat("Error saving file:", comprehensive_file, "\n")
            cat("Error message:", e$message, "\n")
        })
    }
}
