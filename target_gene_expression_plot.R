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
