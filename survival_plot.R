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
