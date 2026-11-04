# =====================================================================
# SCRIPT: 01_eDNA_Functional_Analysis.R
# PROJECT: Urban Forest Resilience & Functional Homogenization (WSL)
# AUTHOR: [Your Name]
# DATE: 2026
# 
# DESCRIPTION: 
# This script processes canopy throughfall eDNA (fungal ITS) data to 
# assign ecological functional guilds (FUNGuild/FungalTraits) and test 
# the hypothesis that non-native urban trees exhibit functional 
# homogenization compared to native congeners.
#
# DEPENDENCIES: microeco, ggplot2, tidyr, dplyr
# INPUTS: 
#   - asv_table.csv (ASVs as rows, samples as columns)
#   - taxonomy_table.csv (ASVs as rows, taxonomic ranks as columns)
#
# OUTPUTS:
#   - guild_comparison_results.csv (Statistical summary table)
#   - Figure_2_Functional_Homogenization.pdf (Paired abundance plot)
# =====================================================================

# 0. Install and load required libraries
# install.packages(c("microeco", "ggplot2", "tidyr", "dplyr"))
library(microeco)
library(ggplot2)
library(tidyr)
library(dplyr)

# Set seed for reproducibility
set.seed(42)

# =====================================================================
# 1. LOAD DATA
# =====================================================================
# Note: Ensure your working directory is set to the folder containing the CSVs.
# For Colab, the path is usually "/content/"

message("Loading ASV and Taxonomy tables...")
asv_data <- read.csv("asv_table.csv", row.names = 1, check.names = FALSE)
tax_data <- read.csv("taxonomy_table.csv", row.names = 1, check.names = FALSE)

# Create the microtable object
dataset <- microtable$new(otu_table = asv_data, tax_table = tax_data)
print(dataset)

# =====================================================================
# 2. ASSIGN FUNCTIONAL GUILDS (FUNGuild)
# =====================================================================
message("Assigning functional guilds using FUNGuild...")
t1 <- trans_func$new(dataset = dataset)

# Assign guilds (using strict confidence levels)
t1$cal_func(fungi_database = "FUNGuild", 
            FUNGuild_confidence = c("Highly Probable", "Probable"))

# Extract the binary guild matrix (1 = present, 0 = absent)
guild_matrix <- t1$res_func

# =====================================================================
# 3. CALCULATE RELATIVE ABUNDANCES
# =====================================================================
message("Calculating guild proportions per sample...")

# Match ASVs between abundance and guild matrices
common_asvs <- intersect(rownames(asv_data), rownames(guild_matrix))
asv_abundance <- as.matrix(asv_data[common_asvs, , drop = FALSE])
guild_matrix <- guild_matrix[common_asvs, , drop = FALSE]

# Multiply abundance by guild presence to get guild abundance per sample
guild_abundance <- t(asv_abundance) %*% as.matrix(guild_matrix)

# Convert to relative abundance (proportion of total assigned fungi)
guild_proportions <- sweep(guild_abundance, 1, rowSums(guild_abundance), "/")
df_prop <- as.data.frame(guild_proportions)

# =====================================================================
# 4. METADATA PARSING & RESHAPING
# =====================================================================
message("Parsing sample metadata...")

# Extract metadata cleanly using string splitting
# Assumes sample names are formatted as "Sample_Native_1" or "Sample_NonNative_1"
df_prop$Sample <- rownames(df_prop)
parts <- strsplit(df_prop$Sample, "_")
df_prop$Origin  <- sapply(parts, function(x) x[2])  
df_prop$Pair_ID <- sapply(parts, function(x) x[3])  

# =====================================================================
# 5. STATISTICAL TESTING (WILCOXON PAIRED LOOP)
# =====================================================================
message("Running paired Wilcoxon tests for all guilds...")

guilds <- setdiff(colnames(df_prop), c("Sample", "Origin", "Pair_ID"))
results_list <- list()

for (g in guilds) {
  # Reshape to wide format for the specific guild
  df_wide <- reshape(df_prop[, c("Pair_ID", "Origin", g)], 
                     idvar = "Pair_ID", 
                     timevar = "Origin", 
                     direction = "wide")
  
  col_native <- paste0(g, ".Native")
  col_non_native <- paste0(g, ".NonNative")
  
  if (col_native %in% colnames(df_wide) & col_non_native %in% colnames(df_wide)) {
    
    # Paired Wilcoxon test
    test <- wilcox.test(df_wide[[col_native]], df_wide[[col_non_native]], paired = TRUE)
    
    # Calculate means
    mean_nat <- mean(df_wide[[col_native]], na.rm = TRUE)
    mean_non <- mean(df_wide[[col_non_native]], na.rm = TRUE)
    
    # Store
    results_list[[g]] <- data.frame(
      Guild = g,
      Mean_Native = round(mean_nat, 4),
      Mean_NonNative = round(mean_non, 4),
      P_value = round(test$p.value, 4)
    )
  }
}

# Combine and process results
final_results <- do.call(rbind, results_list)
final_results$Higher_Abundance <- ifelse(final_results$Mean_Native > final_results$Mean_NonNative, 
                                         "Native", 
                                         ifelse(final_results$Mean_NonNative > final_results$Mean_Native, 
                                                "NonNative", "Equal"))
final_results <- final_results[order(final_results$P_value), ]

# Save results table
write.csv(final_results, "guild_comparison_results.csv", row.names = FALSE)

# =====================================================================
# 6. VISUALIZATION (For Ectomycorrhizal Shift)
# =====================================================================
message("Generating paired abundance plot...")

# Create the plot (focusing on the Ectomycorrhizal guild as the primary example)
p <- ggplot(df_prop, aes(x = Origin, y = Ectomycorrhizal, group = Pair_ID)) +
  geom_line(color = "gray", alpha = 0.5) + 
  geom_point(aes(color = Origin), size = 4) +
  theme_minimal() +
  labs(
    title = "Functional Homogenization: Ectomycorrhizal Fungi",
    subtitle = "Paired comparison of native and non-native urban trees",
    y = "Relative Abundance (Proportion)",
    x = "Tree Origin"
  ) +
  scale_color_manual(values = c("Native" = "forestgreen", "NonNative" = "firebrick")) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

# Save the plot as a high-resolution PDF for publication
ggsave("Figure_2_Functional_Homogenization.pdf", plot = p, width = 6, height = 6, dpi = 300)

message("Analysis complete. Results saved to working directory.")
