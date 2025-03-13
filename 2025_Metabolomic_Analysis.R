##Set your working Directory
setwd("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/")
# Load necessary libraries
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(rio)
library(tidyverse)
library(openxlsx)
library(tibble)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(broom.mixed)
library(performance)
library(lme4)
library(lmerTest)  # Enables p-values for mixed models
library(readr)
# Files Needed;
#1. Metabolomics Excel sheet from Asael
#2. Sample map
#3. Latest Frozen_Data by Kristina
#4. List of Unique_ID for which we have DNA
## Read the DNA unique_IDs
Gene_IDS <- import("Gene_Unique_IDs.csv")
head(Gene_IDS)
# read in the HILIC results and remove the columns that are not samples or pooled QC samples
allData <- as.data.frame(import( "HILIC_raw_from_Asael.xlsx", sheet = "Organized"))
head(allData)
rownames( allData) <- allData$compoundId
dataSub <- allData[ , -c( 1:9)] ##Removing the Blanks
# order the data by injection order (pooled and repeat samples are listed at the end in the original file)
dataSub.ordered <- dataSub[ , c( 618, 1:25, 619, 26:50, 620, 51:75, 621, 76:100, 622, 101:125, 623, 126:150, 624, 151:175, 625, 176:200, 626, 201:239, 629, 240:264, 630, 265:289, 631, 290:314, 632, 315:339, 633, 340:358, 634, 610:612, 359:380, 635, 381:397, 636, 398:422, 637, 423:447, 638, 448:476, 639, 477:501, 640, 502:526, 641, 527:551, 642, 552:560, 613:617, 561:571, 643, 572:596, 644, 597:609)]
## Group the Batches and create an Ordered Batch Map
batchMap <- rbind( data.frame( batch = 1, sample = colnames( dataSub.ordered)[ 1:248]), data.frame( batch = 2, sample = colnames( dataSub.ordered)[ 249:372]), data.frame( batch = 3, sample = colnames( dataSub.ordered)[ 373:416]), data.frame( batch = 4, sample = colnames( dataSub.ordered)[ 417:498]), data.frame( batch = 5.1, sample = colnames( dataSub.ordered)[ 499:586]), data.frame( batch = 5.2, sample = colnames( dataSub.ordered)[ 587:642]))
batchMap$batch <- factor( batchMap$batch, levels = c( 1, 2, 3, 4, 5.1, 5.2))
## Create an additional row with the Injection order into the Machine
order <- batchMap$sample
batchMap = batchMap %>% 
  mutate(Run.ID=row_number())
head(batchMap)
table(batchMap$Run.ID)
## ADD METADATA YOU WANT TO USE - Use Updated Frozen data Each time
All_data <- import("THGP_database_TurkanaOnly_corrected_merged_2025-01-28.txt")
dim(All_data)
table(All_data$MW_scores_market_diet_index)
table(All_data$Coffee.yes.no)
table(All_data$Main.subsistence.activity)
## Choose whatever Metadata you want to work with
Meta_ONLY <- All_data %>% dplyr::select(c("Unique.ID","MW_scores_h_sol_trim", "MW_scores_lifestyle", "Main.subsistence.activity", "Age", "Sex", "MW_scores_market_diet_index"))
Meta_ONLY <- Meta_ONLY %>% filter(Sex %in% c("Male", "Female"))
## Sample Map
sampleMap <- as.data.frame(import("/Users/bm0211/Documents/AyrolesLab/metabolomics_sampleMap_final.csv"))
## Merge Map and Unique.ID data 
batchMap_with_Meta <- merge(sampleMap, batchMap, by=2, all = T)
batchMap_with_Meta <- unique(merge(batchMap_with_Meta, Meta_ONLY, by="Unique.ID", all.x = T))
batchMap_with_Meta <- dplyr::select(batchMap_with_Meta, -(3:8))
dim(batchMap_with_Meta)
colnames(batchMap_with_Meta)
table(batchMap_with_Meta$MW_scores_h_sol_trim)
## Ordered metabolites file
dataSub.ordered = rownames_to_column(dataSub.ordered, var = "Metabolites")
dim(dataSub.ordered)
# Pivot and merge the data WITHOUT normalization
NARROW_dataSub.ordered <- dataSub.ordered %>% 
  pivot_longer(!Metabolites, names_to = "indiv.ID", values_to = "Concentration")
# Join HILIC_Data with Meta_data
NARROW_dataSub.ordered <- left_join(NARROW_dataSub.ordered, batchMap_with_Meta, by = c('indiv.ID' = 'metabRunNum'))
# Separate metadata columns for later merge
metadata_columns <- NARROW_dataSub.ordered %>%
  dplyr::select(-Metabolites, -Concentration) %>%
  distinct(indiv.ID, .keep_all = TRUE)
# Pivot the data back to wide format using raw concentration values (NO normalization)
WIDE_dataSub <- NARROW_dataSub.ordered %>%
  pivot_wider(
    id_cols = indiv.ID,           # Unique identifier for individuals
    names_from = Metabolites,     # Use metabolites as column names
    values_from = Concentration   # Keep raw concentration values
  )
# Merge the metadata back into the wide dataset
WIDE_dataSub_with_metadata <- WIDE_dataSub %>%
  left_join(metadata_columns, by = "indiv.ID")
# Check dimensions of the new wide dataset
dim(WIDE_dataSub_with_metadata)
colnames(WIDE_dataSub_with_metadata)
# Remove duplicates in Unique.ID, keeping the first occurrence
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  distinct(Unique.ID, .keep_all = TRUE)
table(is.na(WIDE_dataSub_with_metadata$MW_scores_lifestyle))
# Grouping lifestyles
##The Marina Score Below
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  mutate(MW_scores_lifestyle = ifelse(MW_scores_lifestyle == "Urban", "Urban", "Non_Urban")) %>%
  filter(!is.na(MW_scores_lifestyle))
table(is.na(WIDE_dataSub_with_metadata$MW_scores_lifestyle))
table(WIDE_dataSub_with_metadata$MW_scores_lifestyle)
colnames(WIDE_dataSub_with_metadata)
# Save the final wide dataset to a CSV
For_writing <- dplyr::select(WIDE_dataSub_with_metadata, c(630:638,1:629))
write.csv(For_writing, "HILIC_RAW_DATA_with_covariates.csv", row.names = FALSE)


# Read the Updated_HMDB_NAMES file # Download from Google Sheets
updated_hmdb_names <- read.csv("Updated_No_match_Metabo.csv", stringsAsFactors = FALSE)
head(updated_hmdb_names)
table(updated_hmdb_names$HMDB)
# Print only the first 20 duplicated HMDB entries
updated_hmdb_names %>% 
  filter(HMDB %in% HMDB[duplicated(HMDB)]) %>%
  head(13) %>%
  print()
# Create a named vector for renaming (excluding NA values)
rename_vector <- setNames(updated_hmdb_names$HMDB, updated_hmdb_names$Query)
rename_vector <- rename_vector[!is.na(rename_vector)]  # Remove NA values
# Rename matching columns in For_writing
colnames(For_writing) <- ifelse(
  colnames(For_writing) %in% names(rename_vector),
  rename_vector[colnames(For_writing)],
  colnames(For_writing)  # Keep original name if no match
)
# Check renamed columns
colnames(For_writing)
table(For_writing$MW_scores_lifestyle)
# Save the final wide dataset with HMDB names
write.csv(For_writing, "HILIC_RAW_DATA_with_covs_HMDB_column_names.csv", row.names = FALSE)

# Pivot to long format for median centering and log transformation
NARROW_dataSub.ordered <- WIDE_dataSub_with_metadata %>%
  pivot_longer(
    cols = -c(indiv.ID, Unique.ID, MW_scores_h_sol_trim, MW_scores_lifestyle, Main.subsistence.activity, Age, Sex, Run.ID, batch),
    names_to = "Metabolites",
    values_to = "Concentration"
  )
head(NARROW_dataSub.ordered)
# Median center and log transform
NARROW_dataSub.ordered_Norm <- NARROW_dataSub.ordered %>%
  group_by(Metabolites) %>%
  mutate(med = median(Concentration, na.rm = TRUE)) %>%
  mutate(normIC = Concentration / med) %>%
  ungroup() %>%
  mutate(min.val = min(abs(normIC[normIC != 0]), na.rm = TRUE) / 10) %>%
  group_by(indiv.ID) %>%
  mutate(normLogIC = log10((normIC + sqrt(normIC^2 + min.val^2)) / 2)) %>%
  ungroup()
head(NARROW_dataSub.ordered_Norm)
# Separate metadata columns for later merging
metadata_columns <- NARROW_dataSub.ordered_Norm %>%
  dplyr::select(-Metabolites, -normLogIC) %>%
  distinct(indiv.ID, .keep_all = TRUE)
# Pivot back to wide format using normLogIC
WIDE_dataSub_normLogIC <- NARROW_dataSub.ordered_Norm %>%
  pivot_wider(
    id_cols = indiv.ID,
    names_from = Metabolites,
    values_from = normLogIC
  )
# Merge the metadata back into the wide dataset
WIDE_dataSub_with_metadata <- WIDE_dataSub_normLogIC %>%
  left_join(metadata_columns, by = "indiv.ID")
colnames(WIDE_dataSub_with_metadata)
Second_write <- dplyr::select(WIDE_dataSub_with_metadata, c(630:638, 1:629))
# Save the transformed dataset
write.csv(Second_write, "HILIC_DATA_NormLog_plus_Metadata.csv", row.names = FALSE)
# Confirm lifestyle groups
table(Second_write$MW_scores_lifestyle)
# Create datasets for each lifestyle
Non_Urban_data <- subset(Second_write, MW_scores_lifestyle == "Non_Urban")
Urban_data <- subset(Second_write, MW_scores_lifestyle == "Urban")
# Save datasets as CSV files
write.csv(Non_Urban_data, "HILIC_DATA_NormLog_Non_Urban_data.csv", row.names = FALSE)
write.csv(Urban_data, "HILIC_DATA_NormLog_Urban_data.csv", row.names = FALSE)


#The above datasets were made for use online in the Metabo-analyst pathway enrichment 
##
## Below, I process the data further and do my analyses. 
# Select and clean up unnecessary columns - FOR THE DOWNSTREAM ANALYSIS FROM HERE 
WIDE_dataSub_with_metadata <- dplyr::select(WIDE_dataSub_with_metadata, -c(639:642))
colnames(WIDE_dataSub_with_metadata)
##Collecting PC1 after doing a PCA on the residuals of the model
## Model; Metabolite ~ MW_scores_h_sol_trim + Age + Sex + batch
# Initialize a data frame to store residuals
residuals_data <- data.frame()
# Define the range of metabolite columns (2:629)
metabolite_columns <- 2:629
# Iterate over each metabolite column to perform regression and extract residuals
for (col_index in metabolite_columns) {
  metabolite_name <- colnames(WIDE_dataSub_with_metadata)[col_index]
  # Extract data for the current metabolite
  current_data <- WIDE_dataSub_with_metadata %>%
    dplyr::select(indiv.ID, all_of(metabolite_name), MW_scores_h_sol_trim, Age, Sex, batch) %>%
    rename(Metabolite = all_of(metabolite_name)) %>%
    filter(!is.na(Metabolite) & !is.na(MW_scores_h_sol_trim) & !is.na(Age) & !is.na(Sex) & !is.na(batch))
  # Check if sufficient data is available for regression
  if (nrow(current_data) < 2 || any(sapply(current_data[, c("Metabolite", "MW_scores_h_sol_trim", "Age", "Sex", "batch")], function(col) length(unique(col)) == 1))) {
    message(paste("Skipping metabolite:", metabolite_name, "- insufficient or constant data"))
    next
  }
  # Fit the regression model with error handling
  model <- tryCatch(
    lm(Metabolite ~ MW_scores_h_sol_trim + Age + Sex + batch, data = current_data),
    error = function(e) NULL
  )
  # Check if the model was successfully fit
  if (is.null(model)) {
    message(paste("Skipping metabolite:", metabolite_name, "- model fitting failed"))
    next
  }
  # Extract residuals and ensure the length matches the data
  if (length(residuals(model)) == nrow(current_data)) {
    current_residuals <- data.frame(
      indiv.ID = current_data$indiv.ID,
      Metabolite = metabolite_name,
      Residual = residuals(model)
    )
    # Combine residuals into a single data frame
    residuals_data <- bind_rows(residuals_data, current_residuals)
  } else {
    message(paste("Skipping metabolite:", metabolite_name, "- residual mismatch"))
  }
}
# Check if residuals data is available
if (nrow(residuals_data) > 0) {
  # Pivot residuals to wide format for PCA
  residuals_wide <- residuals_data %>%
    pivot_wider(names_from = Metabolite, values_from = Residual) %>%
    left_join(WIDE_dataSub_with_metadata %>% dplyr::select(indiv.ID, batch), by = "indiv.ID") %>%
    drop_na() # Remove rows with missing values
  # Ensure only numeric columns are passed to PCA
  metabolite_residuals <- residuals_wide %>%
    dplyr::select(-c(indiv.ID, batch)) %>%
    select_if(is.numeric)
  # Perform PCA on residuals
  pca_result <- prcomp(metabolite_residuals, center = TRUE, scale. = TRUE)
  # Extract PCA scores
  pca_scores <- as.data.frame(pca_result$x)
  pca_scores <- pca_scores %>%
    dplyr::select(PC1, PC2, PC3) %>%
    mutate(indiv.ID = residuals_wide$indiv.ID)
  # Add PCA scores to the original dataset
  WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
    left_join(pca_scores, by = "indiv.ID")
  # Save the updated dataset with PCA scores --> This is to use with Metaboanalyst for regression and volcano plot
  write.csv(WIDE_dataSub_with_metadata, "HILIC_data_with_metadata_and_PCs.csv", row.names = FALSE)
  # Plot the variance explained by PC1 and PC2
  pca_variance <- summary(pca_result)$importance[2, 1:3] * 100
  pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = residuals_wide$batch)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = "PCA of Metabolite Residuals",
      x = paste0("PC1 (", round(pca_variance[1], 2), "% variance)"),
      y = paste0("PC2 (", round(pca_variance[2], 2), "% variance)")
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
  
  # Save the PCA plot as a PDF
  ggsave("PCA_metabolite_residuals.pdf", plot = pca_plot, width = 8, height = 6)
} else {
  message("No valid residuals data for PCA.")
}
####
##########
####
#Regression - mixed effect model with Batch as a random effect
WIDE_dataSub_with_metadata$MW_scores_h_sol_trim <- as.numeric(WIDE_dataSub_with_metadata$MW_scores_h_sol_trim)
WIDE_dataSub_with_metadata$batch <- as.factor(WIDE_dataSub_with_metadata$batch)
table(WIDE_dataSub_with_metadata$batch)
# Rename batch values
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  mutate(batch = recode(batch, 
                        `1` = "one", 
                        `2` = "two", 
                        `3` = "three",
                        `4` = "four",
                        `5.1` = "five", 
                        `5.2` = "six"))
#######
######
# Initialize an empty dataframe to store coefficients
all_coefficients <- data.frame()
# Initialize an empty dataframe to store R-squared values
all_r_squared <- data.frame()
# Define the range of metabolite columns
metabolite_columns <- 2:629  # Ensure correct range
# Loop through each metabolite and run the regression
for (metabolite in colnames(WIDE_dataSub_with_metadata)[metabolite_columns]) {
  # Filter data for the current metabolite
  current_data <- WIDE_dataSub_with_metadata %>%
    dplyr::select(indiv.ID, MW_scores_lifestyle, Age, Sex, batch, PC1, MW_scores_market_diet_index, all_of(metabolite)) %>%
    rename(Metabolite = all_of(metabolite)) %>%
    filter(
      !is.na(Metabolite),
      !is.na(Age),
      !is.na(Sex),
      !is.na(PC1)
    )
  # Remove outliers using the IQR method
  if (nrow(current_data) > 0) {
    Q1 <- quantile(current_data$Metabolite, 0.25, na.rm = TRUE)
    Q3 <- quantile(current_data$Metabolite, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    current_data <- current_data %>%
      filter(Metabolite >= lower_bound & Metabolite <= upper_bound)
  }
  # Ensure enough data remains to fit the model
  if (nrow(current_data) > 6) {
    # Convert batch to a factor
    current_data <- current_data %>%
      mutate(batch = as.factor(batch))
    # Run the mixed-effects model with batch as a random effect
    model <- lmer(Metabolite ~ MW_scores_lifestyle + Sex + Age + PC1 + 
                    MW_scores_lifestyle:Sex + (1 | batch), data = current_data)
    # Extract fixed effect coefficients, including p-values
    temp_results <- broom.mixed::tidy(model, effects = "fixed") %>%
      mutate(
        Metabolite = metabolite,
        N = nrow(current_data)  # Add sample size
      )
    # Extract R-squared values (marginal & conditional R²)
    model_r_squared <- performance::r2(model) %>%
      as.data.frame() %>%
      rename(r.squared = R2_marginal, r.squared_conditional = R2_conditional) %>%
      mutate(Metabolite = metabolite)
    # Append results to master dataframes
    all_coefficients <- bind_rows(all_coefficients, temp_results)
    all_r_squared <- bind_rows(all_r_squared, model_r_squared)
  } else {
    message(paste("Skipping metabolite:", metabolite, "- insufficient data after removing outliers"))
  }
}
# Compute confidence intervals
all_coefficients <- all_coefficients %>%
  mutate(
    lower_CI = estimate - 1.96 * std.error,
    upper_CI = estimate + 1.96 * std.error
  )

# Apply Bonferroni correction to p-values
all_coefficients <- all_coefficients %>%
  mutate(adjusted_p_values = p.adjust(p.value, method = "bonferroni"))
# Save the results to CSV
write_csv(all_coefficients, "all_coefficients_results_lme.csv")
#########
###########
## Plots
## coefficient plot
# Define custom color
color_mapping_b <- c(
  "Age" = "darkorange",
  "MW_scores_lifestyleUrban" = "green",
  "SexMale" = "blue",
  "MW_scores_lifestyleUrban:SexMale" = "cyan"
)
# Filter: Remove intercept & irrelevant terms - Filter for significant coefficient
table(filtered_coefficients$term)
filtered_coefficients <- all_coefficients %>%
  filter(
    !term %in% c("(Intercept)", "PC1"),
    lower_CI > 0 | upper_CI < 0,  # Ensures significance (CI does not include 0)
    abs(estimate) >= 0.3  # Ensures moderate to strong effect
  ) %>%
  arrange(desc(abs(estimate)))  # Sort by effect size for clearer visualization
# Create the coefficient plot
Co_plot <- ggplot(filtered_coefficients, aes(x = estimate, y = reorder(Metabolite, estimate), color = term)) +
  geom_point(size = 3) +  # Points for estimated effects
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2) +  # Confidence intervals
  scale_color_manual(values = color_mapping_b) +  # Custom colors
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Reference line at 0
  labs(
    title = "Significant Moderate to Strong Effects on Metabolites (Mixed Model)",
    x = "Effect Size (Beta Estimate) Filtered for |β| ≥ 0.3 & Significant 95% CI",
    y = "Metabolite",
    color = "Predictor"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
# Print the Co_plot
print(Co_plot)
# Save the combined plot
ggsave("LME_coefficient_plot.png", plot = Co_plot, width = 10, height = 8, dpi = 300)

## Term Count Plot + Volcano Plot
# Filter to find significant results based on P-values
significant_results <- all_coefficients %>%
  filter(adjusted_p_values < 0.05) %>%
  arrange(adjusted_p_values)
# Remove unwanted terms like intercept, PC1
significant_results_no_intercept <- significant_results %>%
  filter(!term %in% c("(Intercept)", "PC1"))
# Summarize data by predictors (MW_scores_h_sol_trim, Age, Sex, interaction terms)
term_counts <- significant_results_no_intercept %>%
  group_by(term) %>%
  summarise(Count = n_distinct(Metabolite), .groups = 'drop')
# Mapping of old term names to new labels
term_labels <- c(
  "Age" = "Age",       
  "MW_scores_lifestyleUrban" = "Lifestyle",
  "MW_scores_lifestyleUrban:SexMale" = "Lifestyle_X_Sex_Interaction",
  "MW_scores_market_diet_index" = "Market_Diet_Index",     
  "SexMale" = "Sex",
  "Non-Significant" = "Non-Significant"
)
# Apply renaming to term_counts and pie_data
term_counts <- term_counts %>%
  mutate(term = recode(term, !!!term_labels))
# Total number of metabolites tested
total_metabolites <- n_distinct(all_coefficients$Metabolite)
# Compute number of non-significant metabolites
non_significant_count <- total_metabolites - sum(term_counts$Count)
# Create pie chart data including non-significant metabolites
pie_data <- term_counts %>%
  mutate(Percentage = (Count / total_metabolites) * 100) %>%
  add_row(term = "Non-Significant", Count = non_significant_count, Percentage = (non_significant_count / total_metabolites) * 100)
# Define custom colors for each predictor
custom_colors <- c(
  "Age" = "darkorange",       
  "Lifestyle" = "Green",
  "Lifestyle_X_Sex_Interaction" = "Red",
  "Market_Diet_Index" = "cyan",     
  "Sex" = "Blue",                       
  "Non-Significant" = "azure2"
)
# Bar plot: Number of significant metabolites per predictor
bar_plot <- ggplot(term_counts, aes(x = reorder(term, Count), y = Count, fill = term)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  labs(
    title = "Count of Significant Metabolites by Predictor",
    x = "Predictor",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
# Calculate the midpoints of each segment for accurate label positioning
pie_data <- pie_data %>%
  arrange(desc(term)) %>%
  mutate(ypos = cumsum(Percentage) - (Percentage / 2))
pie_chart <- ggplot(pie_data, aes(x = "", y = Percentage, fill = term)) +
  geom_bar(stat = "identity", width = 1) +  # White border to separate slices
  coord_polar(theta = "y", clip = "off") +  # Prevents clipping of labels
  geom_text_repel(
    aes(y = ypos, label = paste0(round(Percentage, 1), "%")),  
    nudge_x = 0.8,   # Move labels outward
    size = 5, 
    direction = "y",  # Keeps text aligned to segments
    box.padding = 0.2,
    point.padding = 0.3,
    segment.size = 0.7,
    force = 1,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = custom_colors) +  # Apply same colors as bar plot
  labs(fill = "Predictor") +
  theme_void()
# Combine the bar plot and the pie chart
combined_plot <- ggdraw() +
  draw_plot(bar_plot, 0, 0, 1, 1) +  # Main bar plot
  draw_plot(pie_chart, 0.1, 0.45, 0.5, 0.5)  # Inset pie chart (top left)
# Print the combined plot
print(combined_plot)
# Save the combined plot
ggsave("Significant_Metabolites_BarPlot_with_PieChart_Fixed.png", plot = combined_plot, width = 10, height = 8, dpi = 300)

# Identify top 50 metabolites by effect size for Lifestyle component
top_50_metabolites_h_sol <- significant_results %>%
  filter(term == "MW_scores_lifestyleUrban") %>%
  arrange(desc(abs(estimate))) %>%
  head(10)
# Create a bar plot for the top 50 metabolites by signed effect size with \(N\) on top
top_50_plot_signed <- ggplot(top_50_metabolites_h_sol, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) + # Add \(N\) as text on top of bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites: 
    Ordered by lifestyle's effect size on the Metabolite",
    x = "Metabolite",
    y = "Beta value"
  ) +
  coord_flip() +
  theme_minimal()
# Print the top 50 plot
print(top_50_plot_signed)
# Save the plot
ggsave("Top_Metabolites_Discrete_lifestyle.png", plot = top_50_plot_signed, width = 8, height = 6, dpi = 300)
#
# Volcano Plot with the results from the regression combined with a fold change analysis. 
#
# Calculate log2 fold change
all_coefficients <- all_coefficients %>%
  mutate(log2_fold_change = estimate * log2(10))  # Convert log10-based coefficient to log2 fold change
# Define custom fold change threshold for significant metabolites
fc_threshold <- 1.2  
# Adjust volcano plot with log2 FC
volcano_data <- all_coefficients %>%
  filter(term == "MW_scores_lifestyleUrban") %>%
  mutate(
    log_p_value = -log10(adjusted_p_values),
    significant = adjusted_p_values < 0.00005 & abs(log2_fold_change) > fc_threshold,
    regulation = case_when(
      log2_fold_change > fc_threshold & adjusted_p_values < 0.00005 ~ "Upregulated in Urban",
      log2_fold_change < -fc_threshold & adjusted_p_values < 0.00005 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )
# Volcano Plot with different colors for up/down regulation and labels for significant points
volcano_plot <- ggplot(volcano_data, aes(x = log2_fold_change, y = log_p_value)) +
  geom_point(aes(color = regulation), alpha = 0.8, size = 3) + 
  scale_color_manual(values = c("Upregulated in Urban" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.00005), linetype = "dashed", color = "black") + 
  ggrepel::geom_text_repel(
    data = volcano_data %>% filter(significant), # Only label significant metabolites
    aes(label = Metabolite),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 15  # Adjust the number of labels shown to reduce clutter
  ) +
  labs(
    title = "Log2 Fold Change of Metabolites Urban vs Non_Urban",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme_minimal()
# Print and save the plot
print(volcano_plot)
ggsave("Volcano_Plot_Regression_plus_FC.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
# Identify top 50 metabolites by effect size for Age
top_50_metabolites_age <- significant_results_no_intercept %>%
  filter(term == "Age") %>%
  arrange(desc(abs(estimate))) %>%
  head(10)
# Create a bar plot for the top 50 metabolites by signed effect size for Age with \(N\) displayed
top_50_plot_age <- ggplot(top_50_metabolites_age, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Display \(N\) above the bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites Affected by Age",
    x = "Metabolite",
    y = "Beta value"
  ) +
  coord_flip() +
  theme_minimal()
# Print the plot for Age
print(top_50_plot_age)
# Save the plot for Age
ggsave("Top_Metabolites_Affected_by_Age.png", plot = top_50_plot_age, width = 8, height = 6, dpi = 300)

# Identify top 50 metabolites by effect size for Sex
# Calculate log2 fold change - VOLCANO PLOT FOR SEX DIFFERENCES
all_coefficients <- all_coefficients %>%
  mutate(log2_fold_change = log2(exp(estimate)))  # Convert regression estimate to log2 fold change
# Define fold change threshold
fc_threshold <- 0.4  
# Prepare volcano plot data for Sex comparison
volcano_data_sex <- all_coefficients %>%
  filter(term == "SexMale") %>%
  mutate(
    log_p_value = -log10(adjusted_p_values),
    significant = adjusted_p_values < 0.05 & abs(log2_fold_change) > fc_threshold,
    regulation = case_when(
      log2_fold_change > fc_threshold & adjusted_p_values < 0.05 ~ "Higher in Males",
      log2_fold_change < -fc_threshold & adjusted_p_values < 0.05 ~ "Higher in Females",
      TRUE ~ "Not Significant"
    )
  )

# Create volcano plot for Male vs Female metabolite differences
volcano_plot_sex <- ggplot(volcano_data_sex, aes(x = log2_fold_change, y = log_p_value)) +
  geom_point(aes(color = regulation), alpha = 0.8, size = 3) + 
  scale_color_manual(
    values = c("Higher in Males" = "blue", "Higher in Females" = "red", "Not Significant" = "grey")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  
  ggrepel::geom_text_repel(
    data = volcano_data_sex %>% filter(significant),
    aes(label = Metabolite),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 15  
  ) +
  labs(
    title = "Volcano Plot: Metabolite Differences Between Males and Females",
    x = "Log2 Fold Change (Male/Female)",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme_minimal()
# Print and save the plot
print(volcano_plot_sex)
ggsave("Volcano_Plot_Sex_Corrected.png", plot = volcano_plot_sex, width = 8, height = 6, dpi = 300)

# Filter top 5 metabolites with strongest interaction effects
top_metabs <- significant_results_no_intercept %>%
  filter(term == "MW_scores_lifestyleUrban:SexMale") %>%
  arrange(desc(abs(estimate))) %>%
  head(5)  # Select top 5 metabolites
# Prepare data for plotting
interaction_data <- NARROW_dataSub.ordered_Norm %>%
  filter(Metabolites %in% top_metabs$Metabolite)
# Create scatter plot with regression lines for each metabolite
interaction_plot <- ggplot(interaction_data, aes(x = MW_scores_h_sol_trim, y = normLogIC, color = Sex)) +
  geom_point(alpha = 0.6) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE) +  # Regression lines with confidence intervals
  facet_wrap(~ Metabolites, scales = "free_y") +  # Separate plots for each metabolite
  labs(
    title = "Interaction Effect of MW_scores_h_sol_trim and Sex on Metabolites",
    x = "MW_scores_h_sol_trim",
    y = "Normalized Log Concentration",
    color = "Sex"
  ) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +  # Consistent colors
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
# Print the plot
print(interaction_plot)
# Save the plot as a PNG
ggsave("Interaction_Effect_Sex_H_sol.png", plot = interaction_plot, width = 6, height = 8, dpi = 300)


#########
#############
#################
############
########

# Compiling names to match more metabolites to HMDB using the Metaboanalyst Platform
# Here I match HMDB acsession numbers given to me by Asael, and I match them to my data
# I also have a list of names of those that Metaboanalyst automatically matched and those that it didn't
# The goal was to increase the number of those that have a HMDB Ascession number 

# Import datasets with UTF-8 encoding
Metabolites_Asael <- read_csv("Metabolites_Asael.csv", locale = locale(encoding = "UTF-8"))
No_match_Metabo <- read_csv("Metaboanalyst_match_no_match.csv", locale = locale(encoding = "UTF-8"))

# Function to find exact or close match without modifying characters
find_best_match <- function(compound_name, match_list) {
  # Attempt direct matching
  exact_match <- match_list[match_list == compound_name]
  if (length(exact_match) > 0) {
    return(exact_match[1])
  }
  # If no exact match, use agrepl for fuzzy matching (95% threshold)
  fuzzy_matches <- match_list[agrepl(compound_name, match_list, max.distance = 0.05, ignore.case = TRUE)]
  if (length(fuzzy_matches) > 0) {
    return(fuzzy_matches[1])  # Take the first close match
  } else {
    return(NA)  # No good match found
  }
}
# Identify rows in No_match_Metabo with missing HMDB
missing_hmdb_rows <- No_match_Metabo %>%
  filter(is.na(HMDB) | HMDB == "") %>%
  dplyr::select(Query)
# Track skipped rows
skipped_queries <- c()
# Perform matching for missing HMDB entries
matched_data <- missing_hmdb_rows %>%
  rowwise() %>%
  mutate(
    Matched_Compound = tryCatch(
      find_best_match(Query, Metabolites_Asael$Compound),
      error = function(e) {
        skipped_queries <<- c(skipped_queries, Query)  # Track problematic row
        return(NA)  # Skip problematic row
      }
    )
  ) %>%
  left_join(Metabolites_Asael, by = c("Matched_Compound" = "Compound")) %>%
  dplyr::select(Query, Matched_Compound, HMDB) %>%
  filter(!is.na(HMDB))  # Keep only successful matches
# Update No_match_Metabo with imputed HMDB values
No_match_Metabo <- No_match_Metabo %>%
  left_join(matched_data, by = "Query") %>%
  mutate(HMDB = coalesce(HMDB.x, HMDB.y)) %>%  # Use existing HMDB if available, else imputed
  dplyr::select(-HMDB.x, -HMDB.y)  # Remove extra columns

# Save the updated dataset and google the rest Manually on the HMDB website
write_csv(No_match_Metabo, "Asael_plus_metaboAnalyst.csv")
# Print skipped queries
if (length(skipped_queries) > 0) {
  cat("Skipped rows due to errors in matching:\n")
  print(skipped_queries)
} else {
  cat("No rows were skipped during processing.\n")
}

cat("Successfully created 'Asael_plus_metaboAnalyst.csv' with imputed HMDB values.\n")

###Making violin plots of those of interest 
# Filter data for the selected metabolite
Selected_metab <- WIDE_dataSub_with_metadata %>%
  dplyr::select(indiv.ID, MW_scores_lifestyle, `Ginsenoside Rh3 isomer`) %>%  # Select metabolite of interest
  filter(!is.na(`Ginsenoside Rh3 isomer`))  # Remove missing values
# Create Violin Plot with Jitter
metabolite_violin_plot <- ggplot(Selected_metab, aes(x = MW_scores_lifestyle, y = `Ginsenoside Rh3 isomer`, fill = MW_scores_lifestyle)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with transparency
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +  # Jitter to show individual points
  labs(
    title = "Comparison of Ginsenoside Levels by Lifestyle",
    x = "Lifestyle Category",
    y = "Ginsenoside (units)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend
  scale_fill_brewer(palette = "Set2")  # Use a color palette for better visualization
# Print the plot
print(metabolite_violin_plot)
# Save the plot
ggsave("Metabolite_ViolinPlot.png", plot = metabolite_violin_plot, width = 8, height = 6, dpi = 300)

##
##
colnames(WIDE_dataSub_with_metadata)
table(WIDE_dataSub_with_metadata$Unique.ID)
table(Gene_IDS$IID)
# Subset the dataset into Male and Female groups
female_data <- WIDE_dataSub_with_metadata %>% filter(Sex == "Female")
male_data <- WIDE_dataSub_with_metadata %>% filter(Sex == "Male")

# Function to perform regression analysis for a given dataset
run_regression <- function(data, sex_label) {
  all_coefficients <- data.frame()
  all_r_squared <- data.frame()
  
  metabolite_columns <- 2:629  # Define metabolite range (check column indices)
  
  for (metabolite in colnames(data)[metabolite_columns]) {
    # Filter and clean data
    current_data <- data %>%
      dplyr::select(indiv.ID, MW_scores_lifestyle, Age, batch, PC1, MW_scores_market_diet_index, all_of(metabolite)) %>%
      rename(Metabolite = all_of(metabolite)) %>%
      filter(
        !is.na(Metabolite),
        !is.na(Age),
        !is.na(PC1)
      )
    
    # Remove outliers using IQR method
    if (nrow(current_data) > 0) {
      Q1 <- quantile(current_data$Metabolite, 0.25, na.rm = TRUE)
      Q3 <- quantile(current_data$Metabolite, 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      current_data <- current_data %>%
        filter(Metabolite >= lower_bound & Metabolite <= upper_bound)
    }
    
    # Ensure enough data remains for regression
    if (nrow(current_data) > 6) {
      current_data <- current_data %>%
        mutate(batch = as.factor(batch))
      
      # Run regression model (without Sex)
      model <- lm(Metabolite ~ MW_scores_lifestyle + Age + PC1 + batch + MW_scores_market_diet_index, data = current_data)
      
      # Extract coefficients
      temp_results <- broom::tidy(model) %>%
        mutate(
          Metabolite = metabolite,
          Sex = sex_label,
          N = nrow(current_data)
        )
      
      # Extract R-squared
      model_r_squared <- broom::glance(model) %>%
        dplyr::select(r.squared) %>%
        mutate(Metabolite = metabolite, Sex = sex_label)
      
      # Append results
      all_coefficients <- bind_rows(all_coefficients, temp_results)
      all_r_squared <- bind_rows(all_r_squared, model_r_squared)
    }
  }
  
  # Apply Bonferroni correction
  all_coefficients <- all_coefficients %>%
    mutate(adjusted_p_values = p.adjust(p.value, method = "bonferroni"))
  
  return(all_coefficients)
}

# Run regression separately for males and females
female_results <- run_regression(female_data, "Female")
male_results <- run_regression(male_data, "Male")

# Combine results
all_results <- bind_rows(female_results, male_results)

# Filter significant results
significant_results <- all_results %>%
  filter(adjusted_p_values < 0.05) %>%
  arrange(adjusted_p_values)

# Remove unwanted terms like intercept, PC1, and batch
significant_results_filtered <- significant_results %>%
  filter(!term %in% c("(Intercept)", "PC1", "batch", "batch1", "batch2", "batch3", "batch4", "batch5.1", "batch5.2"))

# Define colors for visualization
custom_colors <- c(
  "Age" = "darkorange",       
  "MW_scores_lifestyleUrban" = "Green",
  "MW_scores_market_diet_index" = "cyan",
  "Non-Significant" = "azure2"
)

# Ensure Fold Change Calculation Is Based on Lifestyle
create_volcano_plot <- function(data, sex_label, remove_legend = FALSE) {
  data <- data %>%
    filter(term == "MW_scores_lifestyleUrban") %>%  # Ensure only lifestyle predictor
    mutate(
      log_p_value = -log10(adjusted_p_values),
      log2_fold_change = estimate * log2(10), # Correct conversion for log10 data
      significant = adjusted_p_values < 0.00005 & abs(log2_fold_change) > 1,
      regulation = case_when(
        log2_fold_change > 1 & adjusted_p_values < 0.00005 ~ "Upregulated in Urban",
        log2_fold_change < -1 & adjusted_p_values < 0.00005 ~ "Downregulated in Urban",
        TRUE ~ "Not Significant"
      )
    )
  
  plot <- ggplot(data, aes(x = log2_fold_change, y = log_p_value)) +
    geom_point(aes(color = regulation), alpha = 0.8, size = 3) + 
    scale_color_manual(values = c("Upregulated in Urban" = "red", 
                                  "Downregulated in Urban" = "blue", 
                                  "Not Significant" = "grey")) +
    geom_hline(yintercept = -log10(0.00005), linetype = "dashed", color = "black") + 
    ggrepel::geom_text_repel(
      data = data %>% filter(significant), 
      aes(label = Metabolite),
      size = 3,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = 8  
    ) +
    labs(
      title = paste("Volcano Plot comparing Lifestyles:", sex_label, "subset data"),
      x = "Log2 Fold Change (Urban vs Non-Urban)",
      y = "-log10(Adjusted P-Value)"
    ) +
    theme_minimal()
  
  # Remove legend if specified
  if (remove_legend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}

# Create separate volcano plots for lifestyle effects in males and females
volcano_plot_female <- create_volcano_plot(female_results, "Female", remove_legend = TRUE)  # Remove legend from female plot
volcano_plot_male <- create_volcano_plot(male_results, "Male")  # Keep legend in male plot

# Combine and save volcano plots
combined_volcano_plot <- volcano_plot_female + volcano_plot_male
ggsave("Volcano_Plot_Lifestyle_Female_vs_Male.png", plot = combined_volcano_plot, width = 12, height = 6, dpi = 300)

# Print the updated combined plot
print(combined_volcano_plot)


# Create bar plots for the number of significant metabolites per predictor
create_bar_plot <- function(data, sex_label) {
  
  # Ensure that we correctly summarize the significant metabolite counts by predictor term
  term_counts <- data %>%
    filter(adjusted_p_values < 0.05) %>%  # Only consider significant metabolites
    filter(!term %in% c("(Intercept)", "PC1", "batch", "batch1", "batch2", "batch3", "batch4", "batch5.1", "batch5.2")) %>%  # Remove unwanted terms
    group_by(term) %>%
    summarise(Count = n_distinct(Metabolite), .groups = 'drop')
  
  # Generate the bar plot
  ggplot(term_counts, aes(x = reorder(term, Count), y = Count, fill = term)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = custom_colors) +  
    labs(
      title = paste("Significant Metabolites by Predictor:", sex_label),
      x = "Predictor",
      y = "Number of Significant Metabolites"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

# Generate and save bar plots separately for males and females
bar_plot_female <- create_bar_plot(female_results, "Female")
bar_plot_male <- create_bar_plot(male_results, "Male")

# Combine the two plots side by side
combined_bar_plot <- bar_plot_female + bar_plot_male

# Save the combined bar plot
ggsave("Bar_Plot_Female_vs_Male.png", plot = combined_bar_plot, width = 12, height = 6, dpi = 300)

# Print the combined bar plot
print(combined_bar_plot)

library(ggplot2)
library(dplyr)
library(patchwork)  # Ensure patchwork is installed for side-by-side plots

# Function to create a bar plot for the top 10 metabolites by effect size for a given predictor
create_top_metabolites_plot <- function(data, predictor, sex_label) {
  top_metabolites <- data %>%
    filter(term == predictor) %>%  # Filter for the specific predictor
    arrange(desc(abs(estimate))) %>%
    head(20)  # Get the top 10 metabolites
  
  ggplot(top_metabolites, aes(x = reorder(Metabolite, estimate), y = estimate, fill = estimate > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
    labs(
      title = paste("Top Metabs:", predictor, "-", sex_label),
      x = "Metabolite",
      y = "Beta Value (Effect Size)"
    ) +
    coord_flip() +
    theme_minimal()
}

# Define the predictor to analyze
predictor <- "MW_scores_lifestyleUrban"

# Generate plots for Male and Female separately
plot_male <- create_top_metabolites_plot(male_results, predictor, "Male")
plot_female <- create_top_metabolites_plot(female_results, predictor, "Female")

# Combine plots side by side
combined_plots <- plot_female + plot_male

# Save the combined plot
ggsave("Top_10_Metabolites_Lifestyle_Female_vs_Male.png", plot = combined_plots, width = 12, height = 6, dpi = 300)

# Print the combined plot
print(combined_plots)

###
# Get unique IDs from both datasets
metadata_ids_unique <- unique(WIDE_dataSub_with_metadata$Unique.ID)

# Find the overlapping IDs
overlapping_ids <- intersect(gene_ids_unique, metadata_ids_unique)

# Count the number of overlapping IDs
num_overlapping <- length(overlapping_ids)

# Print results
print(paste("Number of overlapping Unique IDs:", num_overlapping))

# Load necessary library
library(dplyr)

# Ensure gene_ids_unique is a unique vector
gene_ids_unique <- unique(gene_ids_unique)

# Extract matching rows from WIDE_dataSub_with_metadata based on Unique.ID
matched_metadata <- WIDE_dataSub_with_metadata %>%
  filter(Unique.ID %in% gene_ids_unique)

# Ensure all gene_ids_unique are included in the final dataset
final_data <- data.frame(Unique.ID = gene_ids_unique) %>%
  left_join(matched_metadata, by = "Unique.ID")
dim(final_data)
# Write the final dataset to a CSV file
write.csv(final_data, "final_data_matched.csv", row.names = FALSE)

# Print summary
num_overlapping <- nrow(matched_metadata)
print(paste("Number of overlapping Unique IDs:", num_overlapping))
print(paste("Total rows in final dataset:", nrow(final_data)))



