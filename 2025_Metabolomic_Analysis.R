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
# Files Needed;
#1. Metabolomics Excel sheet from Asael
#2. Sample map
#3. Latest Frozen_Data by Kristina
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
## Choose whatever Metadata you want to work with
Meta_ONLY <- All_data %>% dplyr::select(c("Unique.ID","MW_scores_h_sol_trim", "MW_scores_lifestyle", "Main.subsistence.activity", "Age", "Sex"))
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
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  mutate(MW_scores_lifestyle = ifelse(MW_scores_lifestyle == "Urban", "Urban", "Non_Urban")) %>%
  filter(!is.na(MW_scores_lifestyle))
table(is.na(WIDE_dataSub_with_metadata$MW_scores_lifestyle))
# Save the final wide dataset to a CSV
For_writing <- dplyr::select(WIDE_dataSub_with_metadata, c(630:637,1:629))
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
Second_write <- dplyr::select(WIDE_dataSub_with_metadata, c(630:637, 1:629))
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
WIDE_dataSub_with_metadata <- dplyr::select(WIDE_dataSub_with_metadata, -c(638:641))
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
#Regression
WIDE_dataSub_with_metadata$MW_scores_h_sol_trim <- as.numeric(WIDE_dataSub_with_metadata$MW_scores_h_sol_trim)
##  For the regression, I ran one model at a time, 
## hence you have to switch which terms are in the interacting term of the model each time - do this manually
## This will allow you be sure which analysis you do
## MW_scores_h_sol_trim + Age + Sex + PC1 + batch + MW_scores_h_sol_trim:Sex
## MW_scores_h_sol_trim + Age + Sex + PC1 + batch + MW_scores_h_sol_trim:Age
# Initialize an empty dataframe to store coefficients
all_coefficients <- data.frame()
# Initialize an empty dataframe to store R-squared values
all_r_squared <- data.frame()
# Define the range of metabolite columns
metabolite_columns <- 2:629 ##Always double check with column names
# Loop through each metabolite and run the regression
for (metabolite in colnames(WIDE_dataSub_with_metadata)[metabolite_columns]) {
  # Filter data for the current metabolite
  current_data <- WIDE_dataSub_with_metadata %>%
    dplyr::select(indiv.ID, MW_scores_lifestyle, Age, Sex, batch, PC1, all_of(metabolite)) %>%
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
    # Run the regression model including interaction terms
    model <- lm(Metabolite ~ MW_scores_lifestyle + Sex + Age + PC1 + batch + MW_scores_lifestyle:Sex, data = current_data)
    # Extract coefficients using broom
    temp_results <- broom::tidy(model) %>%
      mutate(
        Metabolite = metabolite,
        N = nrow(current_data) # Add the sample size for the regression
      )
    # Extract R-squared value
    model_r_squared <- broom::glance(model) %>%
      dplyr::select(r.squared) %>%
      mutate(Metabolite = metabolite)
    # Append results to the master dataframes
    all_coefficients <- bind_rows(all_coefficients, temp_results)
    all_r_squared <- bind_rows(all_r_squared, model_r_squared)
  } else {
    message(paste("Skipping metabolite:", metabolite, "- insufficient data after removing outliers"))
  }
}
#view(all_r_squared)
# Apply Bonferroni correction to p-values
all_coefficients <- all_coefficients %>%
  mutate(adjusted_p_values = p.adjust(p.value, method = "bonferroni"))
# Filter to find significant results
significant_results <- all_coefficients %>%
  filter(adjusted_p_values < 0.05) %>%
  arrange(adjusted_p_values)
# Remove unwanted terms like intercept, PC1, and batch
significant_results_no_intercept <- significant_results %>%
  filter(!term %in% c("(Intercept)", "PC1", "batch", "batch1", "batch2", "batch3", "batch4", "batch5.1", "batch5.2"))
# Summarize data by predictors (MW_scores_h_sol_trim, Age, Sex, interaction terms)
term_counts <- significant_results_no_intercept %>%
  group_by(term) %>%
  summarise(Count = n_distinct(Metabolite), .groups = 'drop')
# Create a bar plot for significant metabolites by predictor
p <- ggplot(term_counts, aes(x = reorder(term, Count), y = Count, fill = term)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  labs(
    title = "Model with lifestyle by Sex interaction",
    x = "Predictor",
    y = "Count of Significant Metabolites"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve label readability
# Print the plot
print(p)
# Save the bar plot
ggsave("Number_of_Signif_Metabs_by_Predictor.png", plot = p, width = 8, height = 6, dpi = 300)

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
  mutate(log2_fold_change = log2(exp(estimate)))  # Convert regression estimate to log2 fold change

# Define custom fold change threshold for significant metabolites
fc_threshold <- 0.67  

# Adjust volcano plot with log2 FC
volcano_data <- all_coefficients %>%
  filter(term == "MW_scores_lifestyleUrban") %>%
  mutate(
    log_p_value = -log10(adjusted_p_values),
    significant = adjusted_p_values < 0.05 & abs(log2_fold_change) > fc_threshold,
    regulation = case_when(
      log2_fold_change > fc_threshold & adjusted_p_values < 0.05 ~ "Upregulated in Urban",
      log2_fold_change < -fc_threshold & adjusted_p_values < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Volcano Plot with different colors for up/down regulation and labels for significant points
volcano_plot <- ggplot(volcano_data, aes(x = log2_fold_change, y = log_p_value)) +
  geom_point(aes(color = regulation), alpha = 0.8, size = 3) + 
  scale_color_manual(values = c("Upregulated in Urban" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
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
top_50_metabolites_sex <- significant_results_no_intercept %>%
  filter(term == "SexMale") %>%
  arrange(desc(abs(estimate))) %>%
  head(10)
# Create a bar plot for the top 50 metabolites by signed effect size for Sex with \(N\) displayed
top_50_plot_sex <- ggplot(top_50_metabolites_sex, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Display \(N\) above the bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites Affected by Sex",
    x = "Metabolite",
    y = "Beta value"
  ) +
  coord_flip() +
  theme_minimal()
# Print the plot for Sex
print(top_50_plot_sex)
# Save the plot for Sex
ggsave("Top_Metabolites_Affected_by_Sex.png", plot = top_50_plot_sex, width = 8, height = 6, dpi = 300)

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
ggsave("Interaction_Effect_Sex_H_sol.png", plot = interaction_plot, width = 10, height = 6, dpi = 300)

###
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


