##Set your working Directory
setwd("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/")
# Load necessary libraries
library(tidyr)
library(stringr)
library(ggplot2)
library(rio)
library(tidyverse)
library(MetaboAnalystR)
library(openxlsx)
library(tibble)
library(dplyr)
# Files Needed;
#1. Metabolomics Excel sheet from Asael
#2. Sample map
#3. Latest Frozen_Data by Kristina
# read in the HILIC results and remove the columns that are not samples or pooled QC samples
allData <- as.data.frame(import( "Turkana 613 samples Dual-neg+pos SUM.xlsx", sheet = "Organized"))
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
##ADD METADATA YOU WANT TO USE - Use Updated Frozen data Each time
All_data <- import("/Users/Bm0211/RegEx/Water_Energy/THGP_database_TurkanaOnly_corrected_merged_2025-01-13.txt")
dim(All_data)
##Choose whatever Metadata you want to work with given its a list of ~900 columns by now
colnames(All_data)
Meta_ONLY <- All_data %>% dplyr::select(c("Unique.ID", "Age", "Sex", "Sampling.location", "Tribe", "MW_scores_h_sol_trim", "MW_scores_lifestyle","Chapati.yes.no", "BMI", "Non.HDL.mg.dL", "LDL.HDL.ratio", "TC.HDL.ratio", "LDL.cholesterol.mg.dL", "Triglycerides.mg.dL", "HDL.cholesterol.mg.dL", "Total.cholesterol.mg.dL.", "Occupation", "Main.subsistence.activity"))
head(Meta_ONLY)
##sample Map
sampleMap <- as.data.frame(import("/Users/bm0211/Documents/AyrolesLab/metabolomics_sampleMap_final.csv"))
##Merge Map and Unique.ID data 
batchMap_with_Meta <- merge(sampleMap, batchMap, by=2, all = T)
batchMap_with_Meta <- unique(merge(batchMap_with_Meta, Meta_ONLY, by="Unique.ID", all.x = T))
dim(batchMap_with_Meta)
colnames(batchMap_with_Meta)
table(batchMap_with_Meta$MW_scores_h_sol_trim)
##the ordered metabolites file
dataSub.ordered = rownames_to_column(dataSub.ordered, var = "Metabolites")
dim(dataSub.ordered)
# Pivot and merge the data
NARROW_dataSub.ordered <- dataSub.ordered %>% 
  pivot_longer(!Metabolites,
               names_to = "indiv.ID",
               values_to = "Concentration")
# Join HILIC_Data with Meta_data
NARROW_dataSub.ordered <- left_join(NARROW_dataSub.ordered, batchMap_with_Meta, by = c('indiv.ID' = 'metabRunNum'))
# Median center and log transform - Exact same method suggested by Rabinowiztz lab (Mike)
NARROW_dataSub.ordered_Norm <- NARROW_dataSub.ordered %>%
  group_by(Metabolites) %>%
  mutate(med = median(Concentration, na.rm = TRUE)) %>%
  mutate(normIC = Concentration / med) %>%
  ungroup() %>%
  mutate(min.val = min(abs(normIC[normIC != 0]), na.rm = TRUE) / 10) %>%
  group_by(indiv.ID) %>%
  mutate(normLogIC = log10((normIC + sqrt(normIC^2 + min.val^2)) / 2))
# Remove unwanted rows based on Sex column
NARROW_dataSub.ordered_Norm <- NARROW_dataSub.ordered_Norm[!grepl("Female\\|Female\\|Male|Female\\|Male|Female\\|Male\\|Male", NARROW_dataSub.ordered_Norm$Sex), ]
# Separate metadata columns for later merge
metadata_columns <- NARROW_dataSub.ordered_Norm %>%
  dplyr::select(-Metabolites, -normLogIC) %>%  # Exclude columns used in pivoting
  distinct(indiv.ID, .keep_all = TRUE)  # Ensure one row per individual
# Pivot the data back to wide format using normLogIC
WIDE_dataSub_normLogIC <- NARROW_dataSub.ordered_Norm %>%
  pivot_wider(
    id_cols = indiv.ID,           # Unique identifier for individuals
    names_from = Metabolites,     # Use metabolites as column names
    values_from = normLogIC       # Use normLogIC for the values
  )
# Merge the metadata back into the wide dataset
WIDE_dataSub_with_metadata <- WIDE_dataSub_normLogIC %>%
  left_join(metadata_columns, by = "indiv.ID")

# Check the dimensions of the new wide dataset
dim(WIDE_dataSub_with_metadata)

# Define the outlier detection and marking function
mark_outliers <- function(x) {
  # Calculate the interquartile range (IQR)
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  # Define outlier bounds
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  # Identify outliers
  is_outlier <- x < lower_bound | x > upper_bound
  # Mark outliers with "||OUTLIER"
  x_marked <- ifelse(is_outlier, paste(x, "||OUTLIER"), x)
  return(x_marked)
}

# Apply the outlier marking function to metabolite columns (2:629)
WIDE_dataSub_with_metadata_clean <- WIDE_dataSub_with_metadata %>%
  mutate(across(2:629, ~ mark_outliers(.)))
# Save the final wide dataset to a CSV for further use
write.csv(WIDE_dataSub_with_metadata, "WIDE_dataSub_with_metadata.csv", row.names = FALSE)
#The above dataset is used for analysis on Metabo-analyst pathway enrichment 
# Select and clean up unnecessary columns - FOR THE DOWNSTREAM ANALYSIS FROM HERE 
WIDE_dataSub_with_metadata <- dplyr::select(WIDE_dataSub_with_metadata, -c(630, 632:637, 639, 642, 643, 646:659))

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
    select(indiv.ID, all_of(metabolite_name), MW_scores_h_sol_trim, Age, Sex, batch) %>%
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
    left_join(WIDE_dataSub_with_metadata %>% select(indiv.ID, batch), by = "indiv.ID") %>%
    drop_na() # Remove rows with missing values
  # Ensure only numeric columns are passed to PCA
  metabolite_residuals <- residuals_wide %>%
    select(-c(indiv.ID, batch)) %>%
    select_if(is.numeric)
  # Perform PCA on residuals
  pca_result <- prcomp(metabolite_residuals, center = TRUE, scale. = TRUE)
  # Extract PCA scores
  pca_scores <- as.data.frame(pca_result$x)
  pca_scores <- pca_scores %>%
    select(PC1, PC2, PC3) %>%
    mutate(indiv.ID = residuals_wide$indiv.ID)
  # Add PCA scores to the original dataset
  WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
    left_join(pca_scores, by = "indiv.ID")
  # Save the updated dataset with PCA scores --> This is to use with Metaboanalyst for regression and volcano plot
  write.csv(WIDE_dataSub_with_metadata, "WIDE_dataSub_with_metadata_with_PCs.csv", row.names = FALSE)
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
#Regression with The continuous lifestyle variable h_sol
WIDE_dataSub_with_metadata$MW_scores_h_sol_trim <- as.numeric(WIDE_dataSub_with_metadata$MW_scores_h_sol_trim)
##  For the regression, I ran one model at a time, 
## hence you have to switch which terms are in the interacting term of the model each time - do this manually
## This will allow you be sure which analysis you do
## MW_scores_h_sol_trim * Sex + Age + PC1 + batch
## MW_scores_h_sol_trim * Age + Sex + PC1 + batch
# Initialize an empty dataframe to store coefficients
all_coefficients <- data.frame()
# Define the range of metabolite columns
metabolite_columns <- 2:629 ##Always double check with column names
# Loop through each metabolite and run the regression
for (metabolite in colnames(WIDE_dataSub_with_metadata)[metabolite_columns]) {
  # Filter data for the current metabolite
  current_data <- WIDE_dataSub_with_metadata %>%
    select(indiv.ID, MW_scores_h_sol_trim, Age, Sex, batch, PC1, all_of(metabolite)) %>%
    rename(Metabolite = all_of(metabolite)) %>%
    filter(
      !is.na(Metabolite),
      !is.na(MW_scores_h_sol_trim),
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
    model <- lm(Metabolite ~ MW_scores_h_sol_trim * Sex + Age + PC1 + batch, data = current_data)
    # Extract coefficients using broom
    temp_results <- broom::tidy(model) %>%
      mutate(
        Metabolite = metabolite,
        N = nrow(current_data) # Add the sample size for the regression
      )
    # Append results to the master dataframe
    all_coefficients <- bind_rows(all_coefficients, temp_results)
  } else {
    message(paste("Skipping metabolite:", metabolite, "- insufficient data after removing outliers"))
  }
}
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
    title = "Count of Signif Metabs by Predictor with Sex_X_H_sol_interaction",
    x = "Predictor",
    y = "Count of Significant Metabolites"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve label readability
# Print the plot
print(p)
# Save the bar plot
ggsave("Number_of_Signif_Metabs_by_Predictor_model_with_Sex_X_H_sol_interaction.png", plot = p, width = 8, height = 6, dpi = 300)

# Identify top 50 metabolites by effect size for MW_scores_h_sol_trim
top_50_metabolites_h_sol <- significant_results %>%
  filter(term == "MW_scores_h_sol_trim") %>%
  arrange(desc(abs(estimate))) %>%
  head(50)
# Create a bar plot for the top 50 metabolites by signed effect size with \(N\) on top
top_50_plot_signed <- ggplot(top_50_metabolites_h_sol, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) + # Add \(N\) as text on top of bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites by effect size MW_scores_h_sol_trim_with_sex_interaction",
    x = "Metabolite",
    y = "Effect Size"
  ) +
  coord_flip() +
  theme_minimal()
# Print the top 50 plot
print(top_50_plot_signed)
# Save the plot
ggsave("Top_Metabolites_Affected_by_MW_scores_h_sol_trim_with_sex_interaction.png", plot = top_50_plot_signed, width = 8, height = 6, dpi = 300)

# Identify top 50 metabolites by effect size for Age
top_50_metabolites_age <- significant_results_no_intercept %>%
  filter(term == "Age") %>%
  arrange(desc(abs(estimate))) %>%
  head(50)
# Create a bar plot for the top 50 metabolites by signed effect size for Age with \(N\) displayed
top_50_plot_age <- ggplot(top_50_metabolites_age, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Display \(N\) above the bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites Affected by Age_Model_with sex_X_H_sol_interaction",
    x = "Metabolite",
    y = "Effect Size"
  ) +
  coord_flip() +
  theme_minimal()
# Print the plot for Age
print(top_50_plot_age)
# Save the plot for Age
ggsave("Top_Metabolites_Affected_by_Age_Model_with_sex_interaction.png", plot = top_50_plot_age, width = 8, height = 6, dpi = 300)

# Identify top 50 metabolites by effect size for Sex
top_50_metabolites_sex <- significant_results_no_intercept %>%
  filter(term == "SexMale") %>%
  arrange(desc(abs(estimate))) %>%
  head(50)
# Create a bar plot for the top 50 metabolites by signed effect size for Sex with \(N\) displayed
top_50_plot_sex <- ggplot(top_50_metabolites_sex, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Display \(N\) above the bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites Affected by Sex_Model_with_sex_X_H_sol_interaction",
    x = "Metabolite",
    y = "Effect Size"
  ) +
  coord_flip() +
  theme_minimal()
# Print the plot for Sex
print(top_50_plot_sex)
# Save the plot for Sex
ggsave("Top_Metabolites_Affected_by_Sex_model_with_sex_interaction.png", plot = top_50_plot_sex, width = 8, height = 6, dpi = 300)

# Identify metabolites by effect size for SEX_H_SOL INteraction
table(significant_results_no_intercept$term)
Interact_metabs <- significant_results_no_intercept %>%
  filter(term == "MW_scores_h_sol_trim:SexMale") %>%
  arrange(desc(abs(estimate))) %>%
  head(50)
# Create a bar plot for the top 50 metabolites by signed effect size for Age with \(N\) displayed
Interact_metabs_plot <- ggplot(Interact_metabs, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Display \(N\) above the bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Metabolites Affected by Sex__H_sol_interaction",
    x = "Metabolite",
    y = "Effect Size"
  ) +
  coord_flip() +
  theme_minimal()
# Print the plot for Age
print(Interact_metabs_plot)
# Save the plot for Age
ggsave("Metabolites Affected by Sex__H_sol_interaction.png", plot = Interact_metabs_plot, width = 8, height = 6, dpi = 300)
###
#########
#############
#################
############
########
#Working with Lifestyle as dicrete
table(WIDE_dataSub_with_metadata$MW_scores_lifestyle)
# Combine Pastoralist and PeriUrban into Non_Urban
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  mutate(
    MW_scores_lifestyle = case_when(
      MW_scores_lifestyle %in% c("Pastoralist", "PeriUrban") ~ "Non_Urban",
      MW_scores_lifestyle == "Urban" ~ "Urban",
      TRUE ~ as.character(MW_scores_lifestyle) # Ensure no loss of other categories
    )
  )
table(WIDE_dataSub_with_metadata$MW_scores_lifestyle)
# Relevel MW_scores_lifestyle to make Urban the reference group
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  mutate(MW_scores_lifestyle = factor(MW_scores_lifestyle, levels = c("Urban", "Non_Urban")))

#Regression
#
## Initialize an empty dataframe to store coefficients
all_coefficients <- data.frame()
# Loop through each metabolite (columns 2:630) and run the regression
for (metabolite in colnames(WIDE_dataSub_with_metadata)[2:629]) {
  # Filter data for the current metabolite
  current_data <- WIDE_dataSub_with_metadata %>%
    select(indiv.ID, MW_scores_lifestyle, Age, Sex, batch, PC1, all_of(metabolite)) %>%
    rename(Metabolite = all_of(metabolite)) %>%
    filter(
      !is.na(Metabolite),
      !is.na(MW_scores_lifestyle),
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
  # Ensure enough data to fit the model
  if (nrow(current_data) > 6) {
    # Convert categorical variables to factors
    current_data <- current_data %>%
      mutate(
        batch = as.factor(batch),
        MW_scores_lifestyle = as.factor(MW_scores_lifestyle)
      )
    # Run the regression model including interaction terms
    model <- lm(Metabolite ~ MW_scores_lifestyle * Sex + Age + PC1 + batch, data = current_data)
    # Extract coefficients using broom
    temp_results <- broom::tidy(model) %>%
      mutate(
        Metabolite = metabolite,
        N = nrow(current_data) # Add sample size (N)
      )
    # Append results to the master dataframe
    all_coefficients <- bind_rows(all_coefficients, temp_results)
  }
}
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
table(significant_results_no_intercept$term)
colnames(significant_results_no_intercept)

# Plot for Sex
top_50_metabolites_sex <- significant_results_no_intercept %>%
  filter(term == "SexMale") %>%
  arrange(desc(abs(estimate))) %>%
  head(50)
top_50_plot_sex <- ggplot(top_50_metabolites_sex, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Add N labels above bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites Affected by Sex",
    x = "Metabolite",
    y = "Effect Size"
  ) +
  coord_flip() +
  theme_minimal()
print(top_50_plot_sex)
ggsave("Top_Metabolites_Affected_by_Sex.png", plot = top_50_plot_sex, width = 8, height = 6, dpi = 300)

# Plot for Age
top_50_metabolites_age <- significant_results_no_intercept %>%
  filter(term == "Age") %>%
  arrange(desc(abs(estimate))) %>%
  head(50)
top_50_plot_age <- ggplot(top_50_metabolites_age, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Add N labels above bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites Affected by Age",
    x = "Metabolite",
    y = "Effect Size"
  ) +
  coord_flip() +
  theme_minimal()
print(top_50_plot_age)
ggsave("Top_Metabolites_Affected_by_Age.png", plot = top_50_plot_age, width = 8, height = 6, dpi = 300)

# Filter and reverse coefficients for "Urban" effects
top_50_metabolites_lifestyle <- significant_results_no_intercept %>%
  filter(term == "MW_scores_lifestyleNon_Urban") %>%  # Use the existing term
  mutate(estimate = -estimate) %>%  # Reverse the sign of the estimate
  arrange(desc(abs(estimate))) %>%
  head(50)
# Create the bar plot
top_50_plot_lifestyle <- ggplot(top_50_metabolites_lifestyle, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Add N labels above bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites Affected by Lifestyle (Urban vs Non_Urban)",
    x = "Metabolite",
    y = "Effect Size (Urban vs Non_Urban)"
  ) +
  coord_flip() +
  theme_minimal()
# Print and save the plot
print(top_50_plot_lifestyle)
ggsave("Top_Metabolites_Affected_by_Lifestyle_Urban_Adjusted.png", plot = top_50_plot_lifestyle, width = 8, height = 6, dpi = 300)

# Plot for Interaction (Sex x Lifestyle)
top_50_metabolites_interaction <- significant_results_no_intercept %>%
  filter(term == "MW_scores_lifestyleNon_Urban:SexMale") %>%
  arrange(desc(abs(estimate))) %>%
  head(50)
top_50_plot_interaction <- ggplot(top_50_metabolites_interaction, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +  # Add N labels above bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(
    title = "Top Metabolites Affected by Lifestyle × Sex Interaction",
    x = "Metabolite",
    y = "Effect Size"
  ) +
  coord_flip() +
  theme_minimal()
print(top_50_plot_interaction)
ggsave("Top_Metabolites_Affected_by_sex_lifestyle_Interaction.png", plot = top_50_plot_interaction, width = 8, height = 6, dpi = 300)

# Predictor term count plot 
term_counts <- significant_results_no_intercept %>%
  group_by(term) %>%
  summarise(Count = n_distinct(Metabolite), .groups = 'drop') %>%
  complete(term = unique(significant_results_no_intercept$term), fill = list(Count = 0))

p <- ggplot(term_counts, aes(x = reorder(term, Count), y = Count, fill = term)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  labs(
    title = "Count of Signif Metabos by Model Terms (model_with_Age_lifestyle_interaction)",
    x = "Model Term",
    y = "Count of Significant Metabolites"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
ggsave("Number_of_Significant_Metabolites_by_Model_Terms_model_with_lifestyle_X_Sex_interaction.png", plot = p, width = 8, height = 6, dpi = 300)


##FIGURING OUT THE COMMON ONES ACROSS EVERYTHING
##Run all the Models and pick the common ones across models 
# Initialize empty dataframes to store results for lifestyle and sex
all_lifestyle_results <- list()
all_sex_results <- list()
# Define models
##Filter out the h_sol * Age model because this one has few ones affected by h_sol
models <- list(
  "Lifestyle_Discrete_Sex" = "Metabolite ~ MW_scores_lifestyle * Sex + Age + PC1 + batch",
  "Lifestyle_Discrete_Age" = "Metabolite ~ MW_scores_lifestyle * Age + Sex + PC1 + batch",
  #"Lifestyle_Continuous_Age" = "Metabolite ~ MW_scores_h_sol_trim * Age + Sex + PC1 + batch",
  "Lifestyle_Continuous_Sex" = "Metabolite ~ MW_scores_h_sol_trim * Sex + Age + PC1 + batch"
)
# Function to remove outliers using IQR
remove_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  data <- data %>%
    filter(.data[[column]] >= lower_bound & .data[[column]] <= upper_bound)
  return(data)
}
# Loop through each model
for (model_name in names(models)) {
  # Initialize an empty dataframe to store coefficients
  coefficients <- data.frame()
  for (metabolite in colnames(WIDE_dataSub_with_metadata)[2:629]) {
    # Filter data for the current metabolite
    current_data <- WIDE_dataSub_with_metadata %>%
      select(indiv.ID, MW_scores_lifestyle, MW_scores_h_sol_trim, Age, Sex, batch, PC1, all_of(metabolite)) %>%
      rename(Metabolite = all_of(metabolite)) %>%
      filter(
        !is.na(Metabolite),
        !is.na(MW_scores_lifestyle),
        !is.na(MW_scores_h_sol_trim),
        !is.na(Age),
        !is.na(Sex),
        !is.na(PC1)
      )
    # Remove outliers
    if (nrow(current_data) > 0) {
      current_data <- remove_outliers(current_data, "Metabolite")
    }
    if (nrow(current_data) > 6) {
      # Convert categorical variables to factors
      current_data <- current_data %>%
        mutate(
          batch = as.factor(batch),
          MW_scores_lifestyle = as.factor(MW_scores_lifestyle)
        )
      # Fit the regression model
      model <- lm(as.formula(models[[model_name]]), data = current_data)
      # Extract coefficients
      temp_results <- broom::tidy(model) %>%
        mutate(Metabolite = metabolite)
      coefficients <- bind_rows(coefficients, temp_results)
    }
  }
  
  # Apply Bonferroni correction and filter significant results
  coefficients <- coefficients %>%
    mutate(adjusted_p_values = p.adjust(p.value, method = "bonferroni")) %>%
    filter(adjusted_p_values < 0.05)
  # Separate top 50 for lifestyle and sex
  all_lifestyle_results[[model_name]] <- coefficients %>%
    filter(str_detect(term, "MW_scores")) %>%
    arrange(desc(abs(estimate))) %>%
    head(50)
  all_sex_results[[model_name]] <- coefficients %>%
    filter(term == "SexMale") %>%
    arrange(desc(abs(estimate))) %>%
    head(50)
}

# Identify common metabolites across all models
common_lifestyle <- Reduce(intersect, lapply(all_lifestyle_results, function(x) x$Metabolite))
common_sex <- Reduce(intersect, lapply(all_sex_results, function(x) x$Metabolite))
# Function to remove outliers using IQR
remove_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  data <- data %>%
    filter(.data[[column]] >= lower_bound & .data[[column]] <= upper_bound)
  return(data)
}
# Function to save plots in 3x3 grids over multiple pages
save_violin_plots <- function(plots, filename, ncol = 3, nrow = 3) {
  pdf(filename, width = 11, height = 8.5)
  for (i in seq(1, length(plots), by = ncol * nrow)) {
    gridExtra::grid.arrange(
      grobs = plots[i:min(i + (ncol * nrow - 1), length(plots))],
      ncol = ncol
    )
  }
  dev.off()
}

# Generate violin plots for common lifestyle metabolites
if (length(common_lifestyle) > 0) {
  lifestyle_plots <- lapply(common_lifestyle, function(metabolite) {
    # Remove outliers for the specific metabolite
    data <- WIDE_dataSub_with_metadata %>%
      filter(!is.na(MW_scores_lifestyle), !is.na(!!sym(metabolite))) %>%
      remove_outliers(column = metabolite)  # Use metabolite name directly
    # Calculate sample size (N) for each category
    sample_sizes <- data %>%
      group_by(MW_scores_lifestyle) %>%
      summarise(N = n(), .groups = "drop")
    # Create violin plot
    ggplot(data, aes(x = MW_scores_lifestyle, y = !!sym(metabolite))) +
      geom_violin(trim = FALSE, aes(fill = MW_scores_lifestyle), alpha = 0.5) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.8) +
      geom_text(data = sample_sizes, aes(x = MW_scores_lifestyle, y = Inf, label = paste0("N=", N)),
                vjust = 1.5, size = 3, inherit.aes = FALSE) + # Add N to the plot
      labs(title = metabolite, x = "Lifestyle", y = "Abundance") +
      theme_minimal()
  })
  save_violin_plots(lifestyle_plots, "Lifestyle_Metabolites_Violin_Plots.pdf")
}

# Generate violin plots for common sex metabolites
if (length(common_sex) > 0) {
  sex_plots <- lapply(common_sex, function(metabolite) {
    # Remove outliers for the specific metabolite
    data <- WIDE_dataSub_with_metadata %>%
      filter(!is.na(Sex), !is.na(!!sym(metabolite))) %>%
      remove_outliers(column = metabolite)  # Use metabolite name directly
    # Calculate sample size (N) for each category
    sample_sizes <- data %>%
      group_by(Sex) %>%
      summarise(N = n(), .groups = "drop")
    # Create violin plot
    ggplot(data, aes(x = Sex, y = !!sym(metabolite))) +
      geom_violin(trim = FALSE, aes(fill = Sex), alpha = 0.5) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.8) +
      geom_text(data = sample_sizes, aes(x = Sex, y = Inf, label = paste0("N=", N)),
                vjust = 1.5, size = 3, inherit.aes = FALSE) + # Add N to the plot
      labs(title = metabolite, x = "Sex", y = "Abundance") +
      theme_minimal()
  })
  save_violin_plots(sex_plots, "Sex_Metabolites_Violin_Plots.pdf")
}

library(ggplot2)
library(dplyr)

# Filter the dataset for Ginsenoside Rh3 isomer and non-missing lifestyle data
filtered_data <- WIDE_dataSub_with_metadata %>%
  filter(!is.na(MW_scores_lifestyle), !is.na(`Ginsenoside Rh3 isomer`))

# Remove outliers for Ginsenoside Rh3 isomer using IQR
remove_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  data <- data %>%
    mutate(
      !!sym(column) := ifelse(
        data[[column]] < lower_bound | data[[column]] > upper_bound,
        NA,
        data[[column]]
      )
    )
  return(data)
}

filtered_data <- remove_outliers(filtered_data, "Ginsenoside Rh3 isomer")

# Calculate sample sizes for each lifestyle category
sample_sizes <- filtered_data %>%
  group_by(MW_scores_lifestyle) %>%
  summarise(N = n(), .groups = "drop")

# Create the violin plot with jitter points
violin_plot <- ggplot(filtered_data, aes(x = MW_scores_lifestyle, y = `Ginsenoside Rh3 isomer`)) +
  geom_violin(trim = FALSE, aes(fill = MW_scores_lifestyle), alpha = 0.5) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.8) +
  geom_text(data = sample_sizes, aes(x = MW_scores_lifestyle, y = Inf, label = paste0("N=", N)),
            vjust = 1.5, size = 3, inherit.aes = FALSE) + # Add sample sizes to the plot
  labs(title = "Ginsenoside Rh3 Isomer Abundance by Lifestyle",
       x = "Lifestyle (Urban vs Non-Urban)",
       y = "Abundance") +
  theme_minimal()

# Save the plot as a PDF
ggsave("Ginsenoside_Rh3_isomer_violin_plot.pdf", plot = violin_plot, width = 8, height = 6)

# Display the plot
print(violin_plot)
