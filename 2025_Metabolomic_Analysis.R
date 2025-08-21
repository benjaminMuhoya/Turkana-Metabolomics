##Set your working Directory
setwd("/Users/bm0211/RegEx/Water_Energy/HILIC_ms_analysis/")
# Load necessary libraries
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(rio)
library(tidyverse)
library(tibble)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(broom.mixed)
library(lme4)
library(lmerTest)  # Enables p-values for mixed models
library(readr)
library(modelr)    
library(broom)      
library(patchwork)
library(caret)
library(pls)
# Files Needed;
#1. Metabolomics Excel sheet from Asael
#2. Sample map
#3. Latest Frozen_Data by Kristina
#4. List of Unique_ID for which we have DNA
#5. Key to match Metabolite to HMDB ascension number
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
colnames(All_data)
dim(All_data)
table(All_data$MW_scores_market_diet_index)
table(All_data$Coffee.yes.no)
table(All_data$MW_scores_h_sol)
## Choose whatever Metadata you want to work with
Meta_ONLY <- All_data %>% dplyr::select(c("Unique.ID","MW_scores_h_sol_trim", "Sampling.location","MW_scores_lifestyle", "Main.subsistence.activity", "Age", "Sex", "MW_scores_market_diet_index", "TA_score_TL_108", "MW_scores_h_sol"))
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
table(batchMap_with_Meta$MW_scores_lifestyle)
## Ordered metabolites file
dataSub.ordered = rownames_to_column(dataSub.ordered, var = "Metabolites")
dim(dataSub.ordered)
##########################################################################################################
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
head(WIDE_dataSub)
# Merge the metadata back into the wide dataset
WIDE_dataSub_with_metadata <- WIDE_dataSub %>%
  left_join(metadata_columns, by = "indiv.ID")
# Check dimensions of the new wide dataset
dim(WIDE_dataSub_with_metadata)
colnames(WIDE_dataSub_with_metadata)
head(WIDE_dataSub_with_metadata)
# Remove duplicates in Unique.ID, keeping the first occurrence
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  distinct(Unique.ID, .keep_all = TRUE)
table(is.na(WIDE_dataSub_with_metadata$MW_scores_lifestyle))
# Drop row with indiv.ID "M116" ...This is an outlier identified using PCA
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  filter(indiv.ID != "M116")
##########################################################################################################
# Grouping lifestyles
##The Marina Score Below
WIDE_dataSub_with_metadata <- WIDE_dataSub_with_metadata %>%
  mutate(MW_scores_lifestyle_2 = ifelse(MW_scores_lifestyle == "Urban", "Urban", "Non_Urban")) %>%
  filter(!is.na(MW_scores_lifestyle))
table(is.na(WIDE_dataSub_with_metadata$MW_scores_lifestyle))
table(WIDE_dataSub_with_metadata$MW_scores_lifestyle)
table(WIDE_dataSub_with_metadata$MW_scores_lifestyle_2)
colnames(WIDE_dataSub_with_metadata)
# Save the final wide dataset to a CSV
For_writing <- dplyr::select(WIDE_dataSub_with_metadata, c(630:642,1:629))


write.csv(For_writing, "HILIC_RAW_DATA_with_covariates.csv", row.names = FALSE)

# Read the Updated_HMDB_NAMES file # Download from Google Sheets
updated_hmdb_names <- read.csv("Updated_No_match_Metabo.csv", stringsAsFactors = FALSE)
head(updated_hmdb_names)
table(unique(updated_hmdb_names$HMDB))
# Print only the first 20 duplicated HMDB entries
updated_hmdb_names %>% 
  filter(HMDB %in% HMDB[duplicated(HMDB)]) %>%
  head(20) %>%
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
table(For_writing$MW_scores_lifestyle_2)
# Save the final wide dataset with HMDB names
write.csv(For_writing, "HILIC_HMDB_colnames_RAWdata_with_covs.csv", row.names = FALSE)

# Pivot to long format for median centering and log transformation
NARROW_dataSub.ordered <- WIDE_dataSub_with_metadata %>%
  pivot_longer(
    cols = -c(indiv.ID, Unique.ID, MW_scores_h_sol_trim, MW_scores_lifestyle, MW_scores_lifestyle_2, Main.subsistence.activity, Age, Sex, Run.ID, batch, Sampling.location, TA_score_TL_108),
    names_to = "Metabolites",
    values_to = "Concentration"
  )
head(NARROW_dataSub.ordered)
##########################################################################################################
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
# Pivot back to wide format using the normalized column (norm_LocIC) as the concentration
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
Second_write <- dplyr::select(WIDE_dataSub_with_metadata, c(630:642, 1:629))
# Save the transformed dataset
write.csv(Second_write, "HILIC_DATA_NormLog_plus_Metadata.csv", row.names = FALSE)
# Create datasets for each lifestyle
Non_Urban_data <- subset(Second_write, MW_scores_lifestyle %in% c("PeriUrban", "Pastoralist"))
Urban_data <- subset(Second_write, MW_scores_lifestyle == "Urban")
Pastoralist_data <- subset(Second_write, MW_scores_lifestyle == "Pastoralist")
PeriUrban_data <- subset(Second_write, MW_scores_lifestyle == "PeriUrban")
# Save datasets as CSV files
write.csv(Non_Urban_data, "HILIC_DATA_NormLog_Non_Urban_data.csv", row.names = FALSE)
write.csv(Urban_data, "HILIC_DATA_NormLog_Urban_data.csv", row.names = FALSE)
write.csv(Pastoralist_data, "HILIC_DATA_NormLog_Pastoralist_data.csv", row.names = FALSE)
write.csv(PeriUrban_data, "HILIC_DATA_NormLog_PeriUrban_data.csv", row.names = FALSE)
# Rename matching columns in Second_write
colnames(Second_write) <- ifelse(
  colnames(Second_write) %in% names(rename_vector),
  rename_vector[colnames(Second_write)],
  colnames(Second_write)  # Keep original name if no match
)
# Check renamed columns
colnames(Second_write)
table(Second_write$MW_scores_lifestyle)
table(Second_write$MW_scores_lifestyle_2)
# Save the final wide dataset with HMDB names
write.csv(Second_write, "HILIC_HMDB_col_names_NormLog_plus_Metadata.csv", row.names = FALSE)
# Create datasets for each lifestyle
Non_Urban_data <- subset(Second_write, MW_scores_lifestyle %in% c("PeriUrban", "Pastoralist"))
Urban_data <- subset(Second_write, MW_scores_lifestyle == "Urban")
Pastoralist_data <- subset(Second_write, MW_scores_lifestyle == "Pastoralist")
PeriUrban_data <- subset(Second_write, MW_scores_lifestyle == "PeriUrban")
# Save datasets as CSV files
write.csv(Non_Urban_data, "HILIC_HMDB_col_names_NormLog_Non_Urban_data.csv", row.names = FALSE)
write.csv(Urban_data, "HILIC_HMDB_col_names_NormLog_Urban_data.csv", row.names = FALSE)
write.csv(Pastoralist_data, "HILIC_HMDB_col_names_NormLog_Pastoralist_data.csv", row.names = FALSE)
write.csv(PeriUrban_data, "HILIC_HMDB_col_names_NormLog_PeriUrban_data.csv", row.names = FALSE)

#The above datasets were made for use online in the Metabo-analyst pathway enrichment 
################################################################################################################
## Below, I process the data further and do my analyses. 

# Select and clean up unnecessary columns - FOR THE DOWNSTREAM ANALYSIS FROM HERE 
#for_kristina <- dplyr::select(WIDE_dataSub_with_metadata, -c(643:646))
#colnames(WIDE_dataSub_with_metadata)
#dim(WIDE_dataSub_with_metadata)
#for_kristina <- dplyr::select(WIDE_dataSub_with_metadata, c(Unique.ID))
#for_kristina$passedProcessing <- 1
#for_kristina$passedQC <- NA
#head(for_kristina)
#dim(for_kristina)
#write.csv(for_kristina, "HILIC_First_Round_2022.csv", row.names = FALSE)
################################################################################################################
# The goal is to residualize the data, to remove batch effects
# 1. Create residuals per metabolite after adjusting for batch
residuals_long <- WIDE_dataSub_with_metadata %>% 
  pivot_longer(cols = 2:629,  # Adjust this if your metabolite columns differ
               names_to = "Metabolite",
               values_to = "Abundance") %>% 
  drop_na(Abundance, batch) %>% group_by(Metabolite) %>% 
  filter(all(is.finite(Abundance)) && var(Abundance) > 0) %>%
  mutate(Residual = add_residuals(cur_data(),
      lm(Abundance ~ Age + Sex + batch, data = cur_data())
    )$resid) %>% ungroup() %>% 
  dplyr::select(indiv.ID, Metabolite, Residual)

# 2. Pivot to wide format
residuals_wide <- residuals_long %>% pivot_wider(names_from = Metabolite, values_from = Residual)
# 3. Select metadata columns
metadata_columns <- c(
  "indiv.ID", "MW_scores_market_diet_index", "Unique.ID", "batch", "Run.ID",
  "MW_scores_h_sol_trim", "Sampling.location", "MW_scores_lifestyle", "MW_scores_h_sol",
  "Main.subsistence.activity", "Age", "Sex", "MW_scores_lifestyle_2", "TA_score_TL_108"
)
metadata <- WIDE_dataSub_with_metadata %>% dplyr::select(all_of(metadata_columns))
dim(residuals_wide)
dim(residuals_long)
dim(WIDE_dataSub_with_metadata)
# 4. Merge residuals with metadata
residualized_HILIC <- residuals_wide %>% left_join(metadata, by = "indiv.ID")
dim(residualized_HILIC)
table(residualized_HILIC$MW_scores_lifestyle_2)
table(residualized_HILIC$MW_scores_lifestyle)
# 5. Write final output
write_csv(residualized_HILIC, "residualized_HILIC.csv")
# Ensure your rename vector is already created:
# rename_vector <- setNames(updated_hmdb_names$HMDB, updated_hmdb_names$Query)
# rename_vector <- rename_vector[!is.na(rename_vector)]
# Identify metabolite columns (assumes they come before the metadata columns)
# Adjust range if needed
#metabolite_cols <- colnames(residualized_HILIC)[!colnames(residualized_HILIC) %in% metadata_columns]

# Rename only metabolite columns using HMDB names
#new_colnames <- colnames(residualized_HILIC)
#new_colnames[match(metabolite_cols, new_colnames)] <- ifelse(
 # metabolite_cols %in% names(rename_vector),
 # rename_vector[metabolite_cols],
 # metabolite_cols
#)
#colnames(residualized_HILIC) <- new_colnames

# Check if renaming worked
head(colnames(residualized_HILIC), 20)

# Write the renamed file
#write_csv(residualized_HILIC, "residualized_HILIC_HMDB_colnames_ONLY_LIFESTYLE_EFFECT_LEFT.csv")

################################################################################################################
# PCA on the residuals of models 
## Model; Metabolite ~ MW_scores_h_sol_trim + Age + Sex 
## Model: Metabolite ~ batch
# ----------regress each metabolite & collect residuals ----------
# long format → regress → residuals in one line per individual–metabolite
   # for add_residuals()
residuals_long <- WIDE_dataSub_with_metadata %>% 
  pivot_longer(cols = 2:629,
               names_to = "Metabolite",
               values_to = "Abundance") %>% 
  drop_na(batch) %>%       # keep complete cases
  group_by(Metabolite) %>% 
  filter(all(is.finite(Abundance)) && var(Abundance) > 0) %>%
  mutate(
    Residual = add_residuals(cur_data(),
      lm(Abundance ~ batch, data = cur_data())
    )$resid) %>% ungroup() %>% 
  dplyr::select(indiv.ID, Metabolite, Residual)
# ----------pivot wide, run PCA, merge PCs back----------
# Pivot + Join metadata AFTER residuals
residuals_wide <- residuals_long %>%
  pivot_wider(names_from = Metabolite, values_from = Residual) %>%
  left_join(WIDE_dataSub_with_metadata %>%
              dplyr::select(indiv.ID, batch, MW_scores_lifestyle_2, Sex,
                            MW_scores_lifestyle, MW_scores_h_sol_trim, TA_score_TL_108, MW_scores_h_sol),
            by = "indiv.ID")

dim(residuals_wide)
dim(residuals_long)
dim(WIDE_dataSub_with_metadata)
# keep only residual columns for PCA
metab_mat <- residuals_wide %>% dplyr::select(where(is.numeric), -batch, -MW_scores_h_sol_trim, -TA_score_TL_108, -MW_scores_h_sol)
# Run PCA
pca_res   <- prcomp(metab_mat, center = TRUE, scale. = TRUE)
# pack PC1–PC3 scores
pca_scores <- as_tibble(pca_res$x[, 1:3]) %>% mutate(indiv.ID = residuals_wide$indiv.ID) %>% 
  rename(PC1 = 1, PC2 = 2, PC3 = 3)
head(pca_scores)

table(residuals_wide$Sex)
# add PCs to master data
HILIC_with_PCs <- WIDE_dataSub_with_metadata %>% left_join(pca_scores, by = "indiv.ID")
##Writting
write_csv(HILIC_with_PCs, "HILIC_data_with_metadata_and_PCs.csv")
# ----------  make + save PCA plot ----------
var_exp <- round(summary(pca_res)$importance["Proportion of Variance", 1:2] * 100, 2)

# PCA plot
residuals_wide$MW_scores_h_sol_trim <- as.factor(residuals_wide$MW_scores_h_sol_trim)
# Define plotting variables
factor_vars <- c("MW_scores_lifestyle_2", "MW_scores_lifestyle", "MW_scores_h_sol_trim", "batch", "Sex")
continuous_vars <- c("TA_score_TL_108", "MW_scores_h_sol")
# Color palettes
gradient_colors <- c("#29AB87", "#FFC067", "#FF0000")
wes_colors <- c("1" = "#FF0000", "Urban" = "#FF0000", "2"= "#E1AD01", "Non_Urban"= "#29AB87", "Pastoralist"= "#29AB87","3"= "#29AB87","4"= "#4682B4", "5.1" = "#C1A192", "5" = "#C1A192", "Male" = "#C1A192", "5.2" = "#9A4EAE", "6" = "#9A4EAE","Female" = "#9A4EAE", "PeriUrban" = "#FFC067")
# Loop over variables and plot
for (var in c(factor_vars, continuous_vars)) {
  df_plot <- pca_scores %>%
    left_join(residuals_wide %>% dplyr::select(indiv.ID, !!sym(var)), by = "indiv.ID")
  p <- ggplot(df_plot, aes(PC1, PC2, color = .data[[var]])) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = indiv.ID), size = 3, max.overlaps = 20) +
    labs(title = paste("PCA colored by", var),
      x = paste0("PC1 (", var_exp[1], "% var)"),
      y = paste0("PC2 (", var_exp[2], "% var)"),
      color = var) + theme_minimal(base_size = 12)
  # Handle color scale
  if (var %in% factor_vars) {
    df_plot[[var]] <- as.factor(df_plot[[var]])
    p <- p + scale_color_manual(values = wes_colors, na.value = "black")
  } else {
    p <- p + scale_color_gradientn(colors = gradient_colors, na.value = "black")
  }
  print(p)
}
#ggsave("PCA_metabolite_residuals.png", plot = p, width = 8, height = 6, dpi = 300)
################################################################################################################
################################################################################################################
# Step 1: Extract PC2 loadings
loadings_pc2 <- pca_res$rotation[, "PC2"]
# Step 2: Invert loadings for visualization
pc2_df <- data.frame(Metabolite = names(loadings_pc2),
  Loading = -1 * loadings_pc2) %>% mutate(Abs_Loading = abs(Loading),
    Group = case_when(Loading > 0 ~ "Urban",Loading < 0 ~ "Non_Urban"))
# Step 3: Get top 30 contributors (by absolute value)
top_pc2_df <- pc2_df %>% arrange(desc(Abs_Loading)) %>% slice_head(n = 35)
# Step 4: Plot
ggplot(top_pc2_df, aes(x = reorder(Metabolite, Loading), y = Loading, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Urban" = "#FD6467", "Non_Urban" = "#798234")) +
  coord_flip() + theme_minimal(base_size = 13) +
  labs( title = "Top Metabolites PC2 Loadings",
    x = "Metabolite", y = "PC2 Loading", fill = "More Abundant In")
################################################################################################################
################################################################################################################
######
residualized_HILIC$MW_scores_lifestyle_2 <- as.factor(residualized_HILIC$MW_scores_lifestyle_2)
table(residualized_HILIC$MW_scores_lifestyle_2)
colnames(residualized_HILIC)
# Initialize an empty dataframe to store coefficients
all_coefficients <- data.frame()
all_r_squared <- data.frame()
metabolite_columns <- 2:610  # you know the drill
# Loop over metabolites
for (metabolite in colnames(residualized_HILIC)[metabolite_columns]) {
  message(paste0("🔍 Processing: ", metabolite))
  current_data <- residualized_HILIC %>%
    dplyr::select(indiv.ID, MW_scores_lifestyle, MW_scores_lifestyle_2, Age, Sex, batch, 
                  Sampling.location, all_of(metabolite)) %>%
    rename(Metabolite = all_of(metabolite)) %>%
    filter(!is.na(Metabolite), !is.na(Age), !is.na(Sex))
  # Outlier removal (IQR method)
  if (nrow(current_data) > 0) {
    Q1 <- quantile(current_data$Metabolite, 0.25, na.rm = TRUE)
    Q3 <- quantile(current_data$Metabolite, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    current_data <- current_data %>%
      filter(Metabolite >= (Q1 - 1.5 * IQR), Metabolite <= (Q3 + 1.5 * IQR))
  }
  # Proceed if enough valid data remains
  if (nrow(current_data) > 6 && is.numeric(current_data$Metabolite) && var(current_data$Metabolite, na.rm = TRUE) > 0) {
    # Wrap modeling in tryCatch
    tryCatch({
      model <- lm(Metabolite ~ MW_scores_lifestyle_2 + Sex + Age , data = current_data)
      temp_results <- broom::tidy(model) %>% mutate(Metabolite = metabolite, N = nrow(current_data))
      model_r_squared <- data.frame(r.squared = summary(model)$r.squared,
        r.squared_adj = summary(model)$adj.r.squared,
        Metabolite = metabolite)
      all_coefficients <- bind_rows(all_coefficients, temp_results)
      all_r_squared <- bind_rows(all_r_squared, model_r_squared)
      message(paste0("✅ Model succeeded for: ", metabolite))
    }, error = function(e) {
      message(paste("❌ Skipping", metabolite, "- model failed with error:", e$message))
    })
  } else {
    message(paste("⚠️ Skipping", metabolite, "- not enough data or zero variance"))
  }
}
dim(all_coefficients)
colnames(all_coefficients)
table(all_coefficients$term)
# Compute confidence intervals
all_coefficients <- all_coefficients %>% mutate( lower_CI = estimate - 1.96 * std.error,
    upper_CI = estimate + 1.96 * std.error)
# Apply Bonferroni correction to p-values
all_coefficients <- all_coefficients %>% mutate(adjusted_p_values = p.adjust(p.value, method = "fdr"))
# Save the results to CSV
write_csv(all_coefficients, "all_coefficients_results_Residualized_HILIC_lm_model.csv")
#########
###########################################################################################################################
# Initialize dataframe for Cohen's D
cohen_d_results <- data.frame()
# Loop through each metabolite and compute Cohen's D for lifestyle
for (metabolite in colnames(residualized_HILIC)[metabolite_columns]) {
  current_data <- residualized_HILIC %>%
    dplyr::select(Metabolite = all_of(metabolite), MW_scores_lifestyle_2) %>%
    filter(!is.na(Metabolite), !is.na(MW_scores_lifestyle_2))
  if (n_distinct(current_data$MW_scores_lifestyle_2) == 2) {
    # Split groups
    group_urban <- current_data %>% filter(MW_scores_lifestyle_2 == "Urban") %>% pull(Metabolite)
    group_nonurban <- current_data %>% filter(MW_scores_lifestyle_2 != "Urban") %>% pull(Metabolite)
    # Calculate pooled SD
    sd_pooled <- sqrt(((length(group_urban) - 1) * var(group_urban) + 
                         (length(group_nonurban) - 1) * var(group_nonurban)) /
                        (length(group_urban) + length(group_nonurban) - 2))
    # Cohen's D
    d <- (mean(group_urban) - mean(group_nonurban)) / sd_pooled
    cohen_d_results <- bind_rows(cohen_d_results, data.frame(Metabolite = metabolite, Cohens_d = d))
  }
}
head(cohen_d_results)
#########
###########################################################################################################################
## Plots
## coefficient plot
# Filter: Remove intercept & irrelevant terms - Filter for significant coefficient
filtered_coefficients <- all_coefficients %>%
  filter( !term %in% c("(Intercept)"), lower_CI > 0 | upper_CI < 0,
    abs(estimate) >= 0.35) %>% arrange(desc(abs(estimate)))  # Sort by effect size for clearer visualization
table(filtered_coefficients$term)
# Define custom color
color_mapping_b <- c("Age" = "#D67236","MW_scores_lifestyle_2Urban" = "#8B0000",
  "SexMale" = "cornflowerblue", "MW_scores_lifestyle_2Urban:SexMale" = "#8BAF9F")
# Create the coefficient plot
Co_plot <- ggplot(filtered_coefficients, aes(x = estimate, y = reorder(Metabolite, estimate), color = term)) +
  geom_point(size = 3) +  # Points for estimated effects
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2) +  # Confidence intervals
  scale_color_manual(values = color_mapping_b) +  # Custom colors
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Reference line at 0
  labs(title = "Significant Predictors of Metabolite Abundance",
    subtitle = "Linear model  β estimates adjusted for Age and Sex",
    x = "β  ≥ 0.5, 95% CI Excludes 0", y = "Metabolite",
    color = "Predictor") + theme_minimal() + theme(legend.position = "right")
# Print the Co_plot
print(Co_plot)
# Save the combined plot
ggsave("LME_coefficient_plot_RESIDUALIZED.png", plot = Co_plot, width = 10, height = 8, dpi = 300)

## Term Count Plot + Volcano Plot
# Filter to find significant results based on P-values
significant_results <- all_coefficients %>% filter(adjusted_p_values < 0.05) %>% arrange(adjusted_p_values)
# Remove unwanted terms like intercept, PC1
significant_results_no_intercept <- significant_results %>% filter(!term %in% c("(Intercept)"))
table(significant_results_no_intercept$Metabolite)
# Summarize data by predictors (MW_scores_h_sol_trim, Age, Sex, interaction terms)
term_counts <- significant_results_no_intercept %>%
  group_by(term) %>% summarise(Count = n_distinct(Metabolite), .groups = 'drop')
# Mapping of old term names to new labels
term_labels <- c( "Age" = "Age",       
  "MW_scores_lifestyle_2Urban" = "Lifestyle",
  "MW_scores_lifestyle_2Urban:SexMale" = "Lifestyle_X_Sex_Interaction",
  "SexMale" = "Sex",
  "Non-Significant" = "Non-Significant")

# Apply renaming to term_counts and pie_data
term_counts <- term_counts %>% mutate(term = recode(term, !!!term_labels))
# Total number of metabolites tested
total_metabolites <- n_distinct(all_coefficients$Metabolite)
# Compute number of non-significant metabolites
non_significant_count <- total_metabolites - sum(term_counts$Count)
# Create pie chart data including non-significant metabolites
pie_data <- term_counts %>%
  mutate(Percentage = (Count/total_metabolites) * 100) %>% add_row(term = "Non-Significant", Count = non_significant_count, Percentage = (non_significant_count / total_metabolites) * 100)
# Define custom colors for each predictor
custom_colors <- c("Age" = "#D67236", "Lifestyle" = "#8B0000",
  "Lifestyle_X_Sex_Interaction" = "#8BAF9F",    
  "Sex" = "cornflowerblue", "Non-Significant" = "azure2")
# Bar plot: Number of significant metabolites per predictor
bar_plot <- ggplot(term_counts, aes(x = reorder(term, Count), y = Count, fill = term)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  labs(x = "Predictor", y = "Count") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
print(bar_plot)
# Calculate the midpoints of each segment for accurate label positioning
pie_data <- pie_data %>% arrange(desc(term)) %>%
  mutate(ypos = cumsum(Percentage) - (Percentage / 2))
pie_chart <- ggplot(pie_data, aes(x = "", y = Percentage, fill = term)) +
  geom_bar(stat = "identity", width = 1) +  # White border to separate slices
  coord_polar(theta = "y", clip = "off") +  # Prevents clipping of labels
  geom_text_repel( aes(y = ypos, label = paste0(round(Percentage, 1), "%")),  
    nudge_x = 0.8,   # Move labels outward
    size = 5, direction = "y",  # Keeps text aligned to segments
    box.padding = 0.2, point.padding = 0.3,
    segment.size = 0.7, force = 1,
    show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +  # Apply same colors as bar plot
  labs(fill = "Predictor") + theme_void()
# Combine the bar plot and the pie chart
combined_plot <- ggdraw() + draw_plot(bar_plot, 0, 0, 0.75, 0.75) +  # Main bar plot
  draw_plot(pie_chart, 0.0005, 0.5, 0.5, 0.5)  # Inset pie chart (top left)
# Print the combined plot
print(combined_plot)
# Save the combined plot
ggsave("Significant_Metabolites_BarPlot_with_PieChart_RESIDUALIZED.png", plot = combined_plot, width = 10, height = 8, dpi = 300)
##########################################################################################################
# Identify top 50 metabolites by effect size for Lifestyle component
top_50_metabolites_h_sol <- significant_results %>%
  filter(term == "MW_scores_lifestyle_2Urban") %>% arrange(desc(abs(estimate))) %>% head(30)
# Create a bar plot for the top 50 metabolites by signed effect size with \(N\) on top
top_50_plot_signed <- ggplot(top_50_metabolites_h_sol, aes(x = reorder(Metabolite, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = estimate > 0)) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) + # Add \(N\) as text on top of bars
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  labs(title = "Top Metabolites:Ordered by lifestyle's effect size on the Metabolite",
    x = "Metabolite", y = "Beta value") + coord_flip() + theme_minimal()
# Print the top 50 plot
print(top_50_plot_signed)
# Save the plot
ggsave("Top_Metabolites_Discrete_lifestyle_RESIDUALIZED.png", plot = top_50_plot_signed, width = 8, height = 6, dpi = 300)
#
# Volcano Plot with the results from the regression combined with a fold change analysis. 
# Add log2 fold change and p-value significance
all_coefficients <- all_coefficients %>% mutate(log2_fold_change = estimate * log2(10))  # Convert from log10 scale

fc_threshold <- 1
volcano_data <- all_coefficients %>%
  filter(term == "MW_scores_lifestyle_2Urban") %>% mutate(log_p_value = -log10(adjusted_p_values),
    significant = adjusted_p_values < 0.00005 & abs(log2_fold_change) > fc_threshold,
    regulation = case_when( log2_fold_change > fc_threshold & adjusted_p_values < 0.00005 ~ "Abundant in Urban",
      log2_fold_change < -fc_threshold & adjusted_p_values < 0.00005 ~ "Abundant in Non Urban",
      TRUE ~ "Not Significant"))

# Merge Cohen's D values and define highlight group for top 35
volcano_data <- left_join(volcano_data, cohen_d_results, by = "Metabolite")

top_cohens <- volcano_data %>% filter(abs(Cohens_d) > 0.9) %>% pull(Metabolite)
volcano_data <- volcano_data %>% mutate(highlight_group = case_when(
      Metabolite %in% top_cohens & log2_fold_change > 0 ~ "Abundant in Urban",
      Metabolite %in% top_cohens & log2_fold_change < 0 ~ "Abundant in Non Urban",
      TRUE ~ "Not Top 35"))

highlight_colors <- c("Abundant in Urban" = "#FD6467",
  "Abundant in Non Urban" = "#798234", "Not Top 35" = "gray80")

# Volcano plot colored by Cohen’s D highlight group
volcano_plot_d <- ggplot(volcano_data, aes(x = log2_fold_change, y = log_p_value)) +
  geom_point(aes(color = highlight_group), alpha = 0.85, size = 3) +
  scale_color_manual(values = highlight_colors) +
  geom_hline(yintercept = -log10(0.00005), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(
    data = volcano_data %>% filter(Metabolite %in% top_cohens),
    aes(label = Metabolite, color = highlight_group),
    size = 3, box.padding = 0.3,
    point.padding = 0.2, max.overlaps = 40,
    force = 1, show.legend = FALSE ) + labs(
    title = "Volcano Plot: Top Metabolites by Cohen's D (Urban vs Non-Urban)",
    subtitle = "I colored those where Cohen's D (effect size) above 0.9",
    x = "Log2 Fold Change", y = "-log10(Adjusted P-Value)",
    color = "Legend") +theme_minimal()
print(volcano_plot_d)
ggsave("Volcano_CohensD_Top35.png", plot = volcano_plot_d, width = 8, height = 6, dpi = 300)
##########################################################################################################
# Volcano plot colored by Fold Change regulation ----and names i CARE ABOUT
# Label only pre-selected text targets
text_labels_for <- c("C23:2", "C20:1", "C25:1", "C25:3", "C26:2", "C26:1", "C19:0", "Ecgonine isomer", "Decatrienoyl carnitine C10:3", "decenoylcarnitine", "decadienoylcarnitine", "Decanoyl-L-carnitine", "3-Indolebutyric acid","MG(22:4)", "Trihydroxycoprostanoic acid", "Guanidinosuccinic acid","Ginsenoside Rh3 isomer", "DL cotinine", "Homostachydrine", "PE(18:0/18:0)", "LysoPE(20:0)", "LysoPE(22:1)","8-Amino-7-oxononanoic acid", "Pyrazinemethanethiol", "Indole-3-carboxilic acid-O-sulphate", "6-beta-hydrocortisol", "Phosphocreatinine", "Octanoylcarnitine", "Sphinganine", "dodecenoylcarnitine", "D-gluconate")
####
reg_colors <- c("Abundant in Urban" = "#FD6467", "Abundant in Non Urban" = "#798234","Not Significant" = "grey")
####
volcano_plot_fc <- ggplot(volcano_data, aes(x = log2_fold_change, y = log_p_value)) +
  geom_point(aes(color = regulation), alpha = 0.8, size = 3) +
  scale_color_manual(values = reg_colors) +
  geom_hline(yintercept = -log10(0.00005), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(data = volcano_data %>% filter(significant, Metabolite %in% text_labels_for),
    aes(label = Metabolite), size = 3, box.padding = 0.3, point.padding = 0.2, max.overlaps = 20) +
  labs(title = "Volcano Plot: Fold Change (Urban vs Non-Urban)",
    subtitle = "I colored those where FOLD CHANGE above 1",
    x = "Log2 Fold Change", y = "-log10(Adjusted P-Value)", color = "Regulation") + theme_minimal()
print(volcano_plot_fc)
volcano_plot_fc <- ggplot(volcano_data, aes(x = log2_fold_change, y = log_p_value)) +
  geom_point(aes(color = regulation), alpha = 0.8, size = 3) +
  scale_color_manual(values = reg_colors) +
  geom_hline(yintercept = -log10(0.00005), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(data = volcano_data %>% filter(significant), aes(label = Metabolite),
    size = 3, box.padding = 0.3, point.padding = 0.2, max.overlaps = 50) +
  labs(title = "Volcano Plot: Fold Change (Urban vs Non-Urban)",
    subtitle = "I colored those where FOLD CHANGE above 1",
    x = "Log2 Fold Change", y = "-log10(Adjusted P-Value)", color = "Regulation") + theme_minimal()

print(volcano_plot_fc)
ggsave("Volcano_FoldChange_Labeled.png", plot = volcano_plot_fc, width = 8, height = 6, dpi = 300)
#ggsave("Volcano_Plot_Regression_plus_FCRESIDUALIZED.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
##########################################################################################################
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
  labs(title = "Top Metabolites Affected by Age",
    x = "Metabolite", y = "Beta value") + coord_flip() + theme_minimal()
# Print the plot for Age
print(top_50_plot_age)
# Save the plot for Age
ggsave("Top_Metabolites_Affected_by_AgeRESIDUALIZED.png", plot = top_50_plot_age, width = 8, height = 6, dpi = 300)
################################################################################################################
# Identify top 50 metabolites by effect size for Sex
# Calculate log2 fold change - VOLCANO PLOT FOR SEX DIFFERENCES
all_coefficients <- all_coefficients %>%
  mutate(log2_fold_change = log2(exp(estimate)))  # Convert regression estimate to log2 fold change
# Define fold change threshold
fc_threshold <- 0.2  
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
    values = c("Higher in Males" = "#7294D4", "Higher in Females" = "indianred1", "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  
  ggrepel::geom_text_repel(data = volcano_data_sex %>% filter(significant),
    aes(label = Metabolite),size = 3,box.padding = 0.3,point.padding = 0.2,max.overlaps = 15) +
  labs(title = "Volcano Plot: Metabolite Differences Between Males and Females",
    x = "Log2 Fold Change (Male/Female)",y = "-log10(Adjusted P-Value)") + theme_minimal()
# Print and save the plot
print(volcano_plot_sex)
ggsave("Volcano_Plot_Sex_Corrected_RESIDUALIZED.png", plot = volcano_plot_sex, width = 8, height = 6, dpi = 300)
### Compute Cohen's D for Male vs Female
cohen_d_sex <- data.frame()

for (metabolite in colnames(residualized_HILIC)[2:610]) {
  current_data <- residualized_HILIC %>%
    dplyr::select(Metabolite = all_of(metabolite), Sex) %>%
    filter(!is.na(Metabolite), !is.na(Sex))
  
  if (n_distinct(current_data$Sex) == 2) {
    male_vals <- current_data %>% filter(Sex == "Male") %>% pull(Metabolite)
    female_vals <- current_data %>% filter(Sex == "Female") %>% pull(Metabolite)
    
    if (length(male_vals) > 2 && length(female_vals) > 2) {
      sd_pooled <- sqrt(((length(male_vals) - 1) * var(male_vals) + 
                           (length(female_vals) - 1) * var(female_vals)) /
                          (length(male_vals) + length(female_vals) - 2))
      
      d <- (mean(male_vals) - mean(female_vals)) / sd_pooled
      cohen_d_sex <- bind_rows(cohen_d_sex, data.frame(Metabolite = metabolite, Cohens_d = d))
    }
  }
}

### Step 2: Merge into volcano_data_sex
volcano_data_sex <- volcano_data_sex %>%
  left_join(cohen_d_sex, by = "Metabolite")

### Step 3: Plot Volcano Colored by Cohen's D
volcano_plot_sex_d <- ggplot(volcano_data_sex, aes(x = log2_fold_change, y = log_p_value)) +
  geom_point(aes(color = Cohens_d), alpha = 0.85, size = 3) +
  scale_color_gradient2(low = "blue", mid = "gray80", high = "red", midpoint = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(
    data = volcano_data_sex %>% filter(significant & abs(Cohens_d) > 0.9),
    aes(label = Metabolite), size = 3, box.padding = 0.3, point.padding = 0.2,
    max.overlaps = 25, show.legend = FALSE
  ) +
  labs(
    title = "Volcano Plot: Sex Differences in Metabolite Abundance",
    subtitle = "Color gradient shows Cohen's D effect size (Male vs Female)",
    x = "Log2 Fold Change (Male/Female)",
    y = "-log10(Adjusted P-Value)",
    color = "Cohen's D"
  ) +
  theme_minimal()

print(volcano_plot_sex_d)
ggsave("Volcano_Sex_CohensD.png", plot = volcano_plot_sex_d, width = 8, height = 6, dpi = 300)

##########################################################################################################
# Initialize storage
interaction_coeffs <- data.frame()
interaction_r2 <- data.frame()
interaction_model_data <- list()
# Metabolites to loop over
metabolite_columns <- 2:610
for (metabolite in colnames(residualized_HILIC)[metabolite_columns]) {
  message(paste0("🔍 Processing: ", metabolite))
  current_data <- residualized_HILIC %>%
    dplyr::select(indiv.ID, MW_scores_lifestyle_2, Age, Sex, batch, Sampling.location, all_of(metabolite)) %>%
    rename(Metabolite = all_of(metabolite)) %>%
    filter(!is.na(Metabolite), !is.na(Age), !is.na(Sex))
  # Outlier removal
  if (nrow(current_data) > 0) {
    Q1 <- quantile(current_data$Metabolite, 0.25, na.rm = TRUE)
    Q3 <- quantile(current_data$Metabolite, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    current_data <- current_data %>%
      filter(Metabolite >= (Q1 - 1.5 * IQR), Metabolite <= (Q3 + 1.5 * IQR)) }
  ######
  if (nrow(current_data) > 6 && var(current_data$Metabolite, na.rm = TRUE) > 0) {
    tryCatch({
      # Fit base and interaction models
      model_base <- lm(Metabolite ~ MW_scores_lifestyle_2 + Sex + Age, data = current_data)
      model_inter <- lm(Metabolite ~ MW_scores_lifestyle_2 * Sex + Age, data = current_data)
      ###########
      r2_base <- summary(model_base)$adj.r.squared
      r2_inter <- summary(model_inter)$adj.r.squared
      #########
      # Extract tidy results
      temp_results <- broom::tidy(model_inter) %>%
        mutate(Metabolite = metabolite, N = nrow(current_data))
      #########
      # Find the interaction term dynamically
      interaction_row <- temp_results %>%
        filter(grepl("MW_scores_lifestyle_2.*:Sex.*", term))
      ########
      # Check if R2 improved and interaction term is present
      if (nrow(interaction_row) == 1 && !is.na(r2_inter) && r2_inter > r2_base) {
        interaction_coeffs <- bind_rows(interaction_coeffs, temp_results)
        interaction_r2 <- bind_rows(interaction_r2, data.frame(
          r2_base = r2_base, r2_inter = r2_inter, Metabolite = metabolite
        ))
        interaction_model_data[[metabolite]] <- list(data = current_data,model = model_inter)
        message(paste0("✅ Interaction improves fit: ", metabolite))
      } else {
        message(paste0("ℹ️ Skipping ", metabolite, ": no improvement or missing interaction term"))
      }
    }, error = function(e) {
      message(paste("❌ Error in", metabolite, ":", e$message))
    })
  } else {
    message(paste("⚠️ Skipping", metabolite, "- insufficient data or variance"))
  }
}
# 1. Count of metabolites with improved R²
cat("📈 R² improved for:", nrow(interaction_r2), "metabolites\n")
# 2. Count of metabolites where the interaction term is FDR-significant
interaction_signif <- interaction_coeffs %>% filter(grepl("MW_scores_lifestyle_2.*:Sex.*", term)) %>%
  mutate(fdr_p = p.adjust(p.value, method = "fdr")) %>% filter(fdr_p < 0.05)
cat("⭐ FDR-significant interaction term (p.adj < 0.05) in:", nrow(interaction_signif), "metabolites\n")
# Add CI + FDR-corrected p-values
interaction_coeffs <- interaction_coeffs %>%
  mutate(
    lower_CI = estimate - 1.96 * std.error,
    upper_CI = estimate + 1.96 * std.error
  )

# FDR correction only for interaction term
interaction_term_signif <- interaction_coeffs %>%
  filter(term == "MW_scores_lifestyle_2Urban:SexMale") %>%
  mutate(fdr_p = p.adjust(p.value, method = "fdr")) %>% filter(fdr_p < 0.05)
head(interaction_term_signif)
# Keep only significant ones where model fit improved
significant_metabolites <- interaction_term_signif$Metabolite
######
# Define your color scheme
Sex_colors <- c("Female" = "#E6A0C4", "Male" = "#7294D4")
# Save the plots list
interaction_plots <- list()
# Set fixed width/height in inches
plot_width <- 4
plot_height <- 4
for (metabolite in significant_metabolites) { mod_data <- interaction_model_data[[metabolite]]
  if (!is.null(mod_data$data)) { df <- mod_data$data
    summary_df <- df %>% group_by(MW_scores_lifestyle_2, Sex) %>%
      summarise( mean = mean(Metabolite, na.rm = TRUE),
        sd = sd(Metabolite, na.rm = TRUE),
        n = n(),se = sd / sqrt(n),
        ci_lower = mean - 1.96 * se,
        ci_upper = mean + 1.96 * se,
        .groups = "drop")
    # Build plot
    p <- ggplot(summary_df, aes(x = MW_scores_lifestyle_2, y = mean, group = Sex, color = Sex)) +
      geom_line(size = 1.2, position = position_dodge(0.3)) +
      geom_point(size = 2.5, position = position_dodge(0.3)) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, position = position_dodge(0.3)) +
      scale_color_manual(values = Sex_colors) +
      labs(title = metabolite,
        y = "Estimated Mean ± 95% CI", x = "Lifestyle") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text = element_text(size = 14),
            aspect.ratio = 1)
    # Store the plot
    interaction_plots[[metabolite]] <- p
    # Print to viewer
    print(p)
    message(paste("✅ Plot created and shown for:", metabolite))
  } else {
    message(paste("⚠️ No data found for:", metabolite))
  }
}
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
grp_cols <- c(Urban = "#FD6467", Non_Urban = "#798234","Pastoralist"= "#29AB87")

## ───────────────────────────── 1. Residualize metabolite data ──────────────────────────────
residuals_long <- WIDE_dataSub_with_metadata %>%
  pivot_longer(cols = 2:629, names_to = "Metabolite", values_to = "Abundance") %>%
  drop_na(Abundance, batch, Sex, Age) %>%
  group_by(Metabolite) %>%
  filter(all(is.finite(Abundance)) && var(Abundance) > 0) %>%
  mutate(Residual = add_residuals(cur_data(), lm(Abundance ~ batch + Sex + Age, data = cur_data()))$resid) %>%
  ungroup() %>%
  select(indiv.ID, Metabolite, Residual)

residuals_wide <- residuals_long %>%
  pivot_wider(names_from = Metabolite, values_from = Residual)

metadata <- WIDE_dataSub_with_metadata %>% select(indiv.ID, MW_scores_lifestyle_2, MW_scores_lifestyle)
hilic_data <- left_join(residuals_wide, metadata, by = "indiv.ID")
dim(hilic_data)
head(hilic_data)
colnames(WIDE_dataSub_with_metadata)
metadata_2 <- WIDE_dataSub_with_metadata %>% select(indiv.ID, MW_scores_lifestyle_2,Unique.ID)
hilic_data_2 <- left_join(residuals_wide, metadata, by = "indiv.ID")
write.csv(hilic_data_2, "HILIC_Residualized_Only_Lifestyle_effect_left.csv", row.names = FALSE)

table(hilic_data$MW_scores_lifestyle)
hilic_pas_urb <- hilic_data %>% filter(MW_scores_lifestyle != "PeriUrban") %>% 
  mutate( MW_scores_lifestyle = factor(
      MW_scores_lifestyle, levels = c("Pastoralist", "Urban")))
# sanity check
table(hilic_pas_urb$MW_scores_lifestyle)
# ---------------------------
# 1. Fit PLS-DA Model on Real Labels
# ---------------------------
set.seed(123)
X <- hilic_data %>% select(-indiv.ID, -MW_scores_lifestyle_2) %>% as.matrix() # <— drop PROBLEMATIC METABOLITES here
Y <- factor(hilic_data$MW_scores_lifestyle_2)

train_idx <- createDataPartition(Y, p = 2/3, list = FALSE)
train_df <- data.frame(X[train_idx, ]) %>% mutate(Group = Y[train_idx])
test_df  <- data.frame(X[-train_idx, ])
Y_test   <- Y[-train_idx]

ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3)
pls_fit <- train(Group ~ ., data = train_df,
                 method = "pls", preProcess = c("center", "scale"),
                 trControl = ctrl, tuneLength = 15, metric = "Accuracy")
best_ncomp <- pls_fit$bestTune$ncomp
pls_model <- pls_fit$finalModel

# Compute observed VIP scores
compute_vip <- function(mvr_obj, ncomp) {
  W <- mvr_obj$loading.weights[, 1:ncomp, drop = FALSE]
  T2 <- colSums(mvr_obj$scores[, 1:ncomp, drop = FALSE]^2)
  vip <- sqrt(ncol(W) * rowSums((W^2) %*% diag(T2)) / sum(T2))
  names(vip) <- rownames(W)
  return(vip)
}
obs_vip <- compute_vip(pls_model, best_ncomp)

# ---------------------------
# 2. Permutation Testing for VIP Significance
# ---------------------------
n_perm <- 1000
vip_null <- matrix(NA, nrow = length(obs_vip), ncol = n_perm)
rownames(vip_null) <- names(obs_vip)

for (i in 1:n_perm) {
  perm_group <- sample(train_df$Group)
  perm_fit <- plsda(X[train_idx, ], perm_group, ncomp = best_ncomp, scale = TRUE)
  vip_null[, i] <- compute_vip(perm_fit, best_ncomp)
  if (i %% 50 == 0) cat("Permutation:", i, "\n")
}

# ---------------------------
# 3. Compute Empirical p-values & Adjust
# ---------------------------
p_vals <- sapply(1:length(obs_vip), function(i) {
  mean(vip_null[i, ] >= obs_vip[i])
})
names(p_vals) <- names(obs_vip)

# Adjust for multiple testing
p_adj <- p.adjust(p_vals, method = "fdr")  # or use "BH"

# ---------------------------
# 4. Plot Top with Significant Asterisks
# ---------------------------
top_metabs <- names(sort(obs_vip, decreasing = TRUE))[1:36]
top_loads <- pls_model$loadings[top_metabs, 1:3]
plot_df <- as.data.frame(top_loads) %>%
  rownames_to_column("Metabolite") %>%
  pivot_longer(-Metabolite, names_to = "Component", values_to = "Loading") %>%
  mutate(pval = p_vals[Metabolite],
         padj = p_adj[Metabolite],
         sig = ifelse(pval < 0.05, "*", ""))
head(plot_df)
# ---------------------------
# 5. Heatmap Plot
# ---------------------------
ggplot(plot_df, aes(x = Component, y = fct_reorder(Metabolite, Loading), fill = Loading)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5, vjust = 0.5) +
  scale_fill_gradient2(low = "red", high = "green", mid = "white", midpoint = 0) +
  labs(title = "TopMetabolite Loadings on the first 3 PLS-DA Components",
    subtitle = "Asterix shows p-value significance,\ni.e. testing if VIP score is unusually high",
    x = "PLS Component", y = NULL, fill = "Loading") + theme_minimal(base_size = 14)


# Compute the mean of the permuted VIPs for each metabolite
vip_mean_null <- rowMeans(vip_null, na.rm = TRUE)
# Combine all relevant values into one dataframe
vip_summary_df <- tibble(Metabolite = names(obs_vip),
  VIP_observed   = obs_vip,
  VIP_null_mean  = vip_mean_null,
  P_value        = p_vals,
  P_adj_FDR      = p_adj)
# Extract Comp1 loadings for top metabolites
top_vip_df <- tibble(Metabolite = top_metabs,
  VIP = obs_vip[top_metabs],
  Loading_Comp1 = pls_model$loadings[top_metabs, 1],
  P_value = p_vals[top_metabs])

# Classify based on sign of Component 1 loading
top_vip_df <- top_vip_df %>%
  mutate(Group_Abundance = ifelse(Loading_Comp1 < 0, "Urban", "Non_Urban"))


# Plot
ggplot(top_vip_df, aes(x = VIP, y = fct_reorder(Metabolite, VIP), color = Group_Abundance)) +
  geom_point(size = 4) +
  scale_color_manual(values = grp_cols, name = "Which Group Has More") +
  labs(title = "Top Metabolites by VIP Score",
    subtitle = "Dot color indicates whether metabolite is more\nabundant in Urban or Non-Urban\nAlso, I plot only those with significant p-value",
    x = "VIP Score", y = NULL) + theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 9))

# Choose a metabolite to visualize (e.g., one with significant VIP)
# Filter top 36 metabolites with raw p < 0.05
top36 <- vip_summary_df %>% filter(P_value < 0.05) %>%
  arrange(P_value) %>% slice_head(n = 36) %>% pull(Metabolite)

# Create tidy data for plotting
vip_plot_df <- map_dfr(top36, function(metab) {
  tibble(Metabolite = metab, Permuted_VIP = vip_null[metab, ], Observed_VIP = obs_vip[metab])})

# Plot grid
ggplot(vip_plot_df, aes(x = Permuted_VIP)) +
  geom_histogram(binwidth = 0.0005, fill = "skyblue", color = "skyblue") +
  geom_vline(aes(xintercept = Observed_VIP), color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Metabolite, scales = "free", ncol = 7) +
  labs(title = "VIP Permutation Distributions for Metabolites with p < 0.05",
       subtitle = "Red line = observed VIP; blue = null distribution",
       x = "Permuted VIP Score", y = "Frequency" ) + theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 8), axis.text = element_text(size = 7))
##########################################################################################################
# ---------------------------
# Accuracy on Test Data
# ---------------------------
test_pred <- predict(pls_fit, newdata = test_df)
true_acc <- mean(test_pred == Y_test)

# Print accuracy in percentage
cat(sprintf("✅ Test Set Accuracy: %.2f%%\n", true_acc * 100))

# ---------------------------
# Null Distribution via Permutation
# ---------------------------
perm_acc <- numeric(n_perm)
for (i in 1:n_perm) {
  perm_group <- sample(Y)
  train_perm <- data.frame(X[train_idx, ]) %>% mutate(Group = perm_group[train_idx])
  perm_model <- train(Group ~ ., data = train_perm,
                      method = "pls", preProcess = c("center", "scale"),
                      trControl = ctrl, tuneLength = 15)
  pred <- predict(perm_model, newdata = test_df)
  perm_acc[i] <- mean(pred == Y_test)
}

# ---------------------------
# Plot Accuracy Distribution
# ---------------------------
ggplot(tibble(Permuted_Accuracy = perm_acc), aes(x = Permuted_Accuracy)) +
  geom_histogram(binwidth = 0.02, fill = "grey80", color = "black") +
  
  # Observed accuracy line
  geom_vline(xintercept = true_acc, color = "red", linetype = "dashed", linewidth = 1.2) +
  annotate("text",
           x = true_acc,
           y = max(table(cut(perm_acc, breaks = 30))),
           label = sprintf("Observed Accuracy: %.2f%%", true_acc * 100),
           vjust = -0.5, hjust = 1.05, size = 5, color = "red") +
  
  # Mean null accuracy line
  geom_vline(xintercept = mean(perm_acc), color = "blue", linetype = "dotted", linewidth = 1) +
  annotate("text",
           x = mean(perm_acc),
           y = max(table(cut(perm_acc, breaks = 30))),
           label = sprintf("Mean Null: %.2f%%", mean(perm_acc) * 100),
           vjust = 1.5, hjust = -0.1, size = 4.5, color = "blue") +
  
  labs(
    title = "PLS-DA Accuracy Compared to Permuted Labels",
    subtitle = "Dashed red = observed model accuracy on test set\nDotted blue = mean of null distribution",
    x = "Accuracy on Test Set",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 16)

##########################################################################################################
# Extract performance from PLS model fit
results_df <- pls_fit$results

# Reshape for plotting
plot_df <- results_df %>%
  select(ncomp, Accuracy, Kappa) %>%
  pivot_longer(cols = c(Accuracy, Kappa), names_to = "Metric", values_to = "Value")

# Plot: Accuracy vs Components
ggplot(plot_df, aes(x = ncomp, y = Value, color = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(Accuracy = "#1f77b4", Kappa = "#ff7f0e")) +
  labs(
    title = "PLS-DA Model Performance by Number of Components",
    x = "Number of PLS Components",
    y = "Cross-Validated Metric",
    color = "Metric"
  ) +
  theme_minimal(base_size = 14)

##########################################################################################################
# Extract scores from the final PLS model
# Extract scores from the final PLS model
pls_scores <- pls_model$scores[, 1:2]
score_df <- data.frame(
  LV1 = pls_scores[, 1],
  LV2 = pls_scores[, 2],
  Group = Y[train_idx]
)

# Calculate variance explained
explained_var <- pls_model$Xvar / sum(pls_model$Xvar) * 100
percent_lv1 <- round(explained_var[1], 1)
percent_lv2 <- round(explained_var[2], 1)

# Plot
ggplot(score_df, aes(x = LV1, y = LV2, color = Group)) +
  geom_point(alpha = 0.8, size = 2.5) +
  stat_ellipse(level = 0.95, size = 1, linetype = "solid") +
  scale_color_manual(values = c("Non_Urban" = "#556B2F", "Urban" = "#E34234")) +
  labs(
    title = "PLS-DA Latent Variables Plot \nModel Fit with 1000 Permutations for VIP p-values",
    x = paste0("LV1 (", percent_lv1, "%)"),
    y = paste0("LV2 (", percent_lv2, "%)"),
    color = "Group"
  ) +
  theme_minimal(base_size = 15)

##########################################################################################################
##########################################################################################################
# Load required package
library(stringdist)
# STEP 1: Extract significant metabolite names (from VIP results)
sig_metabs <- vip_summary_df %>%
  filter(P_value < 0.05) %>%
  pull(Metabolite)

# STEP 2: Extract all metabolite column names from hilic_data (excluding metadata)
all_metab_names <- colnames(hilic_data)[!colnames(hilic_data) %in% c("indiv.ID", "MW_scores_lifestyle_2")]

# STEP 3: Use approximate string matching to align metabolite names
matched_indices <- amatch(sig_metabs, all_metab_names, maxDist = 3)
matched_metabs <- all_metab_names[matched_indices]

# STEP 4: Handle unmatched manually
unmatched <- sig_metabs[is.na(matched_indices)]
if (length(unmatched) > 0) {
  cat("Unmatched metabolites (couldn't be mapped):", paste(unmatched, collapse = ", "), "\n")
}

# Clean fuzzy matches
matched_metabs <- matched_metabs[!is.na(matched_metabs)]

# STEP 4B: Manually correct stubborn ones
manual_map <- c(
  "PC.15.0.22.5." = "PC(15:0/22:5)",
  "PC.16.0.19.1." = "PC(16:0/19:1)",
  "PC.17.0.18.2." = "PC(17:0/18:2)",
  "PC.17.0.20.4." = "PC(17:0/20:4)",
  "PC.21.3.16.0." = "PC(21:3/16:0)",
  "PE.18.0.18.0." = "PE(18:0/18:0)",
  "X4.Hydroxyisoleucine.isomer.1" = "4-Hydroxyisoleucine isomer 1",
  "X8.Amino.7.oxononanoic.acid" = "8-Amino-7-oxononanoic acid",
  "C26.0.Cerotic.acid." = "C26:0(Cerotic acid)",
  "n.ethyl.d.serine.isomer" = "n-ethyl-d-serine isomer"
)

# Filter manual mappings that exist in hilic_data
manual_matches <- manual_map[names(manual_map) %in% unmatched]
manual_matches <- manual_matches[manual_matches %in% colnames(hilic_data)]

# STEP 5: Combine both sets
final_metabolites <- unique(c(matched_metabs, manual_matches))

# STEP 6: Subset hilic_data to include only matched metabolites and lifestyle
pca_data <- hilic_data %>%
  select(all_of(final_metabolites), MW_scores_lifestyle_2)

# STEP 7: Perform PCA
pca_result <- prcomp(pca_data %>% select(-MW_scores_lifestyle_2), scale. = TRUE)

# STEP 8: Prepare scores for plotting
pca_scores <- as.data.frame(pca_result$x) %>%
  mutate(Lifestyle = pca_data$MW_scores_lifestyle_2)

# STEP 9: Plot PCA
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Lifestyle)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = c(Urban = "#FD6467", Non_Urban = "#798234")) +
  labs(
    title = "PCA of Metabolites with Significant VIP Scores",
    subtitle = "Metabolites selected from PLS-DA (VIP p < 0.05, permutation-based)",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)"),
    color = "Lifestyle Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")


cat("Matched:", length(final_metabolites), "out of", length(sig_metabs), "significant metabolites\n")
# Confirm all sig_metabs are accounted for
remaining_unmatched <- setdiff(sig_metabs, names(manual_map)[names(manual_map) %in% unmatched] %>% 
                                 union(sig_metabs[!is.na(matched_indices)]))

if (length(remaining_unmatched) == 0) {
  message("✅ All significant metabolites have been matched successfully!")
} else {
  warning("❌ The following significant metabolites were NOT matched:\n",
          paste(remaining_unmatched, collapse = ", "))
}
######################################################################################################################################
# Load plotly
library(plotly)


# Extract top 3 latent variable scores
pls_scores <- pls_model$scores[, 1:3]
score_df_3D <- data.frame(
  LV1 = pls_scores[, 1],
  LV2 = pls_scores[, 2],
  LV3 = pls_scores[, 3],
  Group = Y[train_idx]
)

# Variance explained
explained_var <- pls_model$Xvar / sum(pls_model$Xvar) * 100
percent_lv <- round(explained_var[1:3], 1)

# 3D Scatter Plot
plot_ly(
  score_df_3D,
  x = ~LV1, y = ~LV2, z = ~LV3,
  color = ~Group,
  colors = c("Non_Urban" = "#556B2F", "Urban" = "#E34234"),
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 4, opacity = 0.8)
) %>%
  layout(
    title = "3D PLS-DA Score Plot",
    scene = list(
      xaxis = list(title = paste0("LV1 (", percent_lv[1], "%)")),
      yaxis = list(title = paste0("LV2 (", percent_lv[2], "%)")),
      zaxis = list(title = paste0("LV3 (", percent_lv[3], "%)"))
    )
  )

#######################################################################################################################################
# Load plotly
library(plotly)

# Extract PC scores
pca_scores_3D <- as.data.frame(pca_result$x[, 1:3]) %>%
  mutate(Lifestyle = pca_data$MW_scores_lifestyle_2)

# Calculate % variance explained
expl_var <- summary(pca_result)$importance[2, 1:3]
expl_var <- round(expl_var * 100, 1)

# Plot
plot_ly(
  pca_scores_3D,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~Lifestyle,
  colors = c("Urban" = "#FD6467", "Non_Urban" = "#798234"),
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 4, opacity = 0.85)
) %>%
  layout(
    title = "3D PCA of VIP-Selected Metabolites",
    scene = list(
      xaxis = list(title = paste0("PC1 (", expl_var[1], "%)")),
      yaxis = list(title = paste0("PC2 (", expl_var[2], "%)")),
      zaxis = list(title = paste0("PC3 (", expl_var[3], "%)"))
    ),
    legend = list(orientation = "h", x = 0.3, y = -0.1)
  )

######################################################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
residualized_HILIC$MW_scores_lifestyle <- as.factor(residualized_HILIC$MW_scores_lifestyle)
table(residualized_HILIC$MW_scores_lifestyle_2)
table(residualized_HILIC$MW_scores_lifestyle)
colnames(residualized_HILIC)

# --------------------------
# 1. Run regression per metabolite and store coefficients
# --------------------------
all_coefficients <- data.frame()
all_r_squared <- data.frame()
metabolite_columns <- 2:610  # your metabolite columns

for (metabolite in colnames(residualized_HILIC)[metabolite_columns]) {
  message(paste0("🔍 Processing: ", metabolite))
  
  current_data <- residualized_HILIC %>%
    dplyr::select(indiv.ID, MW_scores_lifestyle, Age, Sex, batch, 
                  Sampling.location, all_of(metabolite)) %>%
    rename(Metabolite = all_of(metabolite)) %>%
    filter(!is.na(Metabolite), !is.na(Age), !is.na(Sex))
  
  # Outlier removal using IQR
  if (nrow(current_data) > 0) {
    Q1 <- quantile(current_data$Metabolite, 0.25, na.rm = TRUE)
    Q3 <- quantile(current_data$Metabolite, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    current_data <- current_data %>%
      filter(Metabolite >= (Q1 - 1.5 * IQR), Metabolite <= (Q3 + 1.5 * IQR))
  }
  
  if (nrow(current_data) > 6 && is.numeric(current_data$Metabolite) && var(current_data$Metabolite, na.rm = TRUE) > 0) {
    tryCatch({
      model <- lm(Metabolite ~ MW_scores_lifestyle + Sex + Age, data = current_data)
      temp_results <- broom::tidy(model) %>%
        mutate(Metabolite = metabolite, N = nrow(current_data))
      
      model_r_squared <- data.frame(
        r.squared = summary(model)$r.squared,
        r.squared_adj = summary(model)$adj.r.squared,
        Metabolite = metabolite
      )
      
      all_coefficients <- bind_rows(all_coefficients, temp_results)
      all_r_squared <- bind_rows(all_r_squared, model_r_squared)
      
      message(paste0("✅ Model succeeded for: ", metabolite))
    }, error = function(e) {
      message(paste("❌ Skipping", metabolite, "- model failed with error:", e$message))
    })
  } else {
    message(paste("⚠️ Skipping", metabolite, "- not enough data or zero variance"))
  }
}

# Confidence intervals
all_coefficients <- all_coefficients %>%
  mutate(lower_CI = estimate - 1.96 * std.error,
         upper_CI = estimate + 1.96 * std.error)

# Extract p-values for MW_scores_lifestyle from lm() output
# These are the two dummy-coded comparisons: PeriUrban vs Pastoralist and Urban vs Pastoralist
pvals <- all_coefficients %>%
  filter(term %in% c("MW_scores_lifestylePeriUrban", "MW_scores_lifestyleUrban")) %>%
  group_by(Metabolite) %>%
  summarise(min_p = min(p.value, na.rm = TRUE), .groups = "drop") %>%
  mutate(FDR_p = p.adjust(min_p, method = "bonferroni")) %>%
  filter(FDR_p < 0.05) %>%
  rename(p_lifestyle = min_p)
head(pvals)
# --------------------------
# 3. Compute pairwise Cohen's d PLUS the p-value using t.test
# --------------------------
cohen_d_results <- data.frame()

pairwise <- combn(c("Pastoralist", "PeriUrban", "Urban"), 2, simplify = FALSE)

for (met in colnames(residualized_HILIC)[metabolite_columns]) {
  df <- residualized_HILIC %>%
    dplyr::select(Metabolite = all_of(met), Group = MW_scores_lifestyle) %>%
    drop_na()
  
  for (pair in pairwise) {
    g1 <- df %>% filter(Group == pair[1]) %>% pull(Metabolite)
    g2 <- df %>% filter(Group == pair[2]) %>% pull(Metabolite)
    
    n1 <- length(g1)
    n2 <- length(g2)
    
    if (n1 > 2 && n2 > 2 && var(g1) > 0 && var(g2) > 0) {
      pooled_sd <- sqrt(((n1 - 1) * var(g1) + (n2 - 1) * var(g2)) / (n1 + n2 - 2))
      d <- (mean(g1) - mean(g2)) / pooled_sd
      
      # t-test
      p_val <- t.test(g1, g2, var.equal = TRUE)$p.value
      
      # 95% CI for Cohen’s d using standard error approximation
      se_d <- sqrt((n1 + n2) / (n1 * n2) + (d^2 / (2 * (n1 + n2))))
      ci_low <- d - 1.96 * se_d
      ci_high <- d + 1.96 * se_d
      
      cohen_d_results <- bind_rows(cohen_d_results,
                                   data.frame(
                                     Metabolite = met,
                                     Group1 = pair[1],
                                     Group2 = pair[2],
                                     Cohens_d = d,
                                     d_CI_lower = ci_low,
                                     d_CI_upper = ci_high,
                                     p_value = p_val
                                   ))
    }
  }
}
# Apply Bonferroni correction (you can switch to "fdr" if preferred)
cohen_d_results <- cohen_d_results %>% mutate(p_adj = p.adjust(p_value, method = "bonferroni"))
# Count where d > 1 and lower CI > 1
n_strong_effects <- cohen_d_results %>% filter(Cohens_d > 1, d_CI_lower > 1) %>% nrow()
cat("Number of pairwise comparisons where d > 1 and lower CI > 1:", n_strong_effects, "\n")
# Filter strong effect comparisons
strong_effects <- cohen_d_results %>% filter(Cohens_d > 1, d_CI_lower > 0.5)
# Just get unique metabolites for plotting (or filter for top pairs if needed)
metabolites_to_plot <- unique(strong_effects$Metabolite)
# Prepare long-format data for boxplot
plot_data <- residualized_HILIC %>%
  pivot_longer(cols = all_of(colnames(residualized_HILIC)[metabolite_columns]),
               names_to = "Metabolite", values_to = "Value") %>%
  rename(Group = MW_scores_lifestyle) %>%
  filter(Metabolite %in% metabolites_to_plot)
# Y-position helper: get 95th percentile per metabolite to place d-labels
y_max <- plot_data %>%
  group_by(Metabolite) %>%
  summarise(y_base = quantile(Value, 0.95, na.rm = TRUE)) %>%
  ungroup()

# Bar coordinates for annotations
annotate_df <- strong_effects %>%
  mutate(
    pair = paste(Group1, "vs", Group2),
    label = paste0("d = ", round(Cohens_d, 2)),
    x_start = case_when(Group1 == "Pastoralist" & Group2 == "PeriUrban" ~ 1,
                        Group1 == "Pastoralist" & Group2 == "Urban"      ~ 1,
                        Group1 == "PeriUrban"   & Group2 == "Urban"     ~ 2),
    x_end = case_when(Group1 == "Pastoralist" & Group2 == "PeriUrban" ~ 2,
                      Group1 == "Pastoralist" & Group2 == "Urban"     ~ 3,
                      Group1 == "PeriUrban"   & Group2 == "Urban"     ~ 3)
  ) %>%
  left_join(y_max, by = "Metabolite") %>%
  group_by(Metabolite) %>%
  mutate(y = y_base + 0.05 * row_number()) %>%
  ungroup()
# Get unique metabolites to plot (from annotate_df or strong_effects)
metabolites_to_plot <- unique(annotate_df$Metabolite)
# Define grid layout (you can change this as needed)
ncol <- 3
nrow <- 3
plots_per_page <- ncol * nrow
# Automatically calculate how many pages are needed
n_pages <- ceiling(length(metabolites_to_plot) / plots_per_page)
# Plot each page
for (i in 1:n_pages) {
  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
    
    # Annotations
    geom_segment(data = annotate_df,
                 aes(x = x_start, xend = x_end, y = y, yend = y),
                 inherit.aes = FALSE) +
    geom_text(data = annotate_df,
              aes(x = (x_start + x_end)/2, y = y + 0.02, label = label),
              size = 3, vjust = 0, inherit.aes = FALSE) +
    
    ggforce::facet_wrap_paginate(~ Metabolite,
                                 ncol, nrow,
                                 page = i,
                                 scales = "free") +
    scale_fill_manual(values = c("Pastoralist" = "#29AB87",
                                 "PeriUrban" = "orange",
                                 "Urban" = "#FD6467")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(size = 10)
    ) +
    labs(x = NULL, y = "Residual abundance")
  
  print(p)
  
  # Save each page (optional)
  ggsave(sprintf("CohensD_boxplot_page_%02d.png", i), p,
         width = 9, height = 9, dpi = 300)
}




#############
#########################################################################################################################################
############
########
###Making violin plots of those of interest 
# Define metabolites to select
selected_metabolites <- c("Trihydroxycoprostanoic acid", "Guanidinosuccinic acid")
selected_metabolites <- c("Ginsenoside Rh3 isomer", "DL cotinine", "Homostachydrine", "Perfluorooctanesulfonic acid")
selected_metabolites <- c("C20:0(Arachidic acid)", "C20:1", "C25:1", "C26:1", "a-ketoglutarate", "ceramide (d18:1/22:1)", "ceramide (d18:1/22:1)", "trihydroxycoprostanoic acid")
selected_metabolites <- c("C20:0(Arachidic acid)", "C20:1", "C25:1", "C26:1", "a-ketoglutarate", "ceramide (d18:1/22:1)", "ceramide (d18:1/22:1)", "trihydroxycoprostanoic acid", "Ginsenoside Rh3 isomer", "DL cotinine", "Homostachydrine", "Perfluorooctanesulfonic acid", "Trihydroxycoprostanoic acid", "Guanidinosuccinic acid")
colnames(WIDE_dataSub_with_metadata$`Perfluorooctanesulfonic acid`)
# Filter metabolomics dataset for selected metabolites and reshape to long format
Selected_metab <- WIDE_dataSub_with_metadata %>%
  dplyr::select(indiv.ID, MW_scores_lifestyle, all_of(selected_metabolites)) %>%
  pivot_longer(cols = all_of(selected_metabolites), names_to = "Metabolite", values_to = "Concentration") %>%
  filter(!is.na(Concentration))  # Remove missing values
# Function to remove outliers using IQR method
remove_outliers <- function(df) {
  df %>%
    group_by(Metabolite) %>%
    mutate(Q1 = quantile(Concentration, 0.25, na.rm = TRUE),
           Q3 = quantile(Concentration, 0.75, na.rm = TRUE),
           IQR = Q3 - Q1,
           Lower_Bound = Q1 - 1.5 * IQR,
           Upper_Bound = Q3 + 1.5 * IQR) %>%
    filter(Concentration >= Lower_Bound & Concentration <= Upper_Bound) %>%
    dplyr::select(-Q1, -Q3, -IQR, -Lower_Bound, -Upper_Bound)  # Remove extra columns
}

# Apply outlier removal
Selected_metab_clean <- remove_outliers(Selected_metab)
head(Selected_metab_clean)
# ---- Define Significance Threshold ----
significance_threshold <- 0.05  

# ---- Filter for Significant Metabolites ----
significant_metabolites <- all_coefficients %>%
  filter(term == "MW_scores_lifestyleUrban" & Metabolite %in% selected_metabolites & p.value < significance_threshold) %>%
  mutate(
    MW_scores_lifestyle = "Urban",  
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  left_join(
    Selected_metab_clean %>%
      group_by(Metabolite, MW_scores_lifestyle) %>%
      summarize(
        y_position = max(Concentration, na.rm = TRUE) + 0.15 * diff(range(Concentration, na.rm = TRUE)),
        .groups = "drop"
      ),
    by = c("Metabolite", "MW_scores_lifestyle")
  )

# ---- Keep Only Significant Metabolites in the Data ----
#Selected_metab_clean_filtered <- Selected_metab_clean %>%
  #filter(Metabolite %in% significant_metabolites$Metabolite)
## Plot all of them
Selected_metab_clean_filtered <- Selected_metab_clean %>%
  filter(Metabolite %in% selected_metabolites)

table(Selected_metab_clean_filtered$Metabolite)
# ---- Create Violin + Boxplot for Each Significant Metabolite Separately ----
metabolite_violin_bar_plot <- ggplot(Selected_metab_clean_filtered, aes(x = MW_scores_lifestyle, y = Concentration, fill = MW_scores_lifestyle)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot for distribution
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  labs(
    title = "Significant Metabolites by Lifestyle",
    x = "Lifestyle Category",
    y = "Metabolite Concentration",
    fill = "Lifestyle"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  scale_fill_manual(values = c("Urban" = "#8B0000", "Non_Urban" = "#95ccba")) +
  facet_wrap(~Metabolite, scales = "free_y")  # Plot each metabolite separately

# ---- Print the Plot ----
print(metabolite_violin_bar_plot)

# ---- Save the Plot ----
ggsave("C20_Metabolites_ViolinBoxPlot_Significance.png", plot = metabolite_violin_bar_plot, width = 10, height = 6, dpi = 300)
##


for_kahumbu_updated <- merge(for_kahumbu, batchMap_with_Meta, by = "Unique.ID", all = T)
colnames(for_kahumbu_updated)
for_kahumbu_updated <- select(for_kahumbu_updated, c(1,2,5,10,11))
names(for_kahumbu_updated)[2] <- "NMR_run_number"
names(for_kahumbu_updated)[3] <- "Mass_spec_run_number"
write.csv(for_kahumbu_updated, "for_kahumbu_updated.csv", row.names = FALSE)
# Load libraries
library(ggplot2)
library(dplyr)
# ---- Step 1: Create Fake Data ----
lipoprotein_data <- data.frame(
  Lifestyle = rep(c("Lifestyle A", "Lifestyle B"), each = 4),
  Subtype = rep(c("S", "M", "L", "XL"), times = 2),
  Value = c(40, 25, 20, 15,   # Lifestyle A
            20, 40, 30, 10)   # Lifestyle B
)
# ---- Step 2: Calculate Proportions ----
lipoprotein_data <- lipoprotein_data %>%
  group_by(Lifestyle) %>%
  mutate(
    Percentage = Value / sum(Value) * 100,
    Label = paste0(Subtype, "\n", round(Percentage, 1), "%"))
# Define a custom purple-green color palette for subtypes
custom_colors <- c(
  "XL" = "#4B0082",   # Indigo (deep purple)
  "L"  = "#9370DB",   # Medium Purple
  "M"  = "#6B8E23",   # Olive Drab (green tone)
  "S"  = "#32CD32"    # Lime Green
)
# ---- Step 3: Plot as Side-by-Side Pie Charts ----
pie_chart <- ggplot(lipoprotein_data, aes(x = "", y = Percentage, fill = Subtype)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4, color = "white") +
  facet_wrap(~Lifestyle) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Lipoprotein X Subtype Proportions by Lifestyle") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#4B0082"),
    strip.text = element_text(size = 14, face = "bold", color = "#006400"),  # Dark green
    legend.position = "none")
# Save plot
ggsave("Lipoprotein_PieChart_PurpleGreen.png", plot = pie_chart, width = 10, height = 6, dpi = 300)


#######Corelation
##Coherence
#Decoherence
# Load required libraries
library(tidyverse)
library(janitor)

# Step 1: Clean column names
WIDE_data_clean <- WIDE_dataSub_with_metadata %>%
  janitor::clean_names()

# Step 2: Extract metabolite data (columns 2 to 629)
metabolite_data <- WIDE_data_clean[, 2:629]
colnames(WIDE_data_clean)
# Step 3: Convert all metabolite columns to numeric (silently handle any coercion warnings)
metabolite_data <- metabolite_data %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))

# Step 4: Remove problematic columns (all NA/Inf or zero variance)
metabolite_data <- metabolite_data %>%
  select(where(~ sum(is.finite(.)) > 1)) %>%
  select(where(~ sd(., na.rm = TRUE) > 0))
# Step 5: Compute correlation matrix across metabolites
cor_matrix <- cor(metabolite_data, use = "pairwise.complete.obs", method = "pearson")

# Step 6: Convert to long format and extract top correlated pairs
cor_df <- as.data.frame(as.table(cor_matrix)) %>%
  rename(met1 = Var1, met2 = Var2, cor = Freq) %>%
  filter(met1 != met2) %>%
  mutate(pair = map2_chr(met1, met2, ~ paste(sort(c(.x, .y)), collapse = "_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(desc(abs(cor))) %>%
  slice_head(n = 30)

# Step 7: Plot top 30 correlated pairs
ggplot(cor_df, aes(x = reorder(pair, abs(cor)), y = cor, fill = cor)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Top 30 Correlated Metabolite Pairs Across Individuals",
       x = "Metabolite Pair", y = "Pearson Correlation") +
  theme_minimal()

# Convert to long format
cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  rename(met1 = Var1, met2 = Var2, cor = Freq)

# Plot full heatmap
ggplot(cor_long, aes(x = met1, y = met2, fill = cor)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       name = "Pearson\nCorrelation") +
  labs(title = "Full Correlation Heatmap of Metabolites",
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Optional: turn off axis text if too many metabolites
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

library(igraph)
library(ggraph)


# Set threshold for strong negative correlation
neg_threshold <- -0.6

# Filter for strongly negative correlations
cor_filtered <- cor_long %>%
  filter(met1 != met2, cor < neg_threshold) %>%
  mutate(pair = map2_chr(met1, met2, ~ paste(sort(c(.x, .y)), collapse = "_"))) %>%
  distinct(pair, .keep_all = TRUE)

# Create graph object
graph_data <- graph_from_data_frame(d = cor_filtered, directed = FALSE)

# Plot negative correlation network
ggraph(graph_data, layout = "fr") +
  geom_edge_link(aes(edge_alpha = abs(cor), edge_width = abs(cor), edge_colour = cor)) +
  geom_node_point(color = "black", size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = paste0("Strong Negative Metabolite Correlations (r < ", neg_threshold, ")")) +
  theme_void()

########@#################Between Lifestyles
library(tidyverse)
library(janitor)
library(reshape2)

library(tidyverse)
library(janitor)
library(reshape2)

# Clean column names
WIDE_data_clean <- WIDE_dataSub_with_metadata %>% janitor::clean_names()

# Metabolite columns
met_cols <- 2:629

# Correlation function
get_correlation_long <- function(df) {
  metabolite_data <- df[, met_cols] %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
    select(where(~ sum(is.finite(.)) > 1)) %>%
    select(where(~ sd(., na.rm = TRUE) > 0))
  
  cor_matrix <- cor(metabolite_data, use = "pairwise.complete.obs", method = "pearson")
  
  as.data.frame(as.table(cor_matrix)) %>%
    rename(met1 = Var1, met2 = Var2, cor = Freq) %>%
    filter(met1 != met2) %>%
    mutate(pair = map2_chr(met1, met2, ~ paste(sort(c(.x, .y)), collapse = "_"))) %>%
    distinct(pair, .keep_all = TRUE)
}

# Split groups
urban_df <- WIDE_data_clean %>% filter(mw_scores_lifestyle == "Urban")
nonurban_df <- WIDE_data_clean %>% filter(mw_scores_lifestyle == "Non_Urban")

# Get correlations
urban_cor <- get_correlation_long(urban_df) %>% rename(cor_urban = cor)
nonurban_cor <- get_correlation_long(nonurban_df) %>% rename(cor_nonurban = cor)

# Merge and fix names
merged_cor <- full_join(urban_cor, nonurban_cor, by = "pair") %>%
  drop_na(cor_urban, cor_nonurban) %>%
  select(pair, met1 = met1.x, met2 = met2.x, cor_urban, cor_nonurban)

# Fisher z-transformation test
fisher_z <- function(r) 0.5 * log((1 + r)/(1 - r))

merged_cor <- merged_cor %>%
  mutate(
    z_urban = fisher_z(cor_urban),
    z_nonurban = fisher_z(cor_nonurban),
    n_urban = nrow(urban_df),
    n_nonurban = nrow(nonurban_df),
    se_diff = sqrt(1 / (n_urban - 3) + 1 / (n_nonurban - 3)),
    z_diff = (z_urban - z_nonurban) / se_diff,
    p_value = 2 * pnorm(-abs(z_diff)),
    p_adj = p.adjust(p_value, method = "BH")  # FDR correction
  ) %>%
  arrange(p_adj)

# View significant results (adjusted p < 0.05)
significant_differences <- merged_cor %>% filter(p_adj < 0.05)

# Output clean results
print(significant_differences %>% select(met1, met2, cor_urban, cor_nonurban, z_diff, p_value, p_adj))

library(ggplot2)

# Volcano plot: difference in correlation (z_diff) vs significance (-log10 p)
ggplot(merged_cor, aes(x = z_diff, y = -log10(p_adj))) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted", color = "darkgrey") +
  labs(
    title = "Volcano Plot: Metabolite Correlation Differences",
    x = "Fisher z Difference (Urban vs Non-Urban)",
    y = "-log10(FDR-adjusted p-value)"
  ) +
  theme_minimal()
library(ggrepel)

top_labels <- merged_cor %>%
  arrange(p_adj) %>%
  slice_head(n = 10) %>%
  mutate(label = paste(met1, met2, sep = "~"))

ggplot(merged_cor, aes(x = z_diff, y = -log10(p_adj))) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_labels, aes(label = label), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted", color = "darkgrey") +
  labs(
    title = "Volcano Plot: Metabolite Correlation Differences",
    x = "Fisher z Difference (Urban vs Non-Urban)",
    y = "-log10(FDR-adjusted p-value)"
  ) +
  theme_minimal()

# Optional: save
write.csv(significant_differences, "significant_cor_differences_by_lifestyle.csv", row.names = FALSE)


library(tidyverse)
library(lme4)
library(broom.mixed)   # tidy.lmerMod
library(performance)   # r2() & AIC for lmer

library(tidyverse)
library(purrr)
library(lme4)
library(performance)

metab_cols  <- 2:629
metab_names <- colnames(WIDE_dataSub_with_metadata)[metab_cols]

fit_compare <- map_dfr(metab_names, function(metab) {
  
  d <- WIDE_dataSub_with_metadata %>%
    select(indiv.ID, batch, Age, Sex,
           MW_scores_lifestyle, MW_scores_lifestyle_2,
           !!sym(metab)) %>%
    filter(if_all(everything(), ~ !is.na(.)))      # complete rows
  
  # skip if too few rows or no/NA variance
  v <- var(d[[metab]], na.rm = TRUE)
  if (nrow(d) < 6 || is.na(v) || v == 0) return(NULL)
  
  # quote metabolite so special chars are legal in formula
  resp <- paste0("`", metab, "`")
  f3   <- paste(resp, "~ MW_scores_lifestyle   + Sex + Age + (1|batch)")
  f2   <- paste(resp, "~ MW_scores_lifestyle_2 + Sex + Age + (1|batch)")
  
  m3 <- lmer(formula = as.formula(f3), data = d, REML = FALSE)
  m2 <- lmer(formula = as.formula(f2), data = d, REML = FALSE)
  
  tibble(
    Metabolite   = metab,
    AIC_3level   = AIC(m3),
    AIC_2level   = AIC(m2),
    delta_AIC    = AIC(m2) - AIC(m3),          # >0 → 3-level fits better
    R2m_3level   = performance::r2(m3)$R2_marginal,
    R2m_2level   = performance::r2(m2)$R2_marginal
  )
})

# quick console summary
fit_compare %>% 
  summarise(n_metab   = n(),
            prefer_3  = sum(delta_AIC >  2),
            prefer_2  = sum(delta_AIC < -2),
            mean_dAIC = mean(abs(delta_AIC))) %>% 
  print(width = Inf)

library(broom)

anova_tbl <- map_dfr(metab_names, function(metab) {
  
  d <- WIDE_dataSub_with_metadata %>% 
    dplyr::select(Metabolite = !!sym(metab),
                  MW_scores_lifestyle, Age, Sex, batch) %>% 
    drop_na()
  
  v <- var(d$Metabolite, na.rm = TRUE)
  if (nrow(d) < 6 || is.na(v) || v == 0) return(NULL)
  
  # mixed model: batch as (random) intercept
  fit <- lmer(Metabolite ~ MW_scores_lifestyle + Age + Sex + (1 | batch), 
              data = d, REML = FALSE)
  
  # lmerTest::anova() supplies Type-III F-tests with Satterthwaite df
  tidy(anova(fit)) %>% 
    filter(term == "MW_scores_lifestyle") %>% 
    transmute(Metabolite = metab,
              p_value    = p.value,
              F_stat     = statistic)
}) %>% arrange(desc(F_stat))

# display the top few hits in console
print(head(anova_tbl, 10), width = Inf)


# ── Bonferroni-corrected hits, keep 50 biggest effects ──────────────
n_tests <- nrow(anova_tbl)                                    # total models

anova_tbl <- anova_tbl %>% 
  mutate(p_adj = pmin(1, p_value * n_tests)) %>%              # Bonferroni
  filter(p_adj < 0.05) %>%                                    # pass threshold
  arrange(desc(F_stat)) %>%                                   # strongest first
  slice_head(n = 100)                                          # top 50 only

sig_metabs <- anova_tbl$Metabolite                            # vector for plotting

if (length(sig_metabs) == 0) {
  stop("No metabolites survive the Bonferroni correction.")
}

library(stringr)   # for str_extract()

get_sig_bars <- function(metab) {
  d <- WIDE_dataSub_with_metadata %>% 
    dplyr::select(Metabolite = !!sym(metab),
                  MW_scores_lifestyle, Age, Sex, batch) %>% 
    drop_na()
  
  fit  <- lm(Metabolite ~ MW_scores_lifestyle + Age + Sex + (1|batch), data = d)
  emms <- emmeans(fit, "MW_scores_lifestyle")
  tuk  <- pairs(emms) %>% as_tibble()
  
  out <- tuk %>%
    filter(p.value < 0.05) %>%
    mutate(
      group1     = str_extract(contrast, "^[^ ]+"),
      group2     = str_extract(contrast, "[^ ]+$"),
      Metabolite = metab,
      label      = case_when(
        p.value < .001 ~ "***",
        p.value < .01  ~ "**",
        TRUE           ~ "*")
    ) %>%
    dplyr::select(Metabolite, group1, group2, p.value, label)
  
  # if no rows, return an empty tibble with correct column classes
  if (nrow(out) == 0) {
    out <- tibble(
      Metabolite = character(),
      group1     = character(),
      group2     = character(),
      p.value    = numeric(),
      label      = character()
    )
  }
  out
}

sig_bars <- purrr::map_dfr(sig_metabs, get_sig_bars)

# y-position for bars (panel-specific max)
y_pos_tbl <- WIDE_dataSub_with_metadata %>% 
  dplyr::select(all_of(sig_metabs), ind = indiv.ID) %>%
  pivot_longer(-ind, names_to = "Metabolite", values_to = "Abundance") %>%
  group_by(Metabolite) %>%
  summarise(y_pos = max(Abundance, na.rm = TRUE) * 1.07, .groups = "drop")

sig_bars <- sig_bars %>% left_join(y_pos_tbl, by = "Metabolite")


wes_cols <- wesanderson::wes_palette("Darjeeling1", 3)

# long data for all significant metabolites
plot_df <- WIDE_dataSub_with_metadata %>%
  dplyr::select(indiv.ID, MW_scores_lifestyle, all_of(sig_metabs)) %>%
  pivot_longer(-c(indiv.ID, MW_scores_lifestyle),
               names_to  = "Metabolite",
               values_to = "Abundance")

# ── open multipage PDF ──────────────────────────────────────────
pdf("Significant_metabolites_boxplots.pdf", width = 14, height = 10)

for (chunk in split(sig_metabs, ceiling(seq_along(sig_metabs) / 16))) {
  
  df_chunk   <- plot_df  %>% filter(Metabolite %in% chunk)
  bars_chunk <- sig_bars %>% filter(Metabolite %in% chunk)
  
  # panel-specific y max (update every chunk)
  ylim_tbl <- df_chunk %>%
    group_by(Metabolite) %>%
    summarise(y_pos = max(Abundance, na.rm = TRUE) * 1.07, .groups = "drop")
  bars_chunk <- bars_chunk %>%
    dplyr::select(-y_pos) %>% left_join(ylim_tbl, by = "Metabolite")
  
  p <- ggplot(df_chunk,
              aes(MW_scores_lifestyle, Abundance,
                  fill = MW_scores_lifestyle)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
    
    # —— significance bars ————————————————
    ggpubr::stat_pvalue_manual(
      bars_chunk,
      label        = "label",      # column with asterisks
      xmin         = "group1",
      xmax         = "group2",
      y.position   = "y_pos",
      tip.length   = 0.01,
      size         = 3,
      bracket.size = 0.4,
      inherit.aes  = FALSE          # <-- KEY: don’t inherit MW_scores_lifestyle
    ) +
    
    facet_wrap(~ Metabolite, scales = "free_y",
               ncol = 4, nrow = 4) +
    scale_fill_manual(values = wes_cols) +
    labs(title = "Metabolites with FDR-adjusted lifestyle effect",
         x = "Lifestyle group", y = "Abundance") +
    theme_minimal(base_size = 10) +
    theme(strip.text      = element_text(size = 8),
          legend.position = "none")
  
  print(p)    #one PDF page per 16 metabolites
}

dev.off()
message("PDF saved: Significant_metabolites_boxplots.pdf")


library(ggplot2)
library(dplyr)
library(broom)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(ggpubr)

# Metabolites of interest and renaming for clarity
ceramides <- c("ceramide (d18:1/22:1)", 
               "ceramide (d18:1/24:1)", 
               "ceramide(d18:1/16:0)")
names_map <- c("ceramide (d18:1/22:1)" = "Ceramide_22_1",
               "ceramide (d18:1/24:1)" = "Ceramide_24_1",
               "ceramide(d18:1/16:0)" = "Ceramide_16_0")

# Subset and pivot data to long format
long_df <- WIDE_dataSub_with_metadata %>%
  dplyr::select(MW_scores_lifestyle, batch, Age, Sex, all_of(ceramides)) %>%
  pivot_longer(cols = all_of(ceramides), 
               names_to = "Metabolite", 
               values_to = "Abundance") %>%
  mutate(Metabolite = recode(Metabolite, !!!names_map))

# Remove outliers using IQR method per metabolite
long_df <- long_df %>%
  group_by(Metabolite) %>%
  filter({
    q1 <- quantile(Abundance, 0.25, na.rm = TRUE)
    q3 <- quantile(Abundance, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    Abundance >= (q1 - 1.5 * iqr) & Abundance <= (q3 + 1.5 * iqr)
  }) %>%
  ungroup()

# Run ANOVA controlling for covariates per metabolite
anova_results <- long_df %>%
  group_by(Metabolite) %>%
  do({
    model <- aov(Abundance ~ MW_scores_lifestyle + batch + Age + Sex, data = .)
    broom::tidy(model) %>%
      filter(term == "MW_scores_lifestyle") %>%
      transmute(F_stat = statistic, p_val = p.value)
  })

# Add stats back to plotting dataframe
long_df <- long_df %>%
  left_join(anova_results, by = "Metabolite")

# Plot violin plots with square facets
violin_plot <- ggplot(long_df, aes(x = MW_scores_lifestyle, y = Abundance, fill = MW_scores_lifestyle)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~ Metabolite, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    aspect.ratio = 1,  # Makes square facets
    panel.spacing = unit(1, "lines")
  ) +
  labs(x = "Lifestyle", y = "Metabolite Abundance") +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "black") +
  geom_text(
    aes(label = paste0("F=", round(F_stat, 2), ", p=", signif(p_val, 2))),
    x = 1.5, y = Inf, vjust = 1.5, hjust = 0,
    inherit.aes = FALSE,
    data = distinct(long_df, Metabolite, F_stat, p_val)
  )

# Show plot
print(violin_plot)

