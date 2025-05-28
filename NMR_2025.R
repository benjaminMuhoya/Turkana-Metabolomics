setwd("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/")
AMAMAMA <- import("Metabs_Python_MCCV.csv")
dim(AMAMAMA)
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rio)
library(gridExtra)
library(grid) 
library(reshape2)
library(tidyverse)
library(tidyr)
library(broom)
library(lme4)
library(broom.mixed)
library(lmerTest)
library(ggpubr)
library(ggfortify)
library(patchwork)
library(wesanderson)


# Read All_data file (Meta_ONLY)
All_data <- import("THGP_database_TurkanaOnly_corrected_merged_2025-01-28.txt")
##Subset data to only the needed columns
columns_to_keep <- c(1, 25, 31, 138, 172, 173, 175, 179, 181, 182, 203, 205, 206, 230, 232, 235, 239, 241, 242, 268:271, 272, 277, 278, 284, 289, 290, 297, 298, 302, 315, 368, 369, 445, 489:492, 499:509, 601:604, 606:609, 617, 620, 621, 639, 644:658, 680:929, 115, 403, 404)
# Subset the data to include only the selected columns
RAW_Data_for_NMR_ANALYSIS <- All_data %>% dplyr::select(all_of(columns_to_keep))
table(RAW_Data_for_NMR_ANALYSIS$MW_scores_h_sol_trim)
# Rename columns by removing "Nightingale_NMR_" prefix
colnames(RAW_Data_for_NMR_ANALYSIS) <- gsub("^Nightingale_NMR_", "", colnames(RAW_Data_for_NMR_ANALYSIS))
# Remove those ending with _pct & .PCT i.e. Percentages. 
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  dplyr::select(-matches("_pct$"),  -matches("(?i)\\.PCT$"))
# -------------------------------- Process the columns for analysis --------------------------------------
colnames(RAW_Data_for_NMR_ANALYSIS)
# Process all the Numeric ones - scale and remove outliers
columns_to_process <- c(2,36:57,59,74:246) # Age and h_sol included
RAW_Data_for_NMR_ANALYSIS$LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml[RAW_Data_for_NMR_ANALYSIS$LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml < 0] <- NA
# Factor Columns
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  mutate(Sex = ifelse(Sex %in% c("Male", "Female"), Sex, NA))
RAW_Data_for_NMR_ANALYSIS$Sex <- as.factor(RAW_Data_for_NMR_ANALYSIS$Sex)
table(RAW_Data_for_NMR_ANALYSIS$Sex)
RAW_Data_for_NMR_ANALYSIS$LBP_ELISA_Batch <- as.factor(RAW_Data_for_NMR_ANALYSIS$LBP_ELISA_Batch)
RAW_Data_for_NMR_ANALYSIS$MW_scores_lifestyle <- as.factor(RAW_Data_for_NMR_ANALYSIS$MW_scores_lifestyle)
#####################################################################################################
# Ensure lifestyle_group is properly classified
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  mutate(lifestyle_group_MW = ifelse(MW_scores_lifestyle %in% c("Pastoralist", "PeriUrban"),
                                  "Non_Urban", "Urban"))
table(RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW)
table(RAW_Data_for_NMR_ANALYSIS$MW_scores_lifestyle)
# Function to convert columns to numeric
convert_to_numeric <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    df[[column_name]] <- as.numeric(as.character(df[[column_name]]))
  }
  return(df)
}
RAW_Data_for_NMR_ANALYSIS <- convert_to_numeric(RAW_Data_for_NMR_ANALYSIS, columns_to_process)
table(RAW_Data_for_NMR_ANALYSIS$MW_scores_h_sol_trim)
# -------- Function to remove outliers RAW data -------------
remove_outliers_iqr <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    if (is.numeric(df[[column_name]])) {
      Q1 <- quantile(df[[column_name]], 0.25, na.rm = TRUE)
      Q3 <- quantile(df[[column_name]], 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      # Replace outliers with NA
      df[[column_name]][df[[column_name]] < lower_bound | df[[column_name]] > upper_bound] <- NA
    }
  }
  return(df)
}
RAW_Data_for_NMR_ANALYSIS <- remove_outliers_iqr(RAW_Data_for_NMR_ANALYSIS, columns_to_process)
# Standardize and normalize columns
standardize_columns <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    df[[column_name]] <- scale(df[[column_name]], center = TRUE, scale = TRUE)
  }
  return(df)
}
Data_for_NMR_ANALYSIS <- log_transform_and_standardize(RAW_Data_for_NMR_ANALYSIS, columns_to_process)
table(Data_for_NMR_ANALYSIS$MW_scores_h_sol_trim)
# -------- Function to remove outliers and replace with NA (transformed data) -------------
remove_outliers_iqr <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    if (is.numeric(df[[column_name]])) {
      Q1 <- quantile(df[[column_name]], 0.25, na.rm = TRUE)
      Q3 <- quantile(df[[column_name]], 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      # Replace outliers with NA
      df[[column_name]][df[[column_name]] < lower_bound | df[[column_name]] > upper_bound] <- NA
    }
  }
  return(df)
}
Data_for_NMR_ANALYSIS <- remove_outliers_iqr(Data_for_NMR_ANALYSIS, columns_to_process)
table(Data_for_NMR_ANALYSIS$MW_scores_h_sol_trim)
dim(Data_for_NMR_ANALYSIS)
#####################################################################################################
# Checking for Normality
# Drop NA from batch column
Data_for_NMR_ANALYSIS_filtered <- Data_for_NMR_ANALYSIS %>% drop_na(LBP_ELISA_Batch)
# ---- Run Shapiro test and store results ----
normality_results <- data.frame()
for (metab in metabolite_cols) {
  current_data <- Data_for_NMR_ANALYSIS_filtered %>%
    dplyr::select(all_of(metab)) %>%
    drop_na()
  values <- current_data[[1]]
  if (length(values) >= 4 && length(unique(values)) > 1) {
    shapiro <- shapiro.test(values)
    
    normality_results <- rbind(normality_results, data.frame(
      Metabolite = metab,
      N = length(values),
      W_statistic = as.numeric(shapiro$statistic),
      p_value = as.numeric(shapiro$p.value),
      Deviates_from_normal = ifelse(shapiro$p.value < 0.05, "Yes", "No")
    ))
  } else {
    reason <- if (length(values) < 4) "Insufficient data" else "Identical values"
    normality_results <- rbind(normality_results, data.frame(
      Metabolite = metab,
      N = length(values),
      W_statistic = NA,
      p_value = NA,
      Deviates_from_normal = reason
    ))
  }
}
# ---- Filter to usable results and order by deviation ----
plot_data <- normality_results %>%
  filter(Deviates_from_normal %in% c("Yes", "No")) %>%
  mutate(order_priority = ifelse(Deviates_from_normal == "Yes", 1, 2)) %>%
  arrange(order_priority, p_value)
# ---- Create annotated histograms ----
pdf("Metabolite_Histograms_Annotated.pdf", width = 14, height = 10)
plot_list <- list()
count <- 0
for (i in seq_len(nrow(plot_data))) {
  metab <- plot_data$Metabolite[i]
  values <- Data_for_NMR_ANALYSIS_filtered[[metab]]
  values <- values[!is.na(values)]
  # Build label
  label_text <- paste0(
    "N = ", plot_data$N[i],
    ", W = ", signif(plot_data$W_statistic[i], 4),
    ", p = ", signif(plot_data$p_value[i], 3),
    "\nDeviates from Normal (Shapiro.test)? ", plot_data$Deviates_from_normal[i])
  
  p <- ggplot(data.frame(x = values), aes(x = x)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = metab, subtitle = label_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.subtitle = element_text(size = 8))
  
  plot_list[[length(plot_list) + 1]] <- p
  count <- count + 1
  if (count %% 25 == 0) {
    do.call("grid.arrange", c(plot_list, ncol = 5))
    plot_list <- list()
  }
}
# Print final page
if (length(plot_list) > 0) {
  do.call("grid.arrange", c(plot_list, ncol = 5))
}
dev.off()
#####################################################################################################
# Define metabolite columns
metabolite_cols <- names(Data_for_NMR_ANALYSIS[, c(37:57, 59, 74:246)]) ##Exclude inflammation columns
# Define non-normal metabolites to skip
non_normal_metabolites <- c(
  "XXL_VLDL_P", "XXL_VLDL_L", "XXL_VLDL_PL", "XXL_VLDL_C", "XXL_VLDL_CE",
  "XXL_VLDL_FC", "XXL_VLDL_TG", "XL_VLDL_P", "XL_VLDL_L", "XL_VLDL_PL",
  "XL_VLDL_FC", "XL_VLDL_TG", "Colombia_C_reactive_protein..mg.L.")
# Filter to only normal metabolites
normal_metabolites <- setdiff(metabolite_cols, non_normal_metabolites)
## Factor Columns
Data_for_NMR_ANALYSIS$Sex <- as.factor(Data_for_NMR_ANALYSIS$Sex)
Data_for_NMR_ANALYSIS$LBP_ELISA_Batch <- as.factor(Data_for_NMR_ANALYSIS$LBP_ELISA_Batch)
Data_for_NMR_ANALYSIS$lifestyle_group <- as.factor(Data_for_NMR_ANALYSIS$lifestyle_group)
Data_for_NMR_ANALYSIS$lifestyle_group_MW <- relevel(factor(Data_for_NMR_ANALYSIS$lifestyle_group_MW), ref = "Non_Urban")
#####################################################################################################
## Trying Different Lifestyle Groups
Data_for_NMR_ANALYSIS_filtered$lifestyle_group <- factor(
  ifelse(Data_for_NMR_ANALYSIS_filtered$MW_scores_h_sol_trim <= 1, "one_and_below", "above_one"),
  levels = c("one_and_below", "above_one")
)
table(Data_for_NMR_ANALYSIS_filtered$lifestyle_group_MW)
Data_for_NMR_ANALYSIS_filtered$lifestyle_group_MW <- relevel(factor(Data_for_NMR_ANALYSIS_filtered$lifestyle_group_MW), ref = "Non_Urban")
table(Data_for_NMR_ANALYSIS_filtered$lifestyle_group)
Data_for_NMR_ANALYSIS$lifestyle_group <- as.factor(Data_for_NMR_ANALYSIS$lifestyle_group)
####################################################################################################
# Calculate percentage of missing values per metabolite
missing_summary <- Data_for_NMR_ANALYSIS_filtered %>%
  select(all_of(normal_metabolites)) %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(cols = everything(),
               names_to = "metabolite",
               values_to = "percent_missing")
# Filter out metabolites with <10% missingness
metabolites_to_keep <- missing_summary %>%
  filter(percent_missing < 10) %>% pull(metabolite)
# Subset data and drop rows with remaining missing values
pca_data <- Data_for_NMR_ANALYSIS_filtered %>%
  select(all_of(metabolites_to_keep)) %>% drop_na()
# Match metadata
metadata <- Data_for_NMR_ANALYSIS_filtered %>%
  filter(rowSums(is.na(select(., all_of(metabolites_to_keep)))) == 0) %>%
  select(Unique.ID, lifestyle_group_MW, Sex, LBP_ELISA_Batch)
# Confirm dimensions
cat("PCA input dimensions:", dim(pca_data)[1], "samples ×", dim(pca_data)[2], "metabolites\n")
# Run PCA
pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)
# Get 5 colors from FantasticFox1 and add black
custom_colors <- c(wes_palette("Darjeeling1", n = 5), "black")
# Now use it in scale_color_manual
wes_colors <- scale_color_manual(values = custom_colors)
# Plot 1: Colored by lifestyle group
p1 <- autoplot(pca_res, data = metadata, colour = 'lifestyle_group_MW') +
  wes_colors +
  theme_minimal() +
  labs(title = "Colored by Lifestyle Group", x = "PC1", y = "PC2") +
  theme(text = element_text(size = 11))
# Plot 2: Colored by sex
p2 <- autoplot(pca_res, data = metadata, colour = 'Sex') +
  wes_colors +
  theme_minimal() +
  labs(title = "Colored by Sex", x = "PC1", y = "PC2") +
  theme(text = element_text(size = 11))
# Plot 3: Colored by batch
p3 <- autoplot(pca_res, data = metadata, colour = 'LBP_ELISA_Batch') +
  wes_colors +
  theme_minimal() +
  labs(title = "Colored by Batch", x = "PC1", y = "PC2") +
  theme(text = element_text(size = 11))
# Combine the three plots
combined_pca_plot <- (p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.position = "top")

# Save as PNG
ggsave("PCA_by_Lifestyle_Sex_Batch.png",plot = combined_pca_plot, width = 16, height = 6, dpi = 300)

#####################################################################################################
# Initialize dataframes to store results
all_coefficients <- data.frame()
# Loop through each metabolite and run the linear model
for (metabolite in normal_metabolites) {
  # Subset relevant columns and drop NAs only from predictors and the metabolite
  current_data <- Data_for_NMR_ANALYSIS %>%
    dplyr::select(lifestyle_group_MW, Age, Sex, all_of(metabolite)) %>%
    rename(Metabolite = all_of(metabolite)) %>%
    drop_na(Metabolite, Age, Sex)
  # Ensure enough data to run the model
  if (nrow(current_data) > 6) {
    model <- lm(Metabolite ~ lifestyle_group_MW + Sex + Age + lifestyle_group_MW:Sex,
                data = current_data)
    temp_results <- broom::tidy(model) %>%
      mutate(Metabolite = metabolite, N = nrow(current_data))
    all_coefficients <- bind_rows(all_coefficients, temp_results)
  } else {
    message(paste("Skipping metabolite:", metabolite, "- insufficient data"))
  }
}

# Compute confidence intervals
all_coefficients <- all_coefficients %>%
  mutate(lower_CI = estimate - 1.96 * std.error,
    upper_CI = estimate + 1.96 * std.error)
# Apply Bonferroni correction to p-values
all_coefficients <- all_coefficients %>%
  mutate(adjusted_p_values = p.adjust(p.value, method = "bonferroni"))
head(all_coefficients)
#########################
# Pivot and Save to CSV
formatted_results <- all_coefficients %>%
  filter(term != "(Intercept)") %>%
  select(Metabolite, term, estimate, std.error, p.value, adjusted_p_values, N) %>%
  pivot_wider(
    names_from = term,
    values_from = c(estimate, std.error, p.value, adjusted_p_values, N),
    names_glue = "{.value}_{term}")
# Save to CSV
write_csv(formatted_results, "NMR_all_coefficients_results_first_model.csv")


# --------------------- Plot with estimates and confidence intervals -------------------
# Compute confidence intervals and adjust p-values
all_coefficients <- all_coefficients %>%
  mutate(
    lower_CI = estimate - 1.96 * std.error,
    upper_CI = estimate + 1.96 * std.error,
    adjusted_p_values = p.adjust(p.value, method = "bonferroni")
  )
# Filter for significant and strong effects (|β| ≥ 0.4 & CI excludes 0)
filtered_coeff_results <- all_coefficients %>%
  filter(abs(estimate) >= 0.5 & (lower_CI > 0 | upper_CI < 0) & term != "(Intercept)") %>%
  arrange(desc(abs(estimate)))
# Generate the coefficient plot
coef_plot <- ggplot(filtered_coeff_results, aes(x = Metabolite, y = estimate, ymin = lower_CI, ymax = upper_CI, color = term)) +
  geom_pointrange(alpha = 0.7, size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(aes(label = N), hjust = -0.3, size = 3.5, color = "black") + theme_minimal() +
  scale_color_manual(values = c(
    "lifestyle_groupabove_one" = "blue",
    "SexMale" = "red",
    "Age" = "green",
    "MW_scores_market_diet_index_trim" = "purple",
    "lifestyle_groupabove_one:SexMale" = "orange",
    "lifestyle_group_MWUrban" = "blue", 
    "lifestyle_group_MWUrban:SexMale" = "orange")) +
  labs(
    title = "Coefficient Plot: Effect of Predictors",
    x = "Metabolites",
    y = "Beta Estimate",
    color = "Predictor",
    caption = "Filtered for |β| ≥ 0.4 & Significant 95% CI") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right") + coord_flip()

# Save and display the plot
ggsave("coefficient_plot.png", plot = coef_plot, device = "png", width = 8, height = 6, dpi = 300)
print(coef_plot)
#####################################################################################################
# ---------------------- Visualize Distributions of specific columns of interest ------------------------
# ---------------------- Inflammation Histogram + ANOVA ------------------------
# Define inflammation biomarkers and display names
inflammation_map <- c(
  "Colombia_C_reactive_protein..mg.L." = "CRP",
  "Colombia_globulin..g.dL." = "Globulin",
  "Colombia_albumin__serum..g.dL." = "Albumin_Kenya_lab",
  "LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml" = "LBP",
  "Albumin" = "Albumin_NMR", "GlycA" = "GlycA"
)

# Define custom fill colors for each marker
custom_colors <- c(
  "CRP" = "green", "Globulin" = "gold", "Albumin_Kenya_lab" = "pink",
  "LBP" = "red", "Albumin_NMR" = "steelblue", "GlycA" = "plum1"
)

# Fix LBP values
RAW_Data_for_NMR_ANALYSIS$LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml[
  RAW_Data_for_NMR_ANALYSIS$LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml < 0
] <- NA

# -------------------- Histogram Loop --------------------
for (raw_col in names(inflammation_map)) {
  display_name <- inflammation_map[[raw_col]]
  if (!raw_col %in% colnames(RAW_Data_for_NMR_ANALYSIS)) next
  
  subset_data <- RAW_Data_for_NMR_ANALYSIS %>%
    dplyr::select(all_of(raw_col)) %>%
    drop_na()
  
  if (nrow(subset_data) < 10) next
  
  p <- ggplot(subset_data, aes(x = .data[[raw_col]])) +
    geom_histogram(bins = 30, fill = custom_colors[display_name], color = "white", alpha = 0.8) +
    labs(
      title = paste("Distribution of", display_name),
      x = display_name,
      y = "Frequency"
    ) +
    scale_x_continuous(limits = c(0, NA)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  ggsave(
    filename = paste0("Histogram_", display_name, ".png"),
    plot = p,
    width = 4, height = 4, dpi = 300
  )
}
#
# -------------------- ANOVA or Kruskal-Wallis Loop --------------------
# Initialize results
inflammation_results <- data.frame()
#########
for (raw_col in names(inflammation_map)) {
  name <- inflammation_map[[raw_col]]
  df <- RAW_Data_for_NMR_ANALYSIS %>%
    select(lifestyle_group_MW, Age, Sex, !!sym(raw_col)) %>%
    drop_na()
  if (nrow(df) < 10) next
  # Run appropriate test
  if (shapiro.test(df[[raw_col]])$p.value < 0.05) {
    test <- kruskal.test(as.formula(paste(raw_col, "~ lifestyle_group_MW")), data = df)
    stat <- unname(test$statistic)
    pval <- test$p.value
    type <- "Kruskal-Wallis"
  } else {
    mod <- lm(as.formula(paste(raw_col, "~ lifestyle_group_MW + Age + Sex")), data = df)
    aov_out <- anova(mod)
    stat <- aov_out["lifestyle_group_MW", "F value"]
    pval <- aov_out["lifestyle_group_MW", "Pr(>F)"]
    type <- "ANOVA"
  }
  inflammation_results <- bind_rows(inflammation_results, tibble(
    Molecule = name, Test = type, N = nrow(df),
    Statistic = signif(stat, 4),
    p_value = signif(pval, 4)
  ))
}
#########
write_csv(inflammation_results, "inflammation_test_summary.csv")
#########
for (i in seq_len(nrow(inflammation_results))) {
  mol <- inflammation_results$Molecule[i]
  raw_col <- names(inflammation_map)[inflammation_map == mol]
  
  df <- Data_for_NMR_ANALYSIS_filtered %>%
    select(lifestyle_group_MW, !!sym(raw_col)) %>%
    drop_na()
  
  test_type <- inflammation_results$Test[i]
  stat <- inflammation_results$Statistic[i]
  pval <- inflammation_results$p_value[i]
  
  subtitle_text <- paste0(test_type, ": stat = ", stat, ", p = ", signif(pval, 3))
  
  p <- ggplot(df, aes(x = lifestyle_group_MW, y = .data[[raw_col]], fill = lifestyle_group_MW)) +
    geom_boxplot(alpha = 0.6) +
    labs(
      title = mol,
      subtitle = subtitle_text,
      y = mol,
      x = "Lifestyle Group"
    ) +
    theme_minimal()
  
  ggsave(paste0("Boxplot_", mol, ".png"), p, width = 6, height = 4, dpi = 300)
}

#####################################################################################################
# ---------Now I focus on HDL and Inflammation -----------------------------
colnames(Data_for_NMR_ANALYSIS)
# First, define the columns needed for the regression
# Rename specific inflammation columns
HDL_subset_data <- Data_for_NMR_ANALYSIS %>%
  rename(
    CRP = `Colombia_C_reactive_protein..mg.L.`,
    Globulin = `Colombia_globulin..g.dL.`,
    Albumin = `Colombia_albumin__serum..g.dL.`,
    LBP = `LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml`,
    Albumin_NMR = `Albumin`
  )
dim(HDL_subset_data)
# Define inflammation biomarkers
inflammation_cols <- c("CRP", "Globulin", "Albumin", "LBP", "GlycA", "Albumin_NMR")
# Define metabolite columns based on column indices
metabolite_indices <- c(37:46, 48:53, 56:57, 74:146, 148:246) ## All of them 
HDL_indices <- c(38,39,49,50,80,81,84,88,92,96,99,103,106,112, 113,114,218:245)# Includes Triglycerides too
LDL_indices <- c(39,40,49,50,75,78,79,81,83,87,91,95,98,102,105,112,113,114,197:217,181,178,186,155,167,187,177,158,163,179,176,104,77,180,182,159,90,86,94,184,162,183,185,188,189) # Includes Triglycerides too + VLDL that are normally distributed
metabolite_cols <- colnames(HDL_subset_data)[metabolite_indices]
# --------------------------- LM Regression with Inflammation phenotypes ------------------------------
# loop to run linear models and collect results
run_lm_models <- function(inflam, metabolite_cols, data) {
  results <- data.frame()
  for (metab in metabolite_cols) {
    if (!metab %in% colnames(data)) next  # Skip if column missing
    # Select only relevant columns and drop NA from them
    df <- data %>%
      dplyr::select(all_of(c(metab, inflam, "lifestyle_group_MW", "Age", "Sex"))) %>%
      drop_na(.data[[inflam]], .data[[metab]])
    if (nrow(df) < 10) next  # Skip if too few observations
    # Fit model and handle possible errors
    model <- tryCatch(
      lm(as.formula(paste(metab, "~", inflam, "+ lifestyle_group_MW + lifestyle_group_MW:Sex + Age + Sex")), data = df),
      error = function(e) NULL)
    if (!is.null(model)) {
      coefs <- broom::tidy(model) %>%
        filter(term == inflam) %>%
        mutate(Metabolite = metab, N = nrow(df)) %>%
        select(Metabolite, Predictor = term, Beta = estimate, p_value = p.value, N)
      results <- bind_rows(results, coefs)
    }
  }
  return(results)
}
##
all_results <- list()
for (inflam in inflammation_cols) {
  cat("Analyzing:", inflam, "\n")
  res <- run_lm_models(inflam, metabolite_cols, HDL_subset_data)
  total_tests <- nrow(res)
  cat("  Total tests:", total_tests, "\n")
  if (total_tests > 0) {
    res$adj_p <- p.adjust(res$p_value, method = "bonferroni")
    cat("  Bonferroni-significant:", sum(res$adj_p < 0.05), "\n")
    all_results[[inflam]] <- res  # ✅ Store *all* results, not just significant ones
  } else {
    all_results[[inflam]] <- data.frame()
  }
}
##
for (inflam in names(all_results)) {
  neg_hits <- all_results[[inflam]] %>%
    filter(Beta < 0 & p_value < 0.05)  # Filter by raw p-value
  
  if (nrow(neg_hits) > 0) {
    cat("\n🔽 Negative associations for", inflam, "(raw p < 0.05):\n")
    print(neg_hits %>% dplyr::select(Metabolite, Predictor, Beta, p_value, adj_p))
  } else {
    cat("\n🔸 No negative associations found for", inflam, "\n")
  }
}
##
##
for (inflam in names(all_results)) {
  neg_hits <- all_results[[inflam]] %>%
    filter(Beta < 0 & adj_p < 0.05)  # Filter by raw p-value
  
  if (nrow(neg_hits) > 0) {
    cat("\n🔽 Negative associations for", inflam, "(raw p < 0.05):\n")
    print(neg_hits %>% dplyr::select(Metabolite, Predictor, Beta, p_value, adj_p))
  } else {
    cat("\n🔸 No negative associations found for", inflam, "\n")
  }
}
# Save each inflammation marker's results as a separate CSV
for (inflamm in names(all_results)) {
  result_df <- all_results[[inflamm]]
  if (nrow(result_df) > 0) {
    file_name <- paste0("LM_results_", inflamm, ".csv")
    write_csv(result_df, file_name)
  }
}
####################################################################################################
####################################################################################################
## -------------------------------------------- Assessing subclass abundances ----------------------
# ---- Select and Clean Relevant Columns ----
## -- Use this line for LDL 
# Replace it downstream the code
##
##
ldl_columns <- grep("^((?!VLDL).)*LDL_", colnames(RAW_Data_for_NMR_ANALYSIS), value = TRUE, perl = TRUE)
##
##
#---- Here I do HDL only ----
hdl_columns <- grep("HDL_", colnames(RAW_Data_for_NMR_ANALYSIS), value = TRUE)
additional_cols <- c(
  "MW_scores_lifestyle", "MW_scores_lifestyle_num", "MW_scores_h_sol", "MW_scores_h_sol_trim", "MW_scores_traditional_diet_score",
  "MW_scores_traditional_diet_score_trim", "MW_scores_market_diet_index", "MW_scores_market_diet_index_trim",
  "Total_C", "Unique.ID", "Age", "Sex", "BMI"
)

# ---- Subset and rename columns ----
HDL_subset_data <- RAW_Data_for_NMR_ANALYSIS %>%
  dplyr::select(all_of(c(ldl_columns, additional_cols)))
colnames(HDL_subset_data)
# ---- Filter and Transform HDL Data ----
hdl_data <- HDL_subset_data %>%
  mutate(MW_scores_lifestyle = ifelse(MW_scores_lifestyle %in% c("Pastoralist", "PeriUrban"), 
                                      "Non_Urban", MW_scores_lifestyle))
# ---- Reshape to Long Format and Annotate ----
hdl_long <- hdl_data %>%
  dplyr::select(MW_scores_lifestyle, matches("^(XL|L|M|S)_L")) %>%
  pivot_longer(cols = -MW_scores_lifestyle, names_to = "HDL_Type", values_to = "Value") %>%
  mutate(
    SubClass = case_when(
      str_starts(HDL_Type, "XL_") ~ "Extra Large",         # HDL
      str_starts(HDL_Type, "L_") & !str_starts(HDL_Type, "L_LDL") ~ "Large",        # HDL
      str_starts(HDL_Type, "M_") & !str_starts(HDL_Type, "M_LDL") ~ "Medium",       # HDL
      str_starts(HDL_Type, "S_") & !str_starts(HDL_Type, "S_LDL") ~ "Small",        # HDL
      str_starts(HDL_Type, "L_LDL") ~ "Large_LDL",         # LDL
      str_starts(HDL_Type, "M_LDL") ~ "Medium_LDL",        # LDL
      str_starts(HDL_Type, "S_LDL") ~ "Small_LDL"          # LDL
    ),
    Particle = case_when(
      str_detect(HDL_Type, "_FC") ~ "Free Cholesterol",
      str_detect(HDL_Type, "_PL") ~ "Phospholipids",
      str_detect(HDL_Type, "_CE") ~ "Cholesterol Ester",
      str_detect(HDL_Type, "_TG") ~ "Triglycerides",
      str_detect(HDL_Type, "_C")  ~ "Cholesterol",
      str_detect(HDL_Type, "_P")  ~ "No. of Particles",
      str_detect(HDL_Type, "_L")  ~ "Total Lipids"
    )
  ) %>%
  filter(!is.na(SubClass))  # Remove non-subclass entries
# ---- Calculate Summary Stats Per Lifestyle ----
hdl_summary <- hdl_long %>%
  group_by(MW_scores_lifestyle, Particle, SubClass) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  group_by(MW_scores_lifestyle, Particle) %>%
  mutate(
    Total_Mean = sum(Mean_Value),
    Percentage = (Mean_Value / Total_Mean) * 100,
    Label = paste0(round(Percentage, 1), "%")
  ) %>%
  ungroup()
# ---- Add Whole Population Summary ----
hdl_whole_summary <- hdl_long %>%
  group_by(Particle, SubClass) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  group_by(Particle) %>%
  mutate(
    Total_Mean = sum(Mean_Value),
    Percentage = (Mean_Value / Total_Mean) * 100,
    Label = paste0(round(Percentage, 1), "%"),
    MW_scores_lifestyle = "WholePopulation"
  ) %>%
  ungroup()
# ---- Combine Summaries and Plot ----
hdl_combined_summary <- bind_rows(hdl_summary, hdl_whole_summary)
# ---- Clean MW_scores_lifestyle ----
hdl_combined_summary <- hdl_combined_summary %>%
  mutate(MW_scores_lifestyle = as.character(MW_scores_lifestyle)) %>%
  filter(!is.na(MW_scores_lifestyle))
# ---- Loop through each lifestyle group and save a pie chart per group ----
unique_lifestyles <- unique(hdl_combined_summary$MW_scores_lifestyle)
for (lifestyle in unique_lifestyles) {
  lifestyle_data <- hdl_combined_summary %>%
    filter(MW_scores_lifestyle == lifestyle)
  # If no data for this group, skip
  if (nrow(lifestyle_data) == 0) next
  # Calculate positions for pie chart
  lifestyle_data <- lifestyle_data %>%
    group_by(Particle) %>%
    arrange(desc(SubClass)) %>%
    mutate(
      ymax = cumsum(Percentage),
      ymin = lag(ymax, default = 0),
      mid = (ymin + ymax) / 2
    ) %>%
    ungroup()
  # Plot
  lifestyle_plot <- ggplot(lifestyle_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = SubClass)) +
    geom_rect(color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(x = 4.5, y = mid, label = Label), size = 3.5)+ scale_fill_manual(values = c(
      "Small" = "#8B5FBF",        # Purple
      "Medium" = "plum1",         # Pink
      "Large" = "#A7BA42",        # Green
      "Extra Large" = "springgreen",  # Light Green
      "Small_LDL" = "#1B9E77",         # Teal
      "Medium_LDL" = "#D95F02",        # Orange
      "Large_LDL" = "cornflowerblue" 
    )) +
    xlim(c(0, 5)) +
    facet_wrap(~Particle) +
    labs(
      title = paste("HDL Subclass Composition by Lifestyle:", lifestyle),
      fill = "SubClass"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 12)
    )
  # Save
  ggsave(
    filename = paste0("LDL_PieChart_", gsub(" ", "_", lifestyle), ".png"),
    plot = lifestyle_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
}
#########################################################################################################
## Test for statistical differences in Subclass abundance
# Filter particle count variables and combine groups
particles_data <- RAW_Data_for_NMR_ANALYSIS %>%
  select(MW_scores_lifestyle, Age, Sex, matches("_P$")) %>%
  mutate(MW_scores_lifestyle = case_when(
    MW_scores_lifestyle %in% c("Pastoralist", "PeriUrban") ~ "Non_Urban",
    MW_scores_lifestyle == "Urban" ~ "Urban",
    TRUE ~ NA_character_
  )) %>%
  pivot_longer(cols = matches("_P$"), names_to = "Subclass", values_to = "Value") %>%
  mutate(
    SubClass = case_when(
      str_starts(Subclass, "XL_") ~ "Extra Large",
      str_starts(Subclass, "L_LDL") ~ "Large_LDL",
      str_starts(Subclass, "M_LDL") ~ "Medium_LDL",
      str_starts(Subclass, "S_LDL") ~ "Small_LDL",
      str_starts(Subclass, "L_") ~ "Large",
      str_starts(Subclass, "M_") ~ "Medium",
      str_starts(Subclass, "S_") ~ "Small",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SubClass), !is.na(MW_scores_lifestyle)) %>%
  drop_na(Value, Age, Sex)
# Define custom subclass order (LDL first)
ordered_levels <- c("Large_LDL", "Medium_LDL", "Small_LDL",
                    "Extra Large", "Large", "Medium", "Small")
# Run ANOVA and generate facet labels
anova_results <- particles_data %>%
  mutate(SubClass = factor(SubClass, levels = ordered_levels)) %>%
  group_by(SubClass) %>%
  summarise(
    model = list(aov(Value ~ MW_scores_lifestyle + Age + Sex, data = cur_data())),
    .groups = "drop") %>%
  mutate(
    summary_stats = map(model, ~ summary(.x)[[1]]),
    F_value = map_dbl(summary_stats, ~ .x["MW_scores_lifestyle", "F value"]),
    p_value = map_dbl(summary_stats, ~ .x["MW_scores_lifestyle", "Pr(>F)"]),
    label_text = sprintf("%s\nF = %.2f, p = %.3g", SubClass, F_value, p_value))

# Merge labels into particles_data using a join
label_map <- anova_results %>%
  select(SubClass, label_text) %>%
  mutate(SubClass = factor(SubClass, levels = ordered_levels))

particles_data <- particles_data %>%
  mutate(SubClass = factor(SubClass, levels = ordered_levels)) %>%
  left_join(label_map, by = "SubClass") %>%
  mutate(label_text = factor(label_text, levels = unique(label_map$label_text)))

# Plot using facet labels
ggplot(particles_data, aes(x = MW_scores_lifestyle, y = Value, fill = MW_scores_lifestyle)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  facet_wrap(~label_text, scales = "free_y") +
  labs(
    title = "Comparison of Particle Counts by Lipoprotein Subclass (Urban vs Non-Urban)",
    x = "Lifestyle Group",
    y = "Number of Particles"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "plain"),
    legend.position = "none"
  )

##Correlation of Field Values and the NMR values
# Step 1: Define the column names for the molecules
## You can interchange LDL and HDL here 
## You can interchange LDL and HDL here 
molecule_columns <- c("S_HDL_C", "S_HDL_L", "S_HDL_P", "M_HDL_C", "M_HDL_L", "M_HDL_P",
                      "L_HDL_C", "L_HDL_L", "L_HDL_P", "XL_HDL_C", "XL_HDL_L", "XL_HDL_P",
                      "HDL_P", "HDL_L", "HDL_C")
# Scale all columns in the dataset
HDL_subset_data$HDL.cholesterol.mg.dL <- as.numeric(HDL_subset_data$HDL.cholesterol.mg.dL)
scaled_data <- HDL_subset_data %>%
  mutate(across(c("HDL.cholesterol.mg.dL", all_of(molecule_columns)), scale))
# Create empty list for plots
plots <- list()
# Loop and plot
for (molecule in molecule_columns) {
  # Determine subclass based on prefix
  subclass <- case_when(
    str_starts(molecule, "S_") ~ "Small",
    str_starts(molecule, "M_") ~ "Medium",
    str_starts(molecule, "L_") ~ "Large",
    str_starts(molecule, "XL_") ~ "Extra Large",
    TRUE ~ "Other")
  # Assign color based on subclass
  subclass_color <- case_when(
    subclass == "Small" ~ "#8B5FBF",      
    subclass == "Medium" ~ "plum1",       # Soft Green
    subclass == "Large" ~ "#A7BA42",        # Bright Pink
    subclass == "Extra Large" ~ "springgreen",  # Deep Purple
    TRUE ~ "#6B8E23")
  # Calculate correlations
  spearman_corr <- cor(scaled_data$HDL.cholesterol.mg.dL, scaled_data[[molecule]], method = "spearman", use = "complete.obs")
  pearson_corr <- cor(scaled_data$HDL.cholesterol.mg.dL, scaled_data[[molecule]], method = "pearson", use = "complete.obs")
  # Create plot
  p <- ggplot(scaled_data, aes_string(x = molecule, y = "HDL.cholesterol.mg.dL")) +
    geom_point(alpha = 0.2, color = subclass_color) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(
      title = paste("Scatter Plot of", molecule),
      x = molecule,
      y = "HDL (Field)") +
    annotate("text", x = -3, y = 4, label = paste0("Spearman: ", round(spearman_corr, 2)), hjust = 0, color = "red", size = 2) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 5),
      axis.title.x = element_text(size = 6),
      axis.title.y = element_text(size = 6),
      axis.text = element_text(size = 5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA))
  # Store the plot
  plots[[molecule]] <- p
}
# Step 5: Arrange and export the plots
png("HDL_Correlation_Plots.png", width = 1200, height = 1800, res = 300)
grid_plot <- ggarrange(plotlist = plots, ncol = 3, nrow = 5)
print(grid_plot)
dev.off()


##
## Working with Distribution of HDL subclasses with another Population
##I found a paper with summary stats from Hong Kong
#The provided data (mean ± SD for each subclass)
# Assign color based on subclass
# Define custom subclass color scheme
subclass_colors <- c(
  "Small" = "#8B5FBF",
  "Medium" = "plum1",        
  "Large" = "#A7BA42",         
  "Extra Large" = "springgreen")
# Input data
hong_kong_data <- data.frame(
  SubClass = c("Extra Large", "Large", "Medium", "Small"),
  Mean_Value = c(0.26, 1.4, 4.9, 15.0),          # Means for each subclass
  SD = c(0.11, 0.9, 1.2, 1.9))
# Calculate percentages and positions
hong_kong_data <- hong_kong_data %>%
  mutate(
    Percentage = Mean_Value / sum(Mean_Value) * 100,
    Label = paste0(round(Percentage, 1), "%"),
    ymax = cumsum(Percentage),
    ymin = c(0, head(ymax, n = -1)),
    mid = (ymax + ymin) / 2)
# Create pie chart with custom colors
pie_chart <- ggplot(hong_kong_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = SubClass)) +
  geom_rect(color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(x = 4.5, y = mid, label = Label), size = 5) +
  xlim(c(0, 5)) +
  scale_fill_manual(values = subclass_colors) +  # Apply custom colors
  labs(title = "HDL Subclass Distribution in the Hong Kong Population") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12))
# Save as PNG
ggsave("HDL_Subclass_Distribution_Hong_Kong_PieChart.png", plot = pie_chart, width = 8, height = 6, dpi = 300)

###### LDL subclasses https://pubmed.ncbi.nlm.nih.gov/30911231/
# Create a dataframe for both groups
ldl_data <- data.frame(
  SubClass = rep(c("LDL-1", "LDL-2", "LDL-3", "LDL-4"), 2),
  Mean_Value = c(34.6, 24.9, 5.5, 0.9, 33.8, 15.6, 2.3, 0.04),
  SD = c(10.6, 10.1, 6.7, 2.9, 8.1, 9.5, 3.3, 0.2),
  Group = rep(c("Patients", "Controls"), each = 4))
# Calculate total LDL per group
ldl_data <- ldl_data %>%
  group_by(Group) %>%
  mutate(Total = sum(Mean_Value),
         Percentage = (Mean_Value / Total) * 100,
         Label = paste0(round(Percentage, 1), "%"))
# Create pie chart plot
plot <- ggplot(ldl_data, aes(x = "", y = Percentage, fill = SubClass)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = Label), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4) +
  coord_polar(theta = "y") +
  facet_wrap(~ Group) +
  labs(title = "Percentage Distribution of LDL Subclasses by Group", x = NULL, y = NULL) + theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    strip.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12))

# Save the plot
ggsave("LDL_Subclass_Pie_Distribution_Patients_vs_Controls.png", plot = plot, width = 10, height = 6, dpi = 300)


################################################################################################
# Diet, Alcohol and Activity level consumption
## Clean columns needed for CVD risk calculation
## Clean the Cigarete and Diabetes responses into 1 or 0 to help in calculations of CVD risk
colnames(RAW_Data_for_NMR_ANALYSIS)
# Cigarretes 
CVD_Diet_study <- RAW_Data_for_NMR_ANALYSIS %>%
  mutate(Cigarettes.yes.no = case_when(
    Cigarettes.yes.no %in% c("Yes", "yes") ~ "1",
    Cigarettes.yes.no %in% c("No", "no") ~ "0",
    TRUE ~ NA_character_))
table(CVD_Diet_study$Cigarettes.yes.no)
# Diabetes
CVD_Diet_study <- CVD_Diet_study %>%
  mutate(Type.I.diabetes = case_when(
    Type.I.diabetes %in% c("Yes", "yes") ~ "1",
    Type.I.diabetes %in% c("No", "no", "Not Sure") ~ "0",
    TRUE ~ NA_character_))
table(CVD_Diet_study$Type.I.diabetes)
######################################
#######Additional cleaning of Files to make a clean file
#### Add a Lifestyle Category
CVD_Diet_study$lifestyle_grouping_benja <- NA
## Label as Rural based on conditions
table(CVD_Diet_study$MW_scores_h_sol_trim)
CVD_Diet_study[which(CVD_Diet_study$MW_scores_h_sol_trim <2),]$lifestyle_grouping_benja <- "Rural"
CVD_Diet_study$lifestyle_grouping_benja <- ifelse(CVD_Diet_study$MW_scores_h_sol_trim >= 2 & 
                                            (CVD_Diet_study$Main.subsistence.activity %in% c("Charcoal Burner", "Fisherman", "Farmer")),
                                          "Rural",
                                          CVD_Diet_study$lifestyle_grouping_benja)
# Label as "Market_Integrated" based on conditions
CVD_Diet_study$lifestyle_grouping_benja <- ifelse(CVD_Diet_study$MW_scores_h_sol_trim >= 2 & 
                                            !(CVD_Diet_study$Main.subsistence.activity %in% c("Charcoal Burner", "Fisherman", "Farmer")),
                                          "Market_Integrated",
                                          CVD_Diet_study$lifestyle_grouping_benja)
##Label as "Pastoralists"
CVD_Diet_study[which(CVD_Diet_study$MW_scores_h_sol_trim <3 & 
                       CVD_Diet_study$Main.subsistence.activity == "Animal Keeping"),]$lifestyle_grouping_benja <- "Pastoralist"

## Remove NA in lifestyle
CVD_Diet_study <-  CVD_Diet_study[complete.cases(CVD_Diet_study$lifestyle_grouping_benja), ]
## Diet and activity level - standardize responses 
## The Goal, compare Activity and Dietary differences between people
CVD_Diet_study <- CVD_Diet_study %>%
  mutate(
    ## Soda Consumption
    Soda.frequency.daily = case_when(
      Soda.frequency.daily %in% c("More Than Once Daily", "Once Daily", "Weekly", "Yes") ~ "2",
      Soda.frequency.daily == "Rarely" ~ "1",
      Soda.frequency.daily %in% c("No", "Never") | is.na(Soda.frequency.daily) ~ "0",
      TRUE ~ NA_character_
    ),
    ## Rice Consumption
    Rice.frequency.daily = case_when(
      Rice.frequency.daily %in% c("More Than Once Daily", "Daily", "Once Daily", "TWICE_weekly", "Yes") ~ "2",
      Rice.frequency.daily == "RARE" ~ "1",
      Rice.frequency.daily %in% c("No", "Never") | is.na(Rice.frequency.daily) ~ "0",
      TRUE ~ NA_character_
    ),
    ## Bread Consumption
    Bread.frequency = case_when(
      Bread.frequency %in% c(">2 Times Per Week", ">2 Times Per Week|Everyday", "1-2 Times Per Week|>2 Times Per Week", 
                             "Everyday", "1-2 Times Per Week", "Yes", "3") ~ "2",
      Bread.frequency == "Rarely" ~ "1",
      Bread.frequency %in% c("Never|Rarely", "Never") ~ "0",
      TRUE ~ NA_character_
    ),
    ## Fetch Water Frequency
    Fetch_water = case_when(
      Fetch_water == "Most_of_the_day" ~ "3",
      Fetch_water == "Few_hours" ~ "2",
      Fetch_water %in% c("Hours_1_and_below", "No") ~ "1",
      TRUE ~ NA_character_
    ),
    ## Chips
    Chips = case_when(
      Chips %in% c("Daily", "About_Twice_Weekly") ~ "2",
      Chips %in% c("RARE") ~ "1",
      Chips %in% c("NEVER", "NO", "No") ~ "0",
      TRUE ~ NA_character_
    ),
    ## Sweets
    Sweets = case_when(
      Sweets %in% c("YES", "Daily", "About_Twice_Weekly") ~ "2",
      Sweets == "RARE" ~ "1",
      Sweets %in% c("NEVER", "NO", "No") ~ "0",
      TRUE ~ NA_character_
    ),
    ## Chapati
    Chapati = case_when(
      Chapati %in% c("YES", "Daily", "About_Twice_Weekly") ~ "2",
      Chapati == "RARE" ~ "1",
      Chapati %in% c("NEVER", "NO", "No") ~ "0",
      TRUE ~ NA_character_))

table(Cleaned_metadata_file$Fetch_water)
table(Cleaned_metadata_file$Sweets)
table(Cleaned_metadata_file$Chapati)
table(Cleaned_metadata_file$Soda)
table(Cleaned_metadata_file$Bread_Frequency)
table(Cleaned_metadata_file$Rice)
table(Cleaned_metadata_file$Chips)
colnames(Cleaned_metadata_file)

##Refined Diet
Cleaned_metadata_file$Rice <- as.numeric(Cleaned_metadata_file$Rice)
Cleaned_metadata_file$Bread_Frequency <- as.numeric(Cleaned_metadata_file$Bread_Frequency)
Cleaned_metadata_file$Soda <- as.numeric(Cleaned_metadata_file$Soda)
Cleaned_metadata_file$Chips <- as.numeric(Cleaned_metadata_file$Chips)
Cleaned_metadata_file$Chapati <- as.numeric(Cleaned_metadata_file$Chapati)
Cleaned_metadata_file$Sweets <- as.numeric(Cleaned_metadata_file$Sweets)

Cleaned_metadata_file$Diet_MI_items <- rowSums(Cleaned_metadata_file[, c("Rice","Bread_Frequency","Soda","Chips","Chapati","Sweets")])

table(Cleaned_metadata_file$Diet_MI_items)

Cleaned_metadata_file <- transform(
  Cleaned_metadata_file, Alcohol = ifelse(Alcohol %in% c("Daily"), "2",
                                          ifelse(Alcohol %in% c("never", "No", "Never"), "0",
                                                 ifelse(Alcohol == "Occasionally", "1",
                                                        ifelse(Alcohol == "NA", "0", Alcohol)))))
table(Cleaned_metadata_file$Alcohol)





# Load and clean DNA and RNA
Dna_data <- import("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/DNA_2025_emma.csv") |>
  rename(UID_DNA = 3) |>
  select(UID_DNA) |>
  mutate(Unique.ID = UID_DNA)
colnames(Dna_data)
#
rna_data <- import("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/RNA_2025_EMMA.xlsx") |>
  rename(UID_RNA = 3) |>
  select(UID_RNA) |>
  mutate(Unique.ID = UID_RNA)
colnames(rna_data)
#
SC_RNA <- import("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/UIDs_passQC.txt")
names(SC_RNA)[1] <- "SC_RNA"
SC_RNA <- SC_RNA |>
  mutate(Unique.ID = SC_RNA)
colnames(SC_RNA)
# Subset and clean NMR
NMR_uniques <- select(Data_for_NMR_ANALYSIS, c(1, 3, 250, 59))  # adjust cols as needed
names(NMR_uniques)[1] <- "nmr_uniques"
NMR_uniques <- NMR_uniques |>
  filter(!is.na(LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml)) |>
  mutate(Unique.ID = nmr_uniques)
colnames(NMR_uniques)
# Stack all Unique.IDs from all sources
all_ids <- full_join(
  full_join(
    full_join(Dna_data, rna_data, by = "Unique.ID", relationship = "many-to-many"),
    SC_RNA, by = "Unique.ID", relationship = "many-to-many"
  ),
  NMR_uniques, by = "Unique.ID", relationship = "many-to-many")

head(all_ids)
#  replace "BLANK" or "" with NA
all_ids <- all_ids |>
  mutate(across(
    c(UID_DNA, UID_RNA, SC_RNA, nmr_uniques, Unique.ID),
    ~ ifelse(. %in% c("BLANK", ""), NA, .)))

dim(all_ids)
filtered_data <- all_ids |>
  filter(!(is.na(UID_DNA) & is.na(UID_RNA) & is.na(SC_RNA) & is.na(nmr_uniques) & is.na(Unique.ID)))
dim(filtered_data)
# Save the fully merged file
write.csv(filtered_data, "Master_Merged_All_Entries.csv", row.names = FALSE)

library(ggvenn)
colnames(filtered_data)
# Create sets of non-NA Unique IDs in each column
venn_input <- list(
  SC_RNA    = filtered_data$SC_RNA[!is.na(filtered_data$SC_RNA)],
  nmr_uniques = filtered_data$nmr_uniques[!is.na(filtered_data$nmr_uniques)],
  UID_DNA   = filtered_data$UID_DNA[!is.na(filtered_data$UID_DNA)],
  UID_RNA   = filtered_data$UID_RNA[!is.na(filtered_data$UID_RNA)]
)

# Make sure all entries are character for compatibility
venn_input <- lapply(venn_input, as.character)

# Plot Venn diagram
ggvenn(venn_input, fill_alpha = 0.5, show_percentage = FALSE)

ggVennDiagram(venn_input, 
              label_alpha = 0,
              label = "count",
              show_intersect = TRUE) +
  theme_void() +
  labs(title = "Unique ID Overlap Across Datasets") +
  theme(plot.title = element_text(hjust = 0.5))

# Ensure all UID vectors are character
venn_input <- list(
  SC_RNA      = as.character(filtered_data$SC_RNA[!is.na(filtered_data$SC_RNA)]),
  nmr_uniques = as.character(filtered_data$nmr_uniques[!is.na(filtered_data$nmr_uniques)]),
  UID_DNA     = as.character(filtered_data$UID_DNA[!is.na(filtered_data$UID_DNA)]),
  UID_RNA     = as.character(filtered_data$UID_RNA[!is.na(filtered_data$UID_RNA)])
)
# 1. Only in SC_RNA
only_SC_RNA <- setdiff(venn_input$SC_RNA, unlist(venn_input[c("nmr_uniques", "UID_DNA", "UID_RNA")]))
# 2. In both SC_RNA and nmr_uniques, but not in UID_DNA or UID_RNA
SC_RNA_and_nmr <- intersect(venn_input$SC_RNA, venn_input$nmr_uniques)
SC_RNA_nmr_only <- setdiff(SC_RNA_and_nmr, union(venn_input$UID_DNA, venn_input$UID_RNA))
# 3. Only in nmr_uniques
only_nmr <- setdiff(venn_input$nmr_uniques, unlist(venn_input[c("SC_RNA", "UID_DNA", "UID_RNA")]))
# Combine all into one data frame
combined_uids <- data.frame(
  UID = c(only_SC_RNA, SC_RNA_nmr_only, only_nmr),
  Category = c(
    rep("only_SC_RNA", length(only_SC_RNA)),
    rep("SC_RNA_nmr_only", length(SC_RNA_nmr_only)),
    rep("only_nmr", length(only_nmr))
  )
)
table(combined_uids$Category)
# Save to a single CSV
write.csv(combined_uids, "Combined_UID_Categories_for_EMMA.csv", row.names = FALSE)

# For Metabolomics
# Step 1–3: Same as before (categorizing SC_RNA UIDs)
uids_SC_RNA      <- unique(venn_input$SC_RNA)
uids_NMR         <- venn_input$nmr_uniques
uids_DNA         <- venn_input$UID_DNA

SC_RNA_DNA_NMR <- intersect(intersect(uids_SC_RNA, uids_DNA), uids_NMR)
SC_RNA_DNA     <- setdiff(intersect(uids_SC_RNA, uids_DNA), SC_RNA_DNA_NMR)
SC_RNA_NMR     <- setdiff(intersect(uids_SC_RNA, uids_NMR), SC_RNA_DNA_NMR)
SC_RNA_only    <- setdiff(uids_SC_RNA, union(union(SC_RNA_DNA_NMR, SC_RNA_DNA), SC_RNA_NMR))

SC_df <- data.frame(
  UID = c(SC_RNA_DNA_NMR, SC_RNA_DNA, SC_RNA_NMR, SC_RNA_only),
  Source = c(
    rep("SC_RNA ∩ DNA ∩ NMR", length(SC_RNA_DNA_NMR)),
    rep("SC_RNA ∩ DNA", length(SC_RNA_DNA)),
    rep("SC_RNA ∩ NMR", length(SC_RNA_NMR)),
    rep("SC_RNA only", length(SC_RNA_only))
  )
)

# Step 4: Remaining slots
n_needed <- 350 - nrow(SC_df)

# Step 5: Available DNA ∩ NMR not used in SC_df
uids_used <- SC_df$UID
uids_DNA_NMR <- setdiff(intersect(uids_DNA, uids_NMR), uids_used)

# Step 6: Adjust to number available
n_available <- length(uids_DNA_NMR)
n_to_add <- min(n_needed, n_available)

UID_DNA_NMR_df <- data.frame(
  UID = head(uids_DNA_NMR, n_to_add),
  Source = rep("UID_DNA ∩ NMR", n_to_add)
)
# Step 7: Combine and save
final_df <- rbind(SC_df, UID_DNA_NMR_df)
# Warning if not full
if (nrow(final_df) < 350) {
  warning(paste0("Only ", nrow(final_df), " UIDs available. Couldn't reach 350 due to insufficient UID_DNA ∩ NMR entries."))
}
# View breakdown
print(table(final_df$Source))
# Save to CSV
write.csv(final_df, "Selected_350_UIDs_detailed.csv", row.names = FALSE)

# Create a presence/absence table
presence_df <- filtered_data |>
  mutate(
    SC_RNA = !is.na(SC_RNA),
    NMR = !is.na(nmr_uniques),
    DNA = !is.na(UID_DNA),
    RNA = !is.na(UID_RNA)
  ) |>
  select(SC_RNA, NMR, DNA, RNA)

# Add dummy ID for plotting
presence_df$ID <- rownames(presence_df)

# Plot
upset(presence_df, intersect = c("SC_RNA", "NMR", "DNA", "RNA"), name = "Entry Count")



