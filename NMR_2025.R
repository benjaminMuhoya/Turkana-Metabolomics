setwd("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/")
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
# Read All_data file (Meta_ONLY)
All_data <- import("THGP_database_TurkanaOnly_corrected_merged_2025-01-28.txt")
# Define the columns to process for cleaning and testing
##Subset data to only the needed columns
columns_to_keep <- c(1, 25, 31, 445, 489:492, 499:509, 601:604, 606:609, 617, 620, 621, 639, 644:658, 680:929, 115)
# Subset the data to include only the selected columns
Data_for_NMR_ANALYSIS <- All_data %>% dplyr::select(all_of(columns_to_keep))
# Check the column names of the subset
# Rename columns by removing "Nightingale_NMR_" prefix
colnames(Data_for_NMR_ANALYSIS) <- gsub("^Nightingale_NMR_", "", colnames(Data_for_NMR_ANALYSIS))
# Remove rows with NA in columns related to NMR - Use the NMR Unique ID
#Data_for_NMR_ANALYSIS <- Data_for_NMR_ANALYSIS %>%
  #filter(if_any("LBP_ELISA_Batch", ~ !is.na(.))) ## Use the NMR Unique ID to help subset to NMR_DATA only
dim(Data_for_NMR_ANALYSIS)
#####Remove those ending with _pct & .PCT i.e. Percentages. 
Data_for_NMR_ANALYSIS <- Data_for_NMR_ANALYSIS %>%
  dplyr::select(-matches("_pct$"),  -matches("(?i)\\.PCT$"))
######### 
######### Process the columns for analysis
#########
######### Process all the Numeric ones - scale and remove outliers
columns_to_process <- c(2,4:25,27,31,38,40,42:214)
#####
## Factor Columns
Data_for_NMR_ANALYSIS$Sex <- as.factor(Data_for_NMR_ANALYSIS$Sex)
Data_for_NMR_ANALYSIS$LBP_ELISA_Batch <- as.factor(Data_for_NMR_ANALYSIS$LBP_ELISA_Batch)
Data_for_NMR_ANALYSIS$MW_scores_lifestyle <- as.factor(Data_for_NMR_ANALYSIS$MW_scores_lifestyle)
#######
#########
# Convert columns to numeric
convert_to_numeric <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    df[[column_name]] <- as.numeric(as.character(df[[column_name]]))
  }
  return(df)
}
Data_for_NMR_ANALYSIS <- convert_to_numeric(Data_for_NMR_ANALYSIS, columns_to_process)
# Standardize and normalize columns
standardize_and_normalize <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    df[[column_name]] <- scale(df[[column_name]], center = TRUE, scale = TRUE)
  }
  return(df)
}
Data_for_NMR_ANALYSIS <- standardize_and_normalize(Data_for_NMR_ANALYSIS, columns_to_process)
colnames(Data_for_NMR_ANALYSIS)
dim(Data_for_NMR_ANALYSIS)
#######
#########
# Function to remove outliers and replace with NA
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
dim(Data_for_NMR_ANALYSIS)

# Ensure lifestyle_group is properly classified
Data_for_NMR_ANALYSIS <- Data_for_NMR_ANALYSIS %>%
  mutate(lifestyle_group = ifelse(MW_scores_lifestyle %in% c("Pastoralist", "PeriUrban"),
                                  "Non_Urban", "Urban"))

# Define metabolite columns
colnames(Data_for_NMR_ANALYSIS)
dim(Data_for_NMR_ANALYSIS)
metabolite_cols <- names(Data_for_NMR_ANALYSIS[, c(5:25, 27, 42:214)])
Data_for_NMR_ANALYSIS[metabolite_cols] <- lapply(Data_for_NMR_ANALYSIS[metabolite_cols], as.numeric)

# Initialize dataframes to store results
all_coefficients <- data.frame()
# Loop through each metabolite and run the mixed-effects model
for (metabolite in metabolite_cols) {
  # Filter and select relevant columns
  current_data <- Data_for_NMR_ANALYSIS %>%
    dplyr::select(lifestyle_group, Age, Sex, MW_scores_market_diet_index_trim, LBP_ELISA_Batch, all_of(metabolite)) %>%
    rename(Metabolite = all_of(metabolite)) %>%
    drop_na(Metabolite, Age, Sex)
  
  # Ensure enough data remains to fit the model
  if (nrow(current_data) > 6) {
    # Convert batch to a factor
    current_data <- current_data %>%
      mutate(LBP_ELISA_Batch = as.factor(LBP_ELISA_Batch))
    
    # Fit the mixed-effects model with batch as a random effect
    model <- lmer(Metabolite ~ lifestyle_group + Sex + Age + MW_scores_market_diet_index_trim + 
                    lifestyle_group:Sex + (1 | LBP_ELISA_Batch), 
                  data = current_data, REML = FALSE)
    
    # Extract fixed effect coefficients, including p-values
    temp_results <- broom.mixed::tidy(model, effects = "fixed") %>%
      mutate(
        Metabolite = metabolite,
        N = nrow(current_data)  # Add sample size
      )
    
    # Append results to master dataframe
    all_coefficients <- bind_rows(all_coefficients, temp_results)
  } else {
    message(paste("Skipping metabolite:", metabolite, "- insufficient data"))
  }
}
head(all_coefficients)
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
write_csv(all_coefficients, "NMR_all_coefficients_results_lme.csv")

# Compute Fold Change from Regression Results
fold_change_results <- all_coefficients %>%
  filter(term == "lifestyle_groupUrban") %>%  # Only take the coefficient for Urban effect
  mutate(
    Log2FC = estimate,  # Since data is scaled, the coefficient directly represents Log2 Fold Change
    Adjusted_PValue = p.value  # Use Bonferroni-corrected p-values ##I am using raw, p-value here, with adjusted it changed
  ) %>%
  select(Metabolite, Log2FC, Adjusted_PValue)

# Classify significance
fold_change_results <- fold_change_results %>%
  mutate(
    Significant = case_when(
      Adjusted_PValue < 0.03 & Log2FC > 0.5 ~ "Higher in Urban",
      Adjusted_PValue < 0.03 & Log2FC < -0.3 ~ "Lower in Urban",
      TRUE ~ "Not Significant"
    )
  )
# Generate the volcano plot
volcano_plot <- ggplot(fold_change_results, aes(x = Log2FC, y = -log10(Adjusted_PValue), color = Significant)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    name = "Regulation",
    values = c(
      "Higher in Urban" = "blue",
      "Lower in Urban" = "red",
      "Not Significant" = "grey"
    )
  ) +
  geom_hline(yintercept = -log10(0.03), linetype = "dashed", color = "black") +
  geom_text_repel(data = fold_change_results %>% filter(Significant != "Not Significant"),
                  aes(label = Metabolite),
                  size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Urban vs. Non-Urban (LME Model) Raw p-value",
    x = "Log2 Fold Change (Urban / Non-Urban)",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Save the volcano plot
ggsave("volcano_plot.png", plot = volcano_plot, device = "png", width = 8, height = 6, dpi = 300)

# Display the plot
print(volcano_plot)

# Compute confidence intervals and exclude intercept
coeff_results <- all_coefficients %>%
  filter(term != "(Intercept)")

# Classify predictors for color coding
coeff_results <- coeff_results %>%
  mutate(Predictor = case_when(
    grepl("lifestyle_group", term) ~ "Lifestyle",
    grepl("Age", term) ~ "Age",
    grepl("Sex", term) ~ "Sex",
    grepl("MW_scores_market_diet_index_trim", term) ~ "Market Diet Index",
    TRUE ~ "Other"
  ))

# Filter for significant and strong effects (|β| ≥ 0.3 and CI excludes 0)
filtered_coeff_results <- all_coefficients %>%
  filter(abs(estimate) >= 0.35 & (lower_CI > 0 | upper_CI < 0) & term != "(Intercept)") %>%
  arrange(desc(abs(estimate)))  # Sort for better visualization

# Generate the coefficient plot
coef_plot <- ggplot(filtered_coeff_results, aes(x = Metabolite, y = estimate, ymin = lower_CI, ymax = upper_CI, color = term)) +
  geom_pointrange(alpha = 0.7, size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical reference line at 0
  geom_text(aes(label = N), hjust = -0.3, size = 3.5, color = "black") +  # Add sample size (N) next to points
  theme_minimal() +
  scale_color_manual(values = c("lifestyle_groupUrban" = "blue", 
                                "SexMale" = "red", 
                                "Age" = "green", 
                                "MW_scores_market_diet_index_trim" = "purple",
                                "lifestyle_groupUrban:SexMale" = "orange")) +
  labs(
    title = "Coefficient Plot: Effect of Predictors. LME model",
    x = "Metabolites",
    y = "Beta Estimate",
    color = "Predictor",
    caption = "Filtered for |β| ≥ 0.35 & Significant 95% CI"
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
  coord_flip()

# Save the coefficient plot
ggsave("coefficient_plot.png", plot = coef_plot, device = "png", width = 8, height = 6, dpi = 300)

# Display the plot
print(coef_plot)


colnames(Data_for_NMR_ANALYSIS)
dim(Data_for_NMR_ANALYSIS)
#######
#######
#########
## Define the columns needed for the regression
# Rename specific inflammation columns
HDL_subset_data <- Data_for_NMR_ANALYSIS %>%
  rename(
    CRP = `Colombia_C_reactive_protein..mg.L.`,
    Globulin = `Colombia_globulin..g.dL.`,
    Albumin = `Colombia_albumin__serum..g.dL.`,
    LBP = `LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml`,
    Albumin_NMR = `Albumin`
  )
# Define inflammation biomarkers
inflammation_cols <- c("CRP", "Globulin", "Albumin", "LBP", "GlycA", "Albumin_NMR")
# Define metabolite columns based on column indices
metabolite_indices <- c(5:14, 16:21, 24, 25, 42:113, 116:214) 
metabolite_cols <- colnames(HDL_subset_data)[metabolite_indices]

## LME Regression with Inflammation phenotypes
# Initialize a list to store results
all_results <- list()
# Loop through inflammation biomarkers
for (inflam in inflammation_cols) {
  cat("=== Running mixed-effects models for inflammation variable:", inflam, "===\n")
  # Data frame to store results
  beta_results <- data.frame(Metabolite = character(),
                             Predictor = character(),
                             Beta = numeric(),
                             p_value = numeric(),
                             N = numeric(),
                             stringsAsFactors = FALSE)
  for (metab in metabolite_cols) {
    formula_str <- paste(metab, "~", inflam, "+ lifestyle_group + lifestyle_group:Sex + Age + Sex + (1 | LBP_ELISA_Batch)")
    cat("Running model:", formula_str, "\n")
    
    # Subset data for relevant columns and remove NA rows
    sub_data <- HDL_subset_data %>%
      dplyr::select(all_of(c(metab, inflam, "lifestyle_group", "Age", "Sex", "LBP_ELISA_Batch"))) %>%
      na.omit()
    if (nrow(sub_data) == 0) {
      cat("  No data for", metab, "with inflammation", inflam, "\n")
      next
    }
    # Record sample size
    n_used <- nrow(sub_data)
    if (n_used < 10) {
      cat("  Too few data points for", metab, "with inflammation", inflam, "\n")
      next
    }
    # Fit the mixed-effects model
    model <- tryCatch(
      lmer(as.formula(formula_str), data = sub_data, REML = FALSE),
      error = function(e) {
        cat("  Skipping model for", metab, "with", inflam, "due to error:", e$message, "\n")
        return(NULL)
      }
    )
    if (is.null(model)) next
    # Extract the coefficients for all fixed effects using broom.mixed::tidy
    tidy_model <- broom.mixed::tidy(model, effects = "fixed")
    # Keep only relevant predictors
    fixed_effects <- tidy_model %>%
      filter(term %in% c(inflam, "lifestyle_group", "lifestyle_group:Sex", "Age", "Sex")) %>%
      select(term, estimate, p.value)
    # Append results
    if (nrow(fixed_effects) > 0) {
      fixed_effects <- fixed_effects %>%
        mutate(Metabolite = metab, N = n_used) %>%
        rename(Predictor = term, Beta = estimate, p_value = p.value)
      
      beta_results <- rbind(beta_results, fixed_effects)
    }
  }
  
  # Adjust p-values for multiple testing using Bonferroni correction
  if (nrow(beta_results) > 0) {
    beta_results$adj_p_value <- p.adjust(beta_results$p_value, method = "bonferroni")
  } else {
    beta_results$adj_p_value <- numeric(0)
  }
  
  # Keep only significant associations (adjusted p < 0.05)
  significant_results <- beta_results %>% filter(adj_p_value < 0.05)
  # Store results
  all_results[[inflam]] <- significant_results
  # Display summary of significant associations
  if (nrow(significant_results) > 0) {
    cat("Significant associations found for", inflam, "\n")
  } else {
    cat("No significant associations found for", inflam, "\n")
  }
}

head(significant_results)
# Define predictors to plot
predictors_to_plot <- c("GlycA", "Albumin_NMR", "Albumin")

# Define the number of top metabolites to display
top_n_metabolites <- 20  # Change this number to adjust how many metabolites to plot

# Define a function to create coefficient plots
plot_coefficients <- function(results_df, predictor_name) {
  if (nrow(results_df) == 0) {
    cat("No significant results for", predictor_name, "- skipping plot.\n")
    return(NULL)
  }
  
  # Calculate confidence intervals
  results_df <- results_df %>%
    mutate(
      CI_lower = Beta - 1.96 * (abs(Beta) / sqrt(N)),  # Approximate standard error
      CI_upper = Beta + 1.96 * (abs(Beta) / sqrt(N))
    ) %>%
    # Filter to include only significant results: Beta > 0.3 and CI excluding zero
    filter(abs(Beta) > 0.3, CI_lower > 0 | CI_upper < 0) %>%
    # Select top N metabolites sorted by significance (p-value)
    arrange(p_value) %>% 
    head(top_n_metabolites)
  
  if (nrow(results_df) == 0) {
    cat("No significant filtered results for", predictor_name, "- skipping plot.\n")
    return(NULL)
  }
  
  # Assign colors to different predictors
  predictor_colors <- c("GlycA" = "grey", "Albumin_NMR" = "grey", "Albumin" = "grey",
                        "lifestyle_group" = "purple", "lifestyle_group:Sex" = "yellow",
                        "Age" = "orange", "Sex" = "green")
  
  # Create the plot
  p <- ggplot(results_df, aes(x = Beta, y = reorder(Metabolite, Beta), color = Predictor)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
    geom_text(aes(label = N), vjust = -0.5, hjust = -0.2, size = 4) +
    scale_color_manual(values = predictor_colors) +
    labs(title = paste("Effect of", predictor_name, "on Metabolites"),
         x = "Effect Size (Beta)", y = "Metabolites",
         color = "Predictor") +
    theme_minimal()
  
  # Save plot
  plot_path <- file.path(paste0(predictor_name, "_coefficients.png"))
  ggsave(plot_path, p, width = 8, height = 6, dpi = 300)
  cat("Plot saved for", predictor_name, "at", plot_path, "\n")
}

# Loop through selected predictors and generate plots
for (predictor in predictors_to_plot) {
  if (predictor %in% names(all_results)) {
    plot_coefficients(all_results[[predictor]], predictor)
  } else {
    cat("Predictor", predictor, "not found in results - skipping.\n")
  }
}








