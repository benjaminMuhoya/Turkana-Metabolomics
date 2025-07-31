#!/usr/bin/env python
import pandas as pd
import numpy as np
import statsmodels.api as sm
from tqdm import tqdm
import matplotlib.pyplot as plt
# === 1. Load & clean column names ===
df = pd.read_csv("HILIC_RAW_DATA_with_covariates.csv")
df = df.loc[:, df.columns.notna()]  # drop any columns whose name is NaN

# === 2. Load mapping & build rename dict ===
map_df = pd.read_csv("Updated_No_match_Metabo.csv")
rename_dict = {
    q: hmdb
    for q, hmdb in zip(map_df["Query"], map_df["HMDB"])
    if pd.notna(q) and pd.notna(hmdb)
}

# === 3. Rename columns & purge residual NaNs/duplicates ===
df = df.rename(columns=rename_dict)
df = df.loc[:, df.columns.notna()]       # drop any “nan” column names
df = df.loc[:, ~df.columns.duplicated()] # keep only the first of any duplicate names

# === 4. Filter rows & define grouping variable ===
df = df[df["MW_scores_lifestyle"].notna()].copy()
df["MW_scores_lifestyle_2"] = df["MW_scores_lifestyle"].apply(
    lambda x: "Urban" if x == "Urban" else "Non_Urban"
)

# === 5. Identify metadata vs. metabolite columns ===
meta_cols = [
    "Unique.ID", "indiv.ID", "Age", "Sex", "batch","Run.ID","MW_scores_market_diet_index","TA_score_TL_108",
    "MW_scores_lifestyle", "MW_scores_lifestyle_2"
]
potential_mets = [c for c in df.columns if c not in meta_cols]

# === 6. Coerce metabolite columns to numeric ===
for col in potential_mets:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# === 7. Outlier removal by IQR per metabolite ===
for col in potential_mets:
    s = df[col]
    q1, q3 = s.quantile(0.25), s.quantile(0.75)
    iqr = q3 - q1
    lower, upper = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    df.loc[(s < lower) | (s > upper), col] = np.nan

# === 8. Median‐center & log‐transform (normLogIC) ===
# 8a) Compute medians per metabolite
medians = df[potential_mets].median(skipna=True)
# 8b) Normalize by median
normIC = df[potential_mets].divide(medians, axis=1)
# 8c) Global small value to stabilize zeros
min_nonzero = normIC.replace(0, np.nan).abs().min().min()
min_val = min_nonzero / 10
# 8d) Compute log transform
normLogIC = np.log10((normIC + np.sqrt(normIC**2 + min_val**2)) / 2)
# 8e) Replace original metabolite values
df[potential_mets] = normLogIC

# === 9. Prepare for residualization ===
metabolite_cols = [
    col for col in potential_mets
    if df[col].notna().sum() > 0 and df[col].nunique() > 1
]
resid_df = df[["Unique.ID", "MW_scores_lifestyle_2"]].copy()
skipped, success = [], 0

# === 10. Residualize each metabolite ===
for met in tqdm(metabolite_cols, desc="Residualizing"):
    y = df[met].astype(float)
    X = df[["Age", "Sex", "batch"]].copy()
    X["Age"] = pd.to_numeric(X["Age"], errors="coerce")
    X["Sex"] = X["Sex"].astype("category")
    X["batch"] = X["batch"].astype("category")
    X = pd.get_dummies(X, drop_first=True)
    X = sm.add_constant(X).astype(float)

    try:
        model = sm.OLS(y, X, missing="drop").fit()
        resid_df[met] = model.resid
        success += 1
    except Exception as e:
        print(f"❌ Failed regression for {met}: {e}")
        skipped.append(met)

# === 11. Save outputs & print summary ===
resid_df.to_csv("HILIC_HMDB_Residualized_FINAL.csv", index=False)
pd.DataFrame(skipped, columns=["Skipped"]).to_csv("Skipped_Metabolites_FINAL.csv", index=False)

print("\n📊 Residualization Summary:")
print(f"Total metabolites: {len(metabolite_cols)}")
print(f"Matched HMDB IDs: {sum(m in rename_dict.values() for m in metabolite_cols)}")
print(f"✅ Successful regressions: {success}")
print(f"❌ Skipped regressions: {len(skipped)}")
# === 12. Compute & save missingness ===
missingness = df[metabolite_cols].isna().mean()
missingness.to_csv("metabolite_missingness.csv", header=["fraction_missing"])

# === 13. Plot missingness histogram with annotation ===
import matplotlib.pyplot as plt

total = len(missingness)
plt.figure(figsize=(8,5))
plt.hist(missingness * 100, bins=20, edgecolor="k")
plt.xlabel("Percent missing per metabolite")
plt.ylabel("Number of metabolites")
plt.title("Distribution of metabolite missingness")

# annotate total n on the plot
plt.text(
    0.98, 0.95,
    f"Total metabolites: {total}",
    transform=plt.gca().transAxes,
    ha="right", va="top",
    fontsize=12, fontweight="bold"
)

plt.tight_layout()
plt.savefig("missingness_histogram.png")
plt.show()

# === 14. Filter metabolites by missingness & report count ===
threshold = 0.10  # 10%
remaining = (missingness <= threshold).sum()
print(f"Metabolites with ≤{threshold*100:.0f}% missingness: {remaining} / {total}")

# === 15. (Optional) Bar chart of the worst offenders ===
plt.figure(figsize=(10,6))
missingness.sort_values(ascending=False).head(30).plot(kind="bar")
plt.ylabel("Fraction missing")
plt.title("Top 30 metabolites by missingness")
plt.tight_layout()
plt.savefig("missingness_bar_top30.png")
plt.show()

