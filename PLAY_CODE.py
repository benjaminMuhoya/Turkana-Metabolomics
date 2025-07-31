#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt

# === 1. Load matched HMDB→pathway file ===
matched_file = "user_matched_HMDB_pathways.tsv"
um = pd.read_csv(matched_file, sep="\t")

# === 2. Quick summary ===
n_rows, n_cols = um.shape
print(f"🧪 Matched file: {n_rows} metabolites × {n_cols} cols")
print("📌 Columns:", list(um.columns))
print("\n🔍 Sample:")
print(um.head(5).to_string(index=False))

# === 3. Count empty annotations ===
empty = um["Pathways"].isna() | (um["Pathways"].str.strip() == "")
n_empty = empty.sum()
print(f"\n⚠️ {n_empty}/{n_rows} metabolites have no pathways")

# === 4. Explode into (metabolite, pathway) rows ===
exploded_all = (
    um.assign(Pathway=um["Pathways"].str.split("; "))
      .explode("Pathway")
      .dropna(subset=["Pathway"])
)
exploded_all["Pathway"] = exploded_all["Pathway"].str.strip()
exploded_all = exploded_all[exploded_all["Pathway"] != ""]

# === 5. Pathway counts BEFORE removing singletons ===
pathway_counts_all = exploded_all["Pathway"].value_counts()

# === 6. Bubble plot BEFORE removing singletons ===
dist_all = pathway_counts_all.value_counts().sort_index()
x_all = dist_all.index.astype(int)
y_all = dist_all.values
sizes_all = y_all * 40

plt.figure(figsize=(8, 5))
plt.scatter(x_all, y_all, s=sizes_all, alpha=0.6, edgecolors="w")
plt.axvline(x=5, color="gray", linestyle="--", label="Threshold k = 5")
plt.xlabel("Metabolites mapped to a pathway")
plt.ylabel("Number of pathways in database")
plt.title("Bubble Plot BEFORE Removing Singletons")
plt.legend()
plt.tight_layout()
plt.savefig("bubble_plot_before_removing_singletons.png")
plt.show()

# === 7. Remove singletons (k=1) ===
pathway_counts = pathway_counts_all[pathway_counts_all > 1]

# === 8. Distribution summary AFTER filtering ===
print("After removing singletons:")
print(f"  • Pathways remaining: {len(pathway_counts)}")
print(f"  • Min k: {pathway_counts.min()}")
print(f"  • Median k: {pathway_counts.median()}")
print(f"  • Mean k: {pathway_counts.mean():.2f}")
print(f"  • Max k: {pathway_counts.max()}")

summary_text = (
    f"⚠️ {n_empty}/{n_rows} metabolites have no pathways\n\n"
    f"After removing singletons:\n"
    f"• Pathways remaining: {len(pathway_counts)}\n"
    f"• Min k: {pathway_counts.min()}\n"
    f"• Median k: {pathway_counts.median():.1f}\n"
    f"• Mean k: {pathway_counts.mean():.2f}\n"
    f"• Max k: {pathway_counts.max()}"
)

# === 9. Bubble plot AFTER removing singletons ===
dist_filt = pathway_counts.value_counts().sort_index()
x_filt = dist_filt.index.astype(int)
y_filt = dist_filt.values
sizes_filt = y_filt * 40

plt.figure(figsize=(8, 5))
plt.scatter(x_filt, y_filt, s=sizes_filt, alpha=0.6, edgecolors="w")
plt.axvline(x=5, color="gray", linestyle="--", label="Threshold k = 5")
for xi, yi in zip(x_filt, y_filt):
    if xi >= 5:
        plt.text(xi, yi+0.3, str(int(yi)), ha="center", va="bottom")
plt.xlabel("Metabolites mapped to a pathway")
plt.ylabel("Number of pathways in database")
plt.title("Bubble Plot AFTER Removing Singletons")

# Add summary text box
plt.text(1.02, 0.5, summary_text, transform=plt.gca().transAxes,
         fontsize=9, verticalalignment='center',
         bbox=dict(boxstyle="round", alpha=0.05))

plt.legend()
plt.tight_layout()
plt.savefig("bubble_plot_after_removing_singletons.png")
plt.show()

