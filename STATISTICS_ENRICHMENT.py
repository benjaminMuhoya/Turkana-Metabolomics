#!/usr/bin/env python
#   STATISTICS_ENRICHMENT.py  ◇  2025-07-30  (fully-patched)
#   --------------------------------------------------------------
#   Feature testing → Fisher enrichment → dual-panel plot
#   → tidy box-plots for best pathways.
#   --------------------------------------------------------------

import subprocess, sys, textwrap, warnings
subprocess.run(["python", "PLAY_CODE.py"], check=True)

import pandas as pd, numpy as np, seaborn as sns, matplotlib.pyplot as plt
from scipy.stats import f_oneway, ttest_ind, fisher_exact
from statsmodels.stats.multitest import fdrcorrection
from adjustText import adjust_text
import warnings
warnings.filterwarnings(
    "ignore",
    message="Looks like you are using a tranform that doesn't support FancyArrowPatch"
)
warnings.filterwarnings(
    "ignore",
    message="This figure includes Axes that are not compatible with tight_layout"
)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

# ════════════════ CONFIG ════════════════
DATA_FILE      = "HILIC_HMDB_Residualized_FINAL.csv"
MAPPING_FILE   = "user_matched_HMDB_pathways.tsv"
ID2NAME_FILE   = "hmdb_id_to_name.tsv"       # optional look-up

OUT_ENRICH_TSV = "pathway_enrichment_results.tsv"
RAW_P_CUT      = 0.05   # raw-p threshold  ( left panel – dashed black)
FDR_Q_CUT      = 0.05   # FDR  threshold   (right panel – solid  red)
K_THRESHOLD    = 4      # min metabolites per pathway
# ════════════════════════════════════════
saved_figs = []          # <-- collect file-names here

# ───────────────────────────────── STEP 1 ─────────────────────────────────
df   = pd.read_csv(DATA_FILE)
grp  = "MW_scores_lifestyle_2"
if grp not in df.columns:
    sys.exit(f"✘  column “{grp}” not found.")

metabolites = df.columns[2:]
groups      = df[grp].dropna().unique()
print(f"⚙️  Detected {df.shape[0]} samples  –   groups: {groups}")

# ───────────────────────────────── STEP 2 ─────────────────────────────────
pvals = []
for met in metabolites:
    samples = [df.loc[df[grp]==g, met].dropna() for g in groups]
    if all(len(s) > 1 for s in samples):
        p = ttest_ind(*samples, nan_policy="omit")[1] if len(groups)==2 else f_oneway(*samples)[1]
    else:
        p = 1.0
    pvals.append(p)

_, qvals = fdrcorrection(pvals, alpha=FDR_Q_CUT)
stat_df  = pd.DataFrame({"HMDB_ID": metabolites, "pval": pvals, "qval": qvals})

sig_raw  = set(stat_df.loc[stat_df.pval < RAW_P_CUT, "HMDB_ID"])
sig_fdr  = set(stat_df.loc[stat_df.qval < FDR_Q_CUT, "HMDB_ID"])
sig_enr  = sig_fdr if sig_fdr else sig_raw      # for Fisher

print(f"✅ Features p<{RAW_P_CUT}: {len(sig_raw)} / {len(metabolites)}")
print(f"✅ Features q<{FDR_Q_CUT}: {len(sig_fdr)} / {len(metabolites)}")

# ───────────────────────────────── STEP 3 ─────────────────────────────────
map_df = pd.read_csv(MAPPING_FILE, sep="\t").dropna(subset=["Pathways"])
mp = map_df.set_index("HMDB_ID")["Pathways"].str.split("; ").to_dict()

path2mets = {}
for met, pws in mp.items():
    for pw in pws:
        path2mets.setdefault(pw, []).append(met)

pairs = [(pw, m) for pw, mets in path2mets.items() for m in mets]
all_df = pd.DataFrame(pairs, columns=["Pathway", "HMDB_ID"])
total_counts = all_df.Pathway.value_counts()

annotated     = set(mp.keys())
sig_annotated = sig_enr & annotated
sig_counts    = all_df[all_df.HMDB_ID.isin(sig_annotated)].Pathway.value_counts()

print("\n📈 Pathway-size (k) distribution:")
print(f"   min={total_counts.min()}  median={total_counts.median():.1f} "
      f"mean={total_counts.mean():.1f}  max={total_counts.max()}")
print(f"🎯 Keeping pathways with k ≥ {K_THRESHOLD}")

filtered = total_counts[total_counts >= K_THRESHOLD]

# ───────────────────────────────── STEP 4 ─────────────────────────────────
M, N = len(annotated), len(sig_annotated)
results = []
for pw, K in filtered.items():
    k = sig_counts.get(pw, 0)
    if k == 0:           # nothing significant in this pathway
        continue
    table = [[k, N-k], [K-k, (M-K)-(N-k)]]
    _, pval = fisher_exact(table, alternative="greater")

    dirs = pd.Series(path2mets[pw]).map(
        lambda m: df.groupby(grp)[m].mean().idxmax()).value_counts().to_dict()

    results.append(dict(Pathway=pw, Hits=k, Total=K, pval=pval,
                        Hits_Urban=dirs.get("Urban",0),
                        Hits_Non_Urban=dirs.get("Non_Urban",0),
                        Dominant_Group=("Urban" if dirs.get("Urban",0)>dirs.get("Non_Urban",0)
                                        else "Non_Urban" if dirs.get("Non_Urban",0)>dirs.get("Urban",0)
                                        else "Tie")))
enr_df = pd.DataFrame(results)
if enr_df.empty:
    sys.exit("🚫  No enriched pathways after filtering.")

enr_df["qval"] = fdrcorrection(enr_df.pval)[1]
enr_df.sort_values("pval", inplace=True)
enr_df.to_csv(OUT_ENRICH_TSV, sep="\t", index=False)

print("\n📌 Top enriched pathways:")
print(enr_df.head(10).to_string(index=False))

# ───────────────────────────────── STEP 5 ─────────────────────────────────
# ── 5.  VISUALISATION : dual-panel (pretty, bug-fixed) ──────────────────
enr_df["neglog10_p"] = -np.log10(np.clip(enr_df.pval, 1e-300, None))
enr_df["neglog10_q"] = -np.log10(np.clip(enr_df.qval, 1e-300, None))

fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True,
                         gridspec_kw=dict(wspace=0.12, width_ratios=[1, 1]))
fig.suptitle("Pathway enrichment", weight='bold', fontsize=16)

# palette used for both panels
dom_palette = {"Urban": "#1f77b4",
               "Non_Urban": "#ff7f0e",
               "Tie": "#2ca02c"}

def beautify_panel(ax, yvar, cutoff, color, ls, big_label, title):
    # build a plotting frame with a slight x-jitter
    plot_df = enr_df.copy()
    plot_df["Hits_jit"] = plot_df["Hits"] + np.random.uniform(-0.15, 0.15,
                                                              size=len(plot_df))

    # scatter
    sns.scatterplot(data=plot_df,
                    x="Hits_jit",
                    y=yvar,
                    hue="Dominant_Group",
                    palette=dom_palette,
                    s=90, edgecolor="w", linewidth=0.3,
                    ax=ax, legend=False)

    thr = -np.log10(cutoff)
    ax.axhline(thr, color=color, ls=ls, lw=2 if big_label else 1.4)

    ax.text(0.01, 0.97, f"{'q' if yvar=='neglog10_q' else 'p'} = {cutoff}",
            transform=ax.transAxes,
            color=color,
            fontweight='bold',
            fontsize=18 if big_label else 11,
            va='top', ha='left')

    # label only points above the threshold
    label_df = plot_df[plot_df[yvar] >= thr]
    txts = [ax.text(r.Hits_jit, r[yvar], r.Pathway,
                    fontsize=8, weight='bold', ha="left", va="center",
                    clip_on=True)
            for _, r in label_df.iterrows()]

    adjust_text(txts, ax=ax,
                only_move={'texts': 'y'},
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.4),
                expand_text=(1.2, 1.3), expand_points=(1.15, 1.25),
                force_text=(0.6, 0.6), force_points=(0.3, 0.3))

    ax.set_title(title, fontsize=13, pad=10)
    ax.set_xlabel("Hits per pathway", fontsize=11)
    if ax is axes[0]:
        ax.set_ylabel("-log10(p-value)", fontsize=11)
    ax.spines[['right', 'top']].set_visible(False)

# left: raw-p
beautify_panel(axes[0], "neglog10_p", RAW_P_CUT,
               color="black", ls="--", big_label=False,
               title="Raw-p threshold")

# right: FDR-q
beautify_panel(axes[1], "neglog10_q", FDR_Q_CUT,
               color="red", ls="-", big_label=True,
               title="FDR q threshold")

# legend once, on the right
handles, labels = axes[0].get_legend_handles_labels()
axes[1].legend(handles, labels, title="Dominant Group",
               frameon=True, loc="lower right")

plt.tight_layout(rect=[0, 0, 1, 0.95])
# dual-panel plot
panel_fname = "enrichment_dual_panel.png"
plt.savefig(panel_fname, dpi=300)
plt.close()
saved_figs.append(panel_fname)          # <-- add here

# ───────────────────────────────── STEP 6 ─────────────────────────────────
try:
    id2name = (pd.read_csv(ID2NAME_FILE, sep="\t", header=None,
                names=["HMDB_ID","Name"])
               .set_index("HMDB_ID")["Name"].to_dict())
except FileNotFoundError:
    id2name = {}

# ───────────────────────────────── STEP 7 ─────────────────────────────────
path_sel = enr_df.loc[enr_df.qval < FDR_Q_CUT, "Pathway"]
if path_sel.empty:
    path_sel = enr_df.loc[enr_df.pval < RAW_P_CUT, "Pathway"]

if path_sel.empty:
    print("⚠️  No pathways pass either threshold – skipping box-plots.")
    sys.exit(0)

df[grp] = df[grp].astype(str)
palette = dict(zip(sorted(df[grp].unique()),
                   sns.color_palette("deep", len(df[grp].unique()))))

for pw in path_sel:
    mets = [m for m in path2mets[pw] if m in sig_raw][:5]
    if not mets: continue

    melt_df = (df[[grp] + mets]
               .melt(id_vars=grp, var_name="HMDB_ID", value_name="Value"))
    melt_df["Label"] = melt_df.HMDB_ID.map(id2name).fillna(melt_df.HMDB_ID)

    plt.figure(figsize=(max(6, 1.7*len(mets)), 6))
    ax = sns.boxplot(data=melt_df, x="Label", y="Value",
                     hue=grp, palette=palette,
                     showcaps=True, fliersize=0, linewidth=1)
    sns.stripplot(data=melt_df, x="Label", y="Value", hue=grp,
                  dodge=True, palette=palette, edgecolor="gray",
                  linewidth=0.4, alpha=0.55, size=3)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),
              title=grp, bbox_to_anchor=(1.02, 1), loc="upper left")

    ax.set_title(f"{pw}  –  top {len(mets)} raw-p metabolites")
    ax.set_xlabel("") ; ax.set_ylabel("Abundance")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    safe_pw = pw.replace(" ", "_").replace("/", "_").replace(":", "")
    fname   = f"boxplot_{safe_pw}.png"
    plt.savefig(fname, dpi=300)
    plt.close()
    saved_figs.append(fname)                # <-- add here

print("\n✅  Generated figures:")
for f in saved_figs:
    print("   •", f)

