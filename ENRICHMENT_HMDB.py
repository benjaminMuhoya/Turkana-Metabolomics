#!/usr/bin/env python3
"""
ENRICHMENT_HMDB.py  –  master runner (2025-07-31 patch 2)

Usage
-----
    python ENRICHMENT_HMDB.py                # classic run
    python ENRICHMENT_HMDB.py --use_ontology # fold in depth-4 subclasses

Pipeline
--------
0.  Rename_to_HMDB.py  →  HILIC_HMDB_Residualized_FINAL.csv
1.  Parse HMDB XML → build pathway + synonym map
2.  Match user-side metabolites to HMDB pathways
3.  (opt) ANNOTATE_THE_REST.py → merge depth-4 subclasses as “CHEM:<label>”
4.  Save user_matched_HMDB_pathways.tsv
5.  STATISTICS_ENRICHMENT.py → ANOVA, enrichment, figures
"""
# ── imports ──────────────────────────────────────────────────────
import argparse, subprocess, sys, pandas as pd, xml.etree.ElementTree as ET
from pathlib import Path

# ── CLI flag ─────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Residualize ➜ (opt) ontology ➜ stats")
parser.add_argument("--use_ontology", action="store_true",
                    help="Merge depth-4 subclasses from ANNOTATE_THE_REST.py")
args = parser.parse_args()

# ── 0. residualization ───────────────────────────────────────────
print("\n➤ Running Rename_to_HMDB.py …")
subprocess.run([sys.executable, "Rename_to_HMDB.py"], check=True)
print("✅ Residualized data ready.\n")

# ── 1. load residualized data ────────────────────────────────────
DATA_FILE = "HILIC_HMDB_Residualized_FINAL.csv"
XML_FILE  = "hmdb_metabolites_ALL_OF_THEM.xml"

df        = pd.read_csv(DATA_FILE)
hmdb_ids  = df.columns[2:]                       # first 2 cols = IDs + group

print(f"📊 {df.shape[0]} samples • {len(hmdb_ids)} metabolites")
print("\n👥 Lifestyle counts:")
print(df['MW_scores_lifestyle_2'].value_counts().to_string())

# ── 2. parse HMDB XML → pathway & synonym maps ──────────────────
print("\n🔍 Building pathway & synonym maps from HMDB XML …")
URI  = "{http://www.hmdb.ca}"                    # namespace URI
ctx  = ET.iterparse(XML_FILE, events=("end",))

hmdb_pathway_map, synonym_map = {}, {}
for _, elem in ctx:
    if elem.tag.endswith("metabolite"):
        acc = elem.findtext(f"{URI}accession")
        if not acc:
            elem.clear(); continue

        # pathways
        pws = [pw.findtext(f"{URI}name").strip()
               for pw in elem.findall(f".//{URI}pathway")
               if pw.findtext(f"{URI}name")]
        hmdb_pathway_map[acc] = list(dict.fromkeys(pws))

        # synonyms (incl. common name)
        common = elem.findtext(f"{URI}name")
        if common:
            synonym_map[common.lower()] = acc
        for syn in elem.findall(f".//{URI}synonym"):
            if syn.text:
                synonym_map[syn.text.strip().lower()] = acc

        elem.clear()

print(f"✅ {len(hmdb_pathway_map)} metabolites parsed")
print(f"✅ {len(synonym_map)} synonyms collected")

# metabolites lacking pathways
no_pw = [acc for acc, pws in hmdb_pathway_map.items() if not pws]
pd.Series(no_pw, name="HMDB_ID").to_csv("hmdb_no_pathways.tsv",
                                        sep="\t", index=False)
print(f"⚠️  {len(no_pw)} metabolites lack pathways (logged).")

# ── 3. match user metabolites to HMDB accessions ────────────────
matched, unmatched = [], []
for h in hmdb_ids:
    if h in hmdb_pathway_map:
        matched.append(h)
    else:
        acc = synonym_map.get(h.lower())
        (matched if acc in hmdb_pathway_map else unmatched).append(acc or h)

print(f"🔗 Matched {len(matched)} / {len(hmdb_ids)} metabolites")
print(f"❌ Unmatched: {len(unmatched)}")

out_df = pd.DataFrame(
    {"HMDB_ID": m,
     "Pathways": "; ".join(hmdb_pathway_map[m])}
    for m in matched
)
out_df.to_csv("user_matched_HMDB_pathways.tsv", sep="\t", index=False)
print("📝 Saved user_matched_HMDB_pathways.tsv")

# ── 4. optional ontology harvest ────────────────────────────────
if args.use_ontology and no_pw:
    print("\n➤ Running ANNOTATE_THE_REST.py …")
    subprocess.run([sys.executable, "ANNOTATE_THE_REST.py"], check=True)

    lvl4 = (pd.read_csv("hmdb_level4_mapping.tsv", sep="\t")
              .dropna(subset=["chem_subclass_lvl4"])
              .rename(columns={"chem_subclass_lvl4": "Pathways"}))
    lvl4["Pathways"] = "CHEM:" + lvl4["Pathways"]          # tag origin
    new_rows = lvl4[~lvl4.HMDB_ID.isin(out_df.HMDB_ID)]

    out_df = pd.concat([out_df, new_rows], ignore_index=True)
    out_df.to_csv("user_matched_HMDB_pathways.tsv", sep="\t", index=False)
    print(f"✅ Added {len(new_rows)} depth-4 subclass rows "
          f"(total {len(out_df)})")
else:
    print("ℹ️  Ontology step skipped.")

# ── 5. fallback classification (optional QC) ────────────────────
if no_pw:
    fb = {}
    for _, elem in ET.iterparse(XML_FILE, events=("end",)):
        if elem.tag.endswith("metabolite"):
            acc = elem.findtext(f"{URI}accession")
            if acc in no_pw:
                cls = [c.text for c in elem.findall(f".//{URI}classification") if c.text]
                fb[acc] = cls or ["<no classification>"]
            elem.clear()
    pd.DataFrame.from_dict(fb, orient="index", columns=["Classifications"])\
      .to_csv("hmdb_fallback_classification.tsv", sep="\t")
    print("🗂  Wrote hmdb_fallback_classification.tsv")

# ── 6. run statistics & enrichment ──────────────────────────────
print("\n➤ Running STATISTICS_ENRICHMENT.py …")
subprocess.run([sys.executable, "STATISTICS_ENRICHMENT.py"], check=True)

print("\n🎉 Pipeline complete – check pathway_enrichment_results.tsv & figures!")

