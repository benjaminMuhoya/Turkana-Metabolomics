#!/usr/bin/env python3
# ── ontology_probe.py  ─────────────────────────────────────────────
"""
Harvest extra ontology information for HMDB accessions that have *no* pathway
annotation, then dump a complete ‘level-by-level’ taxonomy count table and
a plot of the most frequent depth-4 subclasses.

Outputs
-------
1. hmdb_no_pathways_ontology.tsv          – detailed per-metabolite info
2. hmdb_no_pathways_ontology_preview.txt  – quick human-readable preview
3. hmdb_taxonomy_label_counts.tsv         – every taxonomy label at every depth
4. hmdb_level4_mapping.tsv                – HMDB_ID → depth-4 subclass
5. taxonomy_depth4_top20.png              – bar-chart of most common subclasses
"""

import xml.etree.ElementTree as ET
from   collections import defaultdict, Counter
from   pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt   # << new for the figure

XML_FILE  = "hmdb_metabolites_ALL_OF_THEM.xml"
NO_PW_TSV = "hmdb_no_pathways.tsv"
PREVIEW_N = 5     # schema peek
TOP_N     = 40    # how many depth-4 labels to plot

# ──────────────────────────────────────────────────────────────────
# 1. accession list (those without pathways)
# ──────────────────────────────────────────────────────────────────
no_pw = set(pd.read_csv(NO_PW_TSV, sep="\t")["HMDB_ID"])
print(f"🔍 Probing ontology for {len(no_pw):,} metabolites without pathways…")

ns       = {"hmdb": "http://www.hmdb.ca"}
records  = []
explored = defaultdict(set)

def grab(elem, path):
    """Return sorted unique texts for any XPath inside *elem*"""
    return sorted({t.text.strip() for t in elem.findall(path, ns)
                   if t is not None and t.text})

# ──────────────────────────────────────────────────────────────────
# 2. main streaming parse
# ──────────────────────────────────────────────────────────────────
peek = 0
for _, met in ET.iterparse(XML_FILE, events=("end",)):
    if not met.tag.endswith("metabolite"):
        continue
    acc = met.findtext("hmdb:accession", namespaces=ns)
    if acc not in no_pw:
        met.clear()
        continue                      # skip metabolites that *do* have a pathway

    # --- schema peek (first PREVIEW_N) ----------------------------
    if peek < PREVIEW_N:
        onto = met.find("hmdb:ontology", ns)
        if onto is not None:
            for ch in onto.iter():
                explored[ch.tag.split('}')[-1]].add(ch.tag)
        peek += 1

    # --- real harvest --------------------------------------------
    rec = {"HMDB_ID": acc}

    onto = met.find("hmdb:ontology", ns)
    if onto is not None:
        # origins / biofunctions / applications
        rec["origin"]       = "; ".join(grab(onto, ".//hmdb:origins/hmdb:origin"))
        rec["biofunction"]  = "; ".join(grab(onto, ".//hmdb:biofunctions/hmdb:biofunction"))
        rec["application"]  = "; ".join(grab(onto, ".//hmdb:applications/hmdb:application"))

        # disposition
        disp = onto.find("hmdb:disposition", ns)
        if disp is not None:
            rec["biofluid_locations"] = "; ".join(grab(disp, ".//hmdb:biofluid_locations/hmdb:biofluid"))
            rec["tissue_locations"]   = "; ".join(grab(disp, ".//hmdb:tissue_locations/hmdb:tissue"))
            rec["cellular_locations"] = "; ".join(grab(disp, ".//hmdb:cellular_locations/hmdb:cellular_location"))
            rec["organ"]              = "; ".join(grab(disp, ".//hmdb:organs/hmdb:organ"))
            rec["route_of_exposure"]  = "; ".join(grab(disp, ".//hmdb:route_of_exposure/hmdb:route_of_exposure"))
            rec["source"]             = "; ".join(grab(disp, ".//hmdb:sources/hmdb:source"))

        # processes / roles
        rec["process"] = "; ".join(grab(onto, ".//hmdb:processes/hmdb:process"))
        rec["role"]    = "; ".join(grab(onto, ".//hmdb:roles/hmdb:role"))

    # taxonomy
    tax = met.find("hmdb:taxonomy", ns)
    if tax is not None:
        pieces = [tax.findtext(f"hmdb:{tag}", namespaces=ns) for tag
                  in ("kingdom", "super_class", "class", "sub_class")]
        rec["taxonomy"] = " > ".join(p for p in pieces if p)

    records.append(rec)
    met.clear()

# ──────────────────────────────────────────────────────────────────
# 3. schema preview
# ──────────────────────────────────────────────────────────────────
print("\n🗂  First few <ontology> tag names detected:")
for t in sorted(explored)[:12]:
    print(f"   • {t}")

# ──────────────────────────────────────────────────────────────────
# 4. tidy per-metabolite TSV & preview
# ──────────────────────────────────────────────────────────────────
onto_df = pd.DataFrame(records).sort_values("HMDB_ID")
onto_df["tax_list"]  = onto_df["taxonomy"].fillna("").str.split(" > ")
onto_df["tax_depth"] = onto_df["tax_list"].apply(lambda x: 0 if x == [""] else len(x))

onto_df.to_csv("hmdb_no_pathways_ontology.tsv", sep="\t", index=False)
print(f"\n✅ Wrote hmdb_no_pathways_ontology.tsv  [{onto_df.shape[0]}×{onto_df.shape[1]}]")

print("\n📊 Taxonomy hierarchy depth statistics")
for d, n in onto_df["tax_depth"].value_counts().sort_index().items():
    print(f"  Depth {d}: {n:,} metabolites")

print("\n🔝 Top labels per hierarchy level")
for lvl in range(1, onto_df["tax_depth"].max() + 1):
    s = onto_df.loc[onto_df["tax_depth"] >= lvl, "tax_list"].apply(lambda lst: lst[lvl - 1])
    print(f"\n— Level {lvl} —  (n={len(s):,})")
    print(s.value_counts().head(10).to_string())

# preview file (first 15 rows)
prev = Path("hmdb_no_pathways_ontology_preview.txt")
with prev.open("w") as fh:
    for _, row in onto_df.head(15).iterrows():
        fh.write(f"\n{row.HMDB_ID}:\n")
        for col in onto_df.columns[1:]:
            if col == "tax_list":
                continue
            val = row[col]
            if pd.notna(val) and str(val).strip():
                fh.write(f"  • {col}: {val}\n")
print(f"👀 Preview saved to {prev}")

# ──────────────────────────────────────────────────────────────────
# 5. full label-count dump
# ──────────────────────────────────────────────────────────────────
level_counter = defaultdict(Counter)
for lst in onto_df["tax_list"]:
    if lst == [""]:
        continue
    for lvl, label in enumerate(lst, start=1):
        level_counter[lvl][label] += 1

rows = [{"Level": lvl, "Label": label, "Count": count}
        for lvl in level_counter
        for label, count in level_counter[lvl].items()]
counts_df = pd.DataFrame(rows).sort_values(["Level", "Count"],
                                           ascending=[True, False])
counts_df.to_csv("hmdb_taxonomy_label_counts.tsv", sep="\t", index=False)
print(f"📄 Wrote hmdb_taxonomy_label_counts.tsv  [{counts_df.shape[0]} rows]")

# ──────────────────────────────────────────────────────────────────
# 6.  *** NEW ***  annotation-gain summary
# ──────────────────────────────────────────────────────────────────
tot_no_pw = len(no_pw)
any_tax   = (onto_df["tax_depth"] >= 1).sum()
full4     = (onto_df["tax_depth"] >= 4).sum()

print("\n📈 Annotation coverage gained")
print(f"{'step':<35}{'# metabolites':>12}  comment")
print(f"{'input list without pathways':<35}{tot_no_pw:>12,}  these IDs lacked any HMDB pathway")
print(f"{'with any taxonomy string':<35}{any_tax:>12,}  ({any_tax/tot_no_pw:.0%}) at least kingdom→super-class")
print(f"{'with 4-level lineage':<35}{full4:>12,}  ({full4/tot_no_pw:.0%}) kingdom→…→subclass")

# ──────────────────────────────────────────────────────────────────
# 7.  *** NEW ***  depth-4 mapping column & TSV
# ──────────────────────────────────────────────────────────────────
onto_df["chem_subclass_lvl4"] = onto_df["tax_list"].apply(
    lambda lst: lst[3] if len(lst) >= 4 else pd.NA
)
onto_df.to_csv("hmdb_no_pathways_ontology.tsv", sep="\t", index=False)  # overwrite
onto_df[["HMDB_ID", "chem_subclass_lvl4"]].to_csv("hmdb_level4_mapping.tsv",
                                                 sep="\t", index=False)
print("🗂  Wrote hmdb_level4_mapping.tsv  (HMDB_ID → depth-4 subclass)")

# ──────────────────────────────────────────────────────────────────
# 8.  *** NEW ***  bar-plot of top depth-4 subclasses
# ──────────────────────────────────────────────────────────────────
lvl4 = counts_df[counts_df["Level"] == 4].nlargest(TOP_N, "Count")
plt.figure(figsize=(9, 6))
plt.barh(lvl4["Label"][::-1], lvl4["Count"][::-1])
plt.xlabel("Number of metabolites")
plt.title(f"Top {TOP_N} depth-4 chemical subclasses\n"
          "(for metabolites without HMDB pathway annotation)")
plt.tight_layout()
plt.savefig("taxonomy_depth4_top20.png", dpi=300, bbox_inches="tight")
plt.close()
print("🖼  Saved taxonomy_depth4_top20.png")

# ──────────────────────────────────────────────────────────────────
# 8-bis.  SECOND plot – exclude subclasses with >5 000 metabolites
#         and save as taxonomy_depth4_top20_under5k.png
# ──────────────────────────────────────────────────────────────────
lvl4_all   = counts_df[counts_df["Level"] == 4]          # full depth-4 table
lvl4_under = lvl4_all[lvl4_all["Count"] <= 5_000]        # keep only ≤5 000
lvl4_under = lvl4_under.nlargest(TOP_N, "Count")         # top-N after filter

plt.figure(figsize=(9, 6))
plt.barh(lvl4_under["Label"][::-1], lvl4_under["Count"][::-1])
plt.xlabel("Number of metabolites (≤ 5 000)")
plt.title(f"Top {TOP_N} depth-4 subclasses\n"
          "(after removing very large categories)")
plt.tight_layout()
plt.savefig("taxonomy_depth4_top20_under5k.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("🖼  Saved taxonomy_depth4_top20_under5k.png")

# ──────────────────────────────────────────────────────────────────
print("\n🎉 Ontology probe complete.")

