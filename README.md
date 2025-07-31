# Turkana-Metabolomics
An analysis focused on determining differential metabolite abundance in individuals leading divergent lifestyles. 
# Which metabolites are significantly affected by the lifestyle of the individuals? 

# One-Pager: How to Run the Metabolomics Pipeline 🚀
1. The files that you drop in the folder alongside the scripts
      a) HILIC_RAW_DATA_with_covariates.csv	, Raw intensities + sample info (Age / Sex/batch / MW_scores_lifestyle).
      b) Updated_No_match_Metabo.csv	Two columns → Query (raw name) → HMDB (accession) for renaming.
      c) serum_metabolites_HMDB.xml	Full HMDB database to harvest pathways and ontology.
2. Scripts Rename_to_HMDB.py, ENRICHMENT_HMDB.py, STATISTICS_ENRICHMENT.py, PLAY_CODE.py, ANNOTATE_THE_REST.py

3. Fire it up 🔥
# classic: only “real” pathways
ENRICHMENT_HMDB.py
# OPTIONAL: fold in level-4 chemical subclasses for pathway-less metabolites
python ENRICHMENT_HMDB.py --use_ontology
Everything chains automatically—ENRICHMENT_HMDB.py orchestrates the lot.

# What each step does
1. Rename_to_HMDB.py	Renames metabolite columns to HMDB IDs → kills outliers → median-scales & log-transforms → regresses out Age/Sex/batch → writes HILIC_HMDB_Residualized_FINAL.csv.
2. ENRICHMENT_HMDB.py	• Calls the step above -> Parses HMDB XML → ties each metabolite to its pathways.
     If --use_ontology is supplied:
     runs ANNOTATE_THE_REST.py to grab depth-4 subclasses for any IDs missing a pathway;
     tags them as CHEM:<Subclass> and appends them to the mapping.
     Saves user_matched_HMDB_pathways.tsv.
3. ANNOTATE_THE_REST.py (only when invoked)	--> pulls taxonomy for pathway-less IDs, writes hmdb_level4_mapping.tsv, and a bar plot of top subclasses.
4. PLAY_CODE.py	Auto-run sanity check: bubble plots showing pathway sizes vs metabolite counts.
5. STATISTICS_ENRICHMENT.py	ANOVA (Urban vs Non-Urban) → flags sig metabolites → Fisher enrichment → outputs TSV + dual-panel scatter + box-plots for top hits.
6. Outputs you care about
   A) HILIC_HMDB_Residualized_FINAL.csv	Clean, residualised matrix.
   B) user_matched_HMDB_pathways.tsv	HMDB ID → pathway list (plus CHEM: Subclass rows if ontology used).
   C) hmdb_level4_mapping.tsv	(only with --use_ontology) fallback depth-4 subclass map.
   D) pathway_enrichment_results.tsv	p/q values, hit counts, dominant group.
   E) enrichment_dual_panel.png	Raw-p vs FDR-q scatter.
   F) box-plots for individual top pathways/subclasses.

That’s it: drop the inputs, pick your flavour (with or without ontology), run a single command, harvest tables and figures. Happy enriching! 🌟
