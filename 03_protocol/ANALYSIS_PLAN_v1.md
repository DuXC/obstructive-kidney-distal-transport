# Analysis plan v1: distal transport in tumour-nephrectomy tissue

Frozen on 2026-09-09 before donor-level expression tests. Metadata counts and matrix linkage QC had been inspected. This is a local, time-stamped analysis plan, not a prospectively registered protocol.

## Question and population

Compare donor-level distal nephron transcriptional programs in five upper urinary tract urothelial carcinoma nephrectomies with ureteric obstruction and seven renal cell carcinoma nephrectomies without obstruction (GSE254185; Reck et al., 2025, doi:10.1038/s41467-025-59997-4). The estimand is the difference between these two tissue-sampling contexts. Tumour type and obstruction are perfectly confounded and cannot be separated by regression. RNA expression is not a measurement of urinary transport, electrolyte excretion, protein localisation, or reversibility.

## Data and annotation

Use the six official GEX H5 count files and author-retained barcode metadata. Verify SHA-256, H5 structure, nonnegative integer-valued counts, common gene order, one-to-one barcode-library-donor mapping, and conservation of counts after aggregation. The H5 files contain additional barcodes; only the author's final 46,957 nuclei are retained. Keep complete genes for pseudobulk and score reporting. Do not refit clustering from condition labels.

Retain public level-1 TAL (including author-labelled injured/inflammatory TAL and macula densa) and DCT as primary lineages. PC, IC-A and IC-B are exploratory. CNT and other lineages provide annotation context. Audit canonical markers using donor-balanced expression/detection summaries. Preserve author level-2 states and explicitly identify their condition-aware annotation as a limitation. Broad-lineage effects may include changes in constituent cell states, and are not within-state causal changes. Sensitivity analyses use the author's healthy TAL states (cTAL1/cTAL2/mTAL), healthy DCT states (DCT1/DCT2), and TAL excluding macula densa, where coverage permits.

## Biological replication and QC

Sum original RNA counts by donor × lineage, combining technical libraries. No nucleus is an independent patient. Primary eligibility: at least 20 nuclei in a donor-lineage and at least four donors in each group for testing. Donor points and actual n are mandatory. If a group has fewer than four eligible donors, report descriptive expression and no confirmatory test. Sensitivity thresholds are 10 and 50 nuclei. Leave-one-donor-out (LODO) analyses may have three donors in a group, clearly labelled sensitivity. Preserve all small-sample outcomes rather than change cutoffs after seeing significance.

## Prespecified modules

The authoritative list is `modules_v1.json`; `module_gene_table.tsv` records every gene and source. These are manually curated, small classical physiology programs, not complete GO or Reactome memberships and not validated clinical scores. Five primary comparisons: TAL salt transport, TAL paracellular calcium/magnesium, DCT NaCl, DCT magnesium handling, and DCT calcium handling. Four exploratory comparisons: PC water/vasopressin, PC sodium/potassium, IC-A acid secretion, IC-B bicarbonate handling. Two contextual programs (epithelial injury and chemokine/adhesion) are evaluated within each lineage. All memberships are fixed before tests. Negative regulators with opposite physiological direction are not combined into the positive transporter sets.

Normalize each lineage with edgeR TMM using expression-filtered genes while retaining total RNA library sizes. Module score = arithmetic mean of donor log2 CPM (edgeR prior.count=2) over the fixed members present in the RNA matrix; report membership and detection, and require at least 80% of prespecified genes plus at least two genes for a score. Do not remove a member because it is not significant or has low expression. This score describes transcript abundance; a unit difference is not a fold change in physiological flux.

## Primary and exploratory inference

Gene-level analysis: edgeR filterByExpr (default count criteria; design ~group), TMM, limma-voom, and robust empirical Bayes moderation. Contrast UUO minus Control; report log2 fold change, moderated 95% confidence interval, p value and BH FDR. Control FDR over all tested genes across TAL and DCT jointly, and separately over eligible exploratory lineages. Also retain within-lineage FDR. No integrated/scaled single-cell values enter these tests.

Primary module test: two-sided Welch comparison, difference in mean log2 CPM and 95% Welch interval. BH correction jointly over the five primary tests. Correct four exploratory transport tests as a separate family (untestable contrasts retain NA). Contextual injury/inflammatory comparisons form a separate family. Exhaustive label permutations of the difference in means are a sensitivity, with the observational exchangeability assumption stated. Significant = FDR <0.05; all results and uncertainty are reported.

Captured composition: donor fractions of all retained nuclei and, for epithelial lineages, of retained epithelial nuclei. Report percentage-point difference with a Welch interval and exact label-permutation p; BH separately for the five target lineages within each denominator. These are relative capture proportions, not tissue absolute abundance. Author injury-state fractions are descriptive since labels partly used UUO enrichment.

Covariation: fixed transport modules versus the five-gene epithelial injury and three-gene chemokine/adhesion modules in their corresponding lineage. Show donor scatterplots, pooled Pearson r, and group-adjusted partial Pearson r (95% Fisher interval; t approximation with n-3 df), with BH across tested partial correlations. Require at least eight donors and at least three per group. This is exploratory and cannot establish a mechanism or an association independent of all confounding. Do not interpret a pooled correlation driven only by group separation as independent covariation.

## Sensitivity and stopping choices

Repeat normalization, module estimates, and prespecified transport-gene estimates with thresholds 10/50 and each eligible donor omitted. Fit sex-adjusted expression models as sensitivity; add a joint sex + age-band-midpoint sensitivity only if full rank, acknowledging midpoint approximation and small residual df. Do not adjust away fibrosis/eGFR as presumed confounders: these could be related consequences and are systemic/whole-patient measurements. Summarize direction consistency and effect ranges rather than requiring every LODO p value to remain significant.

TAL/DCT healthy-state aggregation and TAL macula-densa exclusion are secondary composition probes with identical eligibility rules. Interpret healthy-only estimates cautiously because selection into an injury state depends on the phenotype being studied. Do not claim that this isolates a pure within-state obstruction effect.

## Planned figures and interpretation

1. Donor coverage, library mixing, marker audit and RNA QC.
2. Relative capture composition, author-state composition and primary module effects with donor points.
3. Fixed transport gene effects and exploratory collecting-duct modules.
4. LODO/threshold stability and group-adjusted transport-injury covariation.

Use Python/matplotlib at 183 mm width, editable SVG/PDF and 300 dpi PNG for review; every panel maps to a tabular source and script. Figures are preliminary research outputs with n, intervals and multiplicity definitions in legends.

After the four-figure package and Methods/Results draft, decide whether GSE282059 adds useful anatomical localisation. It is same-study spatial evidence, not independent external validation. No spatial or animal data are required to rescue a negative primary comparison.

## Literature increment

The original Figure S12 and Supplementary Data 5 already describe TAL/DCT/CNT/PC injury and include SLC12A1 changes. The intended increment is a transparent donor-level assessment of specified physiological programs with effect sizes, uncertainty, relative composition and stability. Targeted searches also found Smmit (2025; GSE254185 integration benchmark), a 2026 SPU obstruction-injury module abstract, and prior animal aquaporin/distal injury work. No universal novelty claim is made. Full search responses and appraisal are kept with sources.
