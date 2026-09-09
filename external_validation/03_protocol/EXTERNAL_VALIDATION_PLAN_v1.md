# External public-data validation plan v1

This local plan is frozen on 2026-09-10 after inspection of source metadata and donor/lineage coverage, before reading any external count matrix or calculating any external gene/program contrast. It extends the existing five-program discovery analysis; it is not a prospectively registered protocol.

## Question and evidence tier

GSE183277 (Lake et al. human kidney atlas, snCv3) is independent of the original Edinburgh GSE254185 cohort. It tests transportability to AKI and CKD, not independent replication of an obstruction-specific effect. The atlas source already describes altered distal states. Neither cross-disease agreement nor significance eliminates the original tumour/obstruction confounding.

Native GDS searches retrieve candidate series by human kidney/injury/obstruction terms; the broad injury query is capped at 100 records and is targeted feasibility work, not an exhaustive systematic review. GSE180394 has one hydronephrosis sample, below the minimum of four affected donors. GSE183276 shares at least 14 patient IDs with GSE183277. Its reference identifiers need source lineage checking before a claim of independent patients.

## Fixed primary cohort and annotation

Use GSE183277 official snCv3 count RDS and accompanying 2022-03-28 metadata. Keep author-retained nuclei; trim whitespace in patient IDs; verify unique barcode mapping and one condition/sex per donor. Sum untransformed RNA counts across technical libraries for each patient and broad author lineage.

The primary anatomical set consists of specimens with at least 50% cortex according to percent.cortex. Exclude every specimen annotated condition.l3=Stone from the normal-reference group. Keep AKI and CKD separate; AKI includes the source COV-AKI category. Labels and tissue procurement remain source-defined. Metadata identifies four non-stone cortical reference donors, seven CKD donors and nine AKI donors before lineage-specific coverage. Reference tissue procurement and disease tissue procurement differ; this is a limitation, not a covariate that can necessarily be identified.

Primary lineages are subclass.l1=TAL and DCT. All author states within each lineage are retained. The source annotation is not an independent assay of injury; no biological causal inference from state labels.

## Fixed programs and normalization

Reuse exactly the five primary_transport programs from the original modules_v1.json, without changing members or sign. Require at least 80% of fixed genes represented plus at least two genes. Show measured and undetected members explicitly.

Require at least 20 nuclei per patient-lineage and at least four eligible donors per comparison group. Use the same lineage normalization convention as discovery: expression-filtered genes for edgeR TMM, original total library sizes retained, log2 CPM with prior.count=2 on all genes. Scores are arithmetic means across fixed members and measure transcript abundance, not physiological flux. Compute each contrast independently using its eligible disease and reference donors.

## Tests and multiplicity

Primary effects are AKI minus reference and CKD minus reference. Two-sided Welch means with 95% Welch intervals; BH across all ten prespecified program-by-disease comparisons (untestable entries remain NA but count toward the planned family size of ten). Report all five programs for both conditions.

Exploratory gene-level follow-up uses limma-voom, robust empirical Bayes, and a single BH family across all tested genes in both lineages and both primary contrasts. Display all fixed members, not a selected significant list. Gene estimates and unweighted program scores are distinct.

## Prespecified sensitivity and interpretation

1. Primary anatomy set at 10 and 50 nuclei per donor-lineage.
2. All non-stone cortex-containing specimens (percent.cortex > 0), keeping the 20-nucleus rule.
3. Author reference-state nuclei only (state.l2=reference), using identical coverage rules; this probes state mixture and is not a causal within-state estimator.
4. Exclusion of COV-AKI from the AKI contrast.
5. Sex-adjusted program linear models only for full-rank designs with adequate residual degrees of freedom.
6. Donor omission for the primary programs; three donors per group permitted only for this sensitivity. Describe direction and ranges rather than requiring significance in every omission.

Sensitivity estimates are diagnostic and not alternative routes to declare validation. Report the primary family first. State composition is descriptive and denominator-specific.

## Optional second assay

GSE183276 scCv3 may be analyzed only after verifying its count representation and donor identifiers against original GSE169285/GEO records. Patient overlap with snCv3 must be removed for a claim of independent-donor support; unresolved donor identities cannot be counted separately. Its two disease contrasts form a separate ten-test secondary family with the same programs and minimum coverage. It does not constitute an independent study if sourced from the same atlas, even when donor sets do not overlap.

## Outputs

Acquisition and count/metadata QA receipts, full donor coverage and membership tables, all primary and sensitivity estimates, effect plots with donor-level n and 95% intervals, source-linked Chinese report and manuscript addendum. Keep discovery results and the published v0.4.0 release immutable; any release update receives a new version.

