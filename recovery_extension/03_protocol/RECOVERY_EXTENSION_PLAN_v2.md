# Public recovery extension analysis plan v2

This local amendment was written on 10 September 2026 after metadata, probe identity, raw-file schemas, quality-flag availability and acquisition dates were inspected, and before any target gene or program group-effect calculation. It completes the earlier public follow-up plan. The discovery and primary human-atlas analyses remain frozen.

## Biological question and scope

GSE96102 provides 75 independent mouse kidney samples: 36 obstruction/reversal samples, 35 time-matched sham samples and four normal controls. The 1-day obstruction history has 6, 5 and 6 affected animals at recovery days 0, 10 and 28, with 6 sham animals at each time. The 5-day history has 6, 8 and 5 affected animals, with 5, 6 and 6 sham animals. The four normal controls provide descriptive QC only and are excluded from hypothesis tests and normalization for the matched comparisons. Animal identity is taken from the source's explicit biological-replicate description and unique sample filename identifier; it is not inferred from molecular similarity.

The estimand is whole-kidney transcript signal. This cannot establish within-TAL/DCT regulation, measured transport function or obstruction-specific effects in humans. Quantile normalization assumes broadly comparable array-wide distributions and does not estimate total RNA or tubular mass.

## Identity and measurement decision

The series specifies two-colour Agilent GPL4134 arrays, with kidney RNA labelled Cy5 and universal reference RNA labelled Cy3. The processed SOFT table uses 44,996 sequential identifiers; direct matching to the feature numbers in the first raw table reproduced only approximately 6.3% of values within rounding tolerance. This does not establish the author's intended processed-table map. Therefore those processed values will not be used for quantitative effects.

The public feature tables retain ProbeName, FeatureNum, red-channel median foreground/background intensities, and core quality flags. All 75 tables were acquired. Forty-eight tables have a reduced 43-column schema and lack surrogate/found flags; 27 have 112 columns. Absent flags remain missing.

To use the directly identified kidney measurement and one consistent preprocessing route across both schemas, analyse the red Cy5 channel as a single-channel experiment. Read rMedianSignal and rBGMedianSignal from each raw feature table. Retain noncontrol feature-number/probe-name pairs present in all 71 matched-comparison samples; require consistent identities and finite foreground/background values. Apply limma normexp background correction with offset 50, then log2 quantile normalization across these 71 arrays. This is a new reprocessing of raw fluorescence, not a reproduction of the author's Cy5/Cy3 pipeline. The universal-reference LogRatio and processed-sheet P values will not supply effect tests.

The background-correction method, offset, normalization and contrasts are fixed here; they will not be selected using target significance. Mask target spots flagged saturated, feature nonuniform/population outliers, or manually flagged in the red channel. Do not impute missing target values. Collapse replicate spots by median within each ProbeName, then distinct unambiguously annotated probes by median within each mouse gene. Report probe counts, quality masks and above-background fractions.

## Orthology and program definition

Query Ensembl human-to-mouse orthology and retain its stable IDs and returned mouse symbols. Seventeen of the 19 fixed primary genes have one-to-one orthologues. CLCNKA and CLCNKB have many-to-many relationships to Clcnka and Clcnkb. Consequently the original TAL salt and DCT NaCl programs are not evaluable as exact one-to-one cross-species programs. Do not silently drop or substitute their ambiguous members. The remaining three fixed programs are eligible only if all original members map and are measured.

Use the platform's probe-name annotation, requiring each retained ProbeName to map to exactly one mouse gene symbol. A program score is the equal-weight mean of its fixed members' normalized log2 signals. A sample must have all members; no available-member mean. Require at least four eligible independent animals in every group contributing to a contrast.

## Fixed comparisons and multiplicity

Analyse the two obstruction histories separately. For each history, use five contrasts: (1) affected day 0 minus its time-matched sham; (2) affected day 10 minus its time-matched sham; (3) affected day 28 minus its time-matched sham; (4) affected day 28 minus affected day 0; (5) the change in affected-versus-sham contrast from day 0 to day 28. Thus there are ten contrasts across five fixed programs, a single 50-slot primary recovery family. Use means, independent-group Welch/Satterthwaite standard errors and 95% confidence intervals. For the four-group difference in differences, add the four independent mean variances with their signed coefficients and use the Satterthwaite degrees of freedom. BH adjustment retains all 50 slots; an untestable slot enters the denominator as p=1 and is displayed as NA, never as an observed test.

Report all 19 target genes across the same ten contrasts in a separate 190-slot family. The two ambiguous orthology genes remain not evaluable. SLC12A1, SLC12A3 and TRPM6 receive figure emphasis because they were central in the existing manuscript, while every other planned gene result is delivered. This gene family was specified before recovery effects were inspected.

## Quality and acquisition sensitivity

Repeat the program and gene analyses after requiring at least one clean red-channel probe spot above background for each contributing gene in each animal. Apply a separate family of the same fixed size. Report untestable comparisons when this leaves fewer than four animals. This addresses signal detectability and does not replace the primary results.

Acute samples were acquired in 2012 and recovery samples in 2014. Direct late-versus-acute comparisons are fully confounded with acquisition era. They remain explicitly descriptive of the measured contrast, with no isolated recovery-effect claim. Difference in differences uses contemporary shams but assumes no group-by-era technical interaction, which these data cannot test.

For each of the six same-time affected-versus-sham contrasts, audit acquisition day overlap. Run a date-adjusted linear sensitivity only on acquisition dates containing both groups, only when at least four animals per group remain and residual degrees of freedom are at least four. Include date fixed effects and use an HC3 robust covariance estimate. Keep separate 30-slot program and 114-slot gene correction families. Do not relabel fully confounded dates as random batch noise or merge dates to make a comparison eligible. The 5-day history/day-28 comparison has no shared acquisition day and therefore has no date-separated group estimate.

## Reporting and release

Retain null estimates, mapping failures, quality failures, acquisition confounding and all eligible sample counts in source tables and the manuscript supplement. Report this as an exploratory public mouse recovery extension, subordinate to the frozen human discovery and primary external results. The RNA and fluorescence scales must not be pooled. Re-run the same locked scripts from verified inputs, cross-check the statistical calculations independently and inspect final figure/document renders before releasing the revised manuscript.
