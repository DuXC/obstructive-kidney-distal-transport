# Reproduction from donor aggregates and public mouse arrays, version 0.6

This code and data release supports the analysis for “Distal nephron transport programs in cancer nephrectomy tissue and experimental obstruction”. It starts from verified donor-by-lineage and donor-by-state RNA count matrices derived from GSE254185. The main comparison is five UTUC nephrectomies with obstruction versus seven RCC reference nephrectomies. Tumour context is inseparable from obstruction.

## Run

Use Python 3.12 and R with the packages recorded in requirements.txt, R_packages.tsv and R_sessionInfo_original.txt. From the extracted package directory, run:

```sh
python reproduce.py
```

RSCRIPT may identify an Rscript executable. KIDNEY_R_LIB may identify an existing R library. Dependencies are not installed automatically. The statistical code is CPU based. A fast statistics-only run is available as `python reproduce.py --skip-figures`.

The entry point checks input and frozen-plan hashes, reruns TMM/limma-voom gene models, all prespecified program sensitivities, composition and contextual associations, then the explicitly post hoc member-influence diagnostic. Marker summaries pool technical libraries within each donor using nucleus weights before calculating equal-donor averages. It checks eleven output tables against expected snapshots using exact strings and numeric tolerance 1e-10. It reports byte-identical files separately. Current Fig1, Fig2, Fig3 and FigS4 are regenerated under 10_manuscript_refinement_20260909/05_figures. Current FigS1, FigS2 and FigS3 are regenerated under 09_manuscript_development_20260909/05_figures, which also retains the earlier layouts. All figures are exported as SVG, PDF, PNG and TIFF. Figure byte identity is not expected across font/rendering environments. The run writes reproduction_receipt.json only after validation succeeds.

## Included data and provenance

- 02_data/derived: donor count matrices, sample mapping and public retained-nucleus QC data.
- 05_results: source-derived clinical characteristics, captured composition, state counts and marker summaries needed to initialize the models and figures.
- 03_protocol: the original locally timestamped, frozen plan and fixed gene membership, with original hashes. This is a local pre-analysis freeze, not prospective registration.
- expected: reference output tables retained for verification, not statistical inputs.
- 04_scripts/02_link_qc_pseudobulk.py: original raw-RNA linkage and aggregation code for source inspection and a separate full-raw rerun.
- raw_data_sources.tsv: official public RNA/metadata URLs, destinations and verified source checksums. Raw source download and reaggregation are outside the default entry point; the original full raw pipeline was verified in the originating workspace. The release acceptance run starts from the included donor aggregates.

Source: Reck M et al. Multiomic analysis of human kidney disease identifies a tractable inflammatory and pro-fibrotic tubular cell phenotype. Nature Communications 16, 4745 (2025). https://doi.org/10.1038/s41467-025-59997-4. Public data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254185. Clinical variables were extracted from the source Supplementary Data 2. Source attribution must accompany reuse. Original project code is provided under the MIT license. See DATA_LICENSE.md for data provenance and reuse terms.

The verified main results are transcript scores, not physiological flux. Author state selection and captured nuclear composition are retained as interpretation limits. Member omission is a post hoc influence diagnostic; it does not redefine the five primary programs. The separate same-study GSE282059 CosMx panel assessment is supplied in the manuscript supplement and cannot validate the three main programs because key members were unmeasured.


## Associated manuscript metadata

Xiancheng Du, Yongkun Zhu, Kaihua Xue, Shuchun Tao, Chunhui Liu and Chao Sun, in that order. Chunhui Liu and Chao Sun are co-corresponding authors; Chao Sun is the final author. The manuscript is supported by China Postdoctoral Science Foundation 2024M750457 (Chunhui Liu) and Jiangsu Provincial Research Project on Traditional Chinese Medicine and Integrated Chinese-Western Medicine ZXFZ2026021 (Chao Sun). The authors declare no competing interests.


## Version and citation

Version 0.6 preserves the discovery and external human analyses and adds the frozen GSE96102 mouse obstruction/recovery extension. The associated manuscript remains under author review.

Repository: https://github.com/DuXC/obstructive-kidney-distal-transport

Versioned release: https://github.com/DuXC/obstructive-kidney-distal-transport/releases/tag/v0.6.0

Use CITATION.cff to cite this software. Cite Reck et al. (2025) and GSE254185 for discovery data, and Lake et al. (2023), GSE183277 and GSE183276 for the external atlas data. The source and supplementary code paths retain their original phase names for reproducibility.


## External public-data extension

Run python reproduce_external.py from the repository root. This checks the included external donor-aggregate and plan hashes, reruns the two assay analyses, independently verifies Welch intervals and ten-test BH correction, compares fourteen regenerated statistical tables against frozen expected snapshots, and redraws Fig4/FigS5/FigS6. Set RSCRIPT and KIDNEY_R_LIB as needed, using the same packages as the discovery analysis. External environment snapshots are included.

The primary snRNA analysis includes four non-stone reference, nine AKI and seven CKD cortical donors. None of its ten tests passes FDR correction. In the secondary scRNA assay, AKI TAL salt and DCT magnesium pass the separate ten-test family (FDR 0.011 and 0.003); CKD DCT has only two affected donors and remains descriptive. Reference-state restriction attenuates scRNA TAL support. These are partial cross-injury findings, not obstruction-specific replication.

The GSE183276 matrix is author-processed count-scale data with fractional values in a source reporting SoupX preprocessing. Values are not rounded. Fourteen reported patient IDs shared with snRNA are excluded; the three PRE019 libraries are combined. Both assays come from the same atlas and are not independent studies. All representation, donor-mapping and inference decisions are recorded in external_validation/03_protocol.

The acceptance run starts from the included donor aggregates. The optional source-download and aggregation scripts require the official source RDS/metadata files and greater memory; the two external RDS archives total about 2.32 GB and are not duplicated here. Download receipts identify the exact official URLs and checksums. The original whole-matrix linkage and aggregation were verified in the originating workspace.

Source: Lake BB et al. An atlas of healthy and injured cell states and niches in the human kidney. Nature 619, 585–594 (2023). https://doi.org/10.1038/s41586-023-05769-3. GEO GSE183277 and GSE183276; reference-biopsy provenance also uses GSE169285. Source participants are identified only by public study pseudonyms.


## Public mouse obstruction/recovery extension

Run `python reproduce_recovery.py` from the repository root. This checks the frozen plan, source metadata, orthology and normalized-fluorescence input hashes, reconstructs animal-level gene/program values, and compares 13 regenerated result files with exact reference hashes. It does not use reference effects as statistical inputs. Add `--figures` to regenerate Fig5/FigS7. See [recovery_extension/README.md](recovery_extension/README.md) for the optional upstream normalization and official raw-download steps.

GSE96102 contains 75 individual mice: 36 affected, 35 matched shams and four normal controls. The 71 matched-comparison arrays enter the fixed single-channel Cy5 pipeline; this is a new raw-fluorescence reanalysis, not a reproduction of the authors' two-colour ratios. Normexp background correction (offset 50) and log2 quantile normalization retain 43,379 common noncontrol features. Seventeen of 19 human target genes have one-to-one mouse orthologues. CLCNKA/CLCNKB are many-to-many, so the two affected complete programs remain untestable.

The DCT magnesium score is lower in both acute histories and in the primary five-day/day-10 comparison. The acute findings survive eligible acquisition-date adjustment. Trpm6 has no FDR-supported matched-sham difference. Some transporter estimates depend on detectable-probe aggregation. Acute arrays were acquired in 2012 and recovery arrays in 2014; several later contrasts lack adequate shared acquisition dates. Direct temporal effects are era-confounded, sham-adjusted changes require an untestable group-by-era assumption, and day-28 uncertainty is not equivalence or proof of physiological recovery. Whole-kidney scores do not identify within-lineage regulation. All negative, untestable and sensitivity results are included.

Independent R calculations reproduced 427 tested contrasts (Welch, four-group Satterthwaite and HC3 models) at absolute tolerance 1e-8. The acceptance run for this extension starts from verified normalized fluorescence, not from a new raw-data acquisition.
