# Reproduction from donor aggregates, version 0.4

This code and data release supports the analysis for “Distal nephron transport programs in obstructed and unobstructed cancer nephrectomy tissue”. It starts from verified donor-by-lineage and donor-by-state RNA count matrices derived from GSE254185. The main comparison is five UTUC nephrectomies with obstruction versus seven RCC reference nephrectomies. Tumour context is inseparable from obstruction.

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

Version 0.4 preserves the verified numerical analysis and figures from version 0.3 and adds public release metadata, a code license and citation information. The associated manuscript remains under author review.

Repository: https://github.com/DuXC/obstructive-kidney-distal-transport

Versioned release: https://github.com/DuXC/obstructive-kidney-distal-transport/releases/tag/v0.4.0

Use CITATION.cff to cite this software and cite Reck et al. (2025) and GSE254185 for the source data. The source and supplementary code paths retain their original phase names for reproducibility.
