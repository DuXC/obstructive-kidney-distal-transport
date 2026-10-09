# Distal nephron transport transcription

Code and data for **Distal Nephron Transport Transcription in Tumour Nephrectomy and Experimental Obstruction**.

Version **0.7.0** combines three reproducible analyses: the GSE254185 tumour-nephrectomy comparison, cross-injury comparisons in the Lake kidney atlas, and the GSE96102 mouse obstruction and recovery series. The unit of replication is the human donor or individual animal. All fixed programs, constituent-gene results, confidence intervals and sensitivity estimates are retained.

## Reproduce the analyses

Use Python 3.12 and R with the dependencies in `requirements.txt`, `R_packages.tsv` and `environment/`. From the repository root:

```sh
python reproduce.py
python reproduce_external.py
python reproduce_recovery.py --figures --check-normalization
```

`RSCRIPT` can point to an Rscript executable; `KIDNEY_R_LIB` can point to an installed R library. The entry points use the supplied biological-sample aggregates and check their hashes. Discovery and external human outputs are compared with expected tables at numerical tolerance 1e-10. Mouse reference hashes are verified and regenerated effects are compared at tolerance 1e-10, with byte identity reported separately; the optional normalization check reruns limma on the included foreground/background matrices. Execution logs and checks are written during the local run. Dependencies are installed separately.

## Data and figure map

| Analysis | Inputs and programs | Results | Current figures |
|---|---|---|---|
| Human nephrectomy | `02_data/derived`, `03_protocol` | `expected/05_results` and `expected/09_manuscript_development_20260909/02_analysis` | Main 1–3 and S4 in `10_manuscript_refinement_20260909/05_figures`; S1–S3 in `09_manuscript_development_20260909/05_figures` after reproduction |
| Human external atlas | `external_validation/02_data`, `external_validation/03_protocol` | `external_validation/expected` | Main 4 and S5–S6 in `external_validation/06_figures` after reproduction |
| Mouse obstruction | `recovery_extension/01_sources`, `recovery_extension/02_audit`, `recovery_extension/03_protocol` | `recovery_extension/05_results` | Main 5 and S7 in `recovery_extension/06_figures` |

Main Figure 4 displays primary single-nucleus and secondary single-cell estimates side by side. Numerical tables remain the source of the plotted effects. Figures use editable SVG/PDF and print-resolution TIFF exports.

## Study populations and measurements

GSE254185 contributes 46,957 nuclei from 12 donors: five obstructed UTUC and seven unobstructed renal tumour nephrectomies. The source recruitment methods describe six renal cell carcinomas and one oncocytoma in the reference group; individual histology is unlinked to donor IDs. `COHORT_ANNOTATION_20261005.md` documents the source discrepancy and its resolution. The original dated analysis plan is retained as scientific provenance. The comparison estimates differences between the two tumour-nephrectomy contexts.

The independent Lake atlas contributes primary single-nucleus GSE183277 and secondary single-cell GSE183276 comparisons. Fourteen reported patients shared across assays were excluded from the secondary assay. Fractional values in the author-processed single-cell count matrix are retained, and technical libraries are combined by reported patient ID. The assays come from one atlas study. Source counts, eligibility and representation are documented in `external_validation/provenance` and `external_validation/03_protocol`.

GSE96102 contributes 71 animals to time-matched comparisons, plus four normal controls documented in the source inventory. The analysis uses Cy5 foreground/background signals, normexp correction and log2 quantile normalization. Three fixed programs have complete one-to-one orthology. The two chloride-dependent programs are untestable because CLCNKA/CLCNKB orthology is many-to-many. Whole-kidney fluorescence is a different measurement scale from lineage-level RNA scores. Acquisition dates are included for interpretation of acute and recovery comparisons.

All twelve current figure exports are also supplied in `publication_figures/` (SVG, PDF, PNG and 600-dpi TIFF).

## Source data and citation

Downloadable [release v0.7.0](https://github.com/DuXC/obstructive-kidney-distal-transport/releases/tag/v0.7.0) includes:

- `Source_Data_1.xlsx`: discovery donor characteristics, program definitions and statistical summaries.
- `Source_Data_2.tsv.gz`: complete discovery gene-level estimates.
- `Source_Data_3.zip`: external human donor coverage, scores, gene/program estimates and source provenance.
- `Source_Data_4.zip`: mouse sample identities, orthology, fluorescence inputs and all planned contrasts.
- `Analysis_Code_and_Data_v0_7.zip`: reproducible code and analysis inputs.

Source studies and accession records:

- Reck et al. Nature Communications (2025), [doi:10.1038/s41467-025-59997-4](https://doi.org/10.1038/s41467-025-59997-4), [GSE254185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254185).
- Lake et al. Nature (2023), [doi:10.1038/s41586-023-05769-3](https://doi.org/10.1038/s41586-023-05769-3), [GSE183277](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183277) and [GSE183276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183276).
- Wu and Brooks, [GSE96102](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96102).

Official raw-source URLs and checksums are supplied in `raw_data_sources.tsv`, `external_validation/provenance` and `recovery_extension/raw_data_sources.tsv`. Raw-source download and reaggregation scripts accompany the aggregate entry points. The external full-matrix RDS files total approximately 2.32 GB and are retrieved separately.

Use `CITATION.cff` to cite this software and cite the original studies when reusing their data. Project code is licensed under MIT; source-derived data follow `DATA_LICENSE.md`.

## Authors and support

Xiancheng Du, Yongkun Zhu, Kaihua Xue, Shuchun Tao, Chunhui Liu and Chao Sun. Chunhui Liu and Chao Sun are co-corresponding authors. Support: China Postdoctoral Science Foundation 2024M750457 (Chunhui Liu), and Jiangsu Provincial Research Project on Traditional Chinese Medicine and Integrated Chinese-Western Medicine ZXFZ2026021 (Chao Sun). All authors declare no competing interests.
