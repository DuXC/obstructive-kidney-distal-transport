# Data provenance and reuse

Original project code is covered by LICENSE. Third-party data are attributed separately and retain their applicable source terms.

The included donor aggregates, source-derived nucleus and clinical annotations, and result tables originate from the public GEO series [GSE254185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254185) and the associated study by Reck et al., *Multiomic analysis of human kidney disease identifies a tractable inflammatory and pro-fibrotic tubular cell phenotype*, Nature Communications 16, 4745 (2025), https://doi.org/10.1038/s41467-025-59997-4. The source article is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Clinical fields were transcribed from its Supplementary Data 2. The source article's ethics approval applies to the original collection.

Changes in this project comprise donor/lineage aggregation, fixed-program scoring, statistical estimation, sensitivity analysis and figure preparation. Original derived tables and documentation contributed by this project are shared under CC BY 4.0 to the extent the project authors hold rights in them. This notice does not relicense third-party content. Cite the source study and accession alongside this release, retain provenance, and identify further modifications.

The default reproduction starts from the included donor aggregates. raw_data_sources.tsv identifies the original downloads and checksums for a separate raw-data reconstruction. Source sample labels are study pseudonyms, and no re-identification should be attempted.


The external_validation extension derives donor aggregates and annotations from public GEO GSE183277 and GSE183276 in Lake BB et al., An atlas of healthy and injured cell states and niches in the human kidney, Nature 619, 585–594 (2023), https://doi.org/10.1038/s41586-023-05769-3. The source article is licensed under CC BY 4.0. Public source data retain their source-specific terms; this release does not relicense third-party content. Changes consist of specified anatomical selection, donor/lineage aggregation, fixed-program scoring, statistical estimation and figure preparation. The scRNA values are author-processed count estimates. Source-data URLs and hashes, patient mapping and exclusions are retained for attribution and reproducibility.
