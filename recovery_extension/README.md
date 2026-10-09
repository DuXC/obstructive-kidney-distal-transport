# GSE96102 mouse obstruction and recovery

This animal-level reanalysis uses 71 matched-comparison arrays and retains four normal controls in the 75-sample inventory. The five human transport programs have fixed memberships. Seventeen genes have one-to-one mouse orthology; the three programs with complete orthology are testable. CLCNKA/CLCNKB map many-to-many, leaving TAL salt and DCT NaCl programs untestable.

## Reproduction

From the repository root:

```sh
python reproduce_recovery.py --figures --check-normalization
```

Use the Python dependencies in `requirements.txt` and R with limma 3.68.5. `RSCRIPT` selects the R executable and `KIDNEY_R_LIB` selects an installed package library. The command verifies input and fixed-plan hashes, reconstructs all reference result files and compares every retained value and text field against the reference results at tolerance 1e-10. Byte identity is recorded separately. The normalization option rebuilds the 43,379 by 71 matrix from foreground/background signals. Runtime checks are written locally by the entry point.

## Source data and measurements

`raw_data_sources.tsv` lists the 75 official source URLs, byte sizes and SHA-256 values. The original compressed arrays total 625,754,441 bytes. GPL4134 annotations, public sample metadata and dated Ensembl responses establish sample, probe and orthology mappings. Animals have unique filename identifiers and are independent biological replicates. `02_audit` contains scientific measurement and mapping records; `03_protocol` contains the dated analysis plan and its hash.

The Cy5 tissue channel receives normexp background correction (offset 50), log2 transformation and quantile normalization. This specifies the single-channel estimand used here; the source study used two-colour ratios. Gene estimates aggregate replicate spots within ProbeName and distinct unambiguous probes by median. Complete program membership is required per animal. Detectability sensitivity can change both contributing probes and eligible animals.

To download the source arrays and reconstruct upstream inputs in a separate working copy:

```sh
python recovery_extension/04_scripts/02_acquire_feature_tables.py
python recovery_extension/04_scripts/03_extract_cy5_inputs.py
python reproduce_recovery.py --check-normalization
```

Verify downloads against `raw_data_sources.tsv`. Source schema and quality-flag definitions are retained in the input records. An independent R implementation, `04_scripts/07_crosscheck_statistics.R`, checks Welch, four-group Satterthwaite and date-adjusted HC3 calculations, including their BH families; run it with this extension as the working directory.

## Statistical interpretation

The main families contain 50 program and 190 gene slots; date sensitivity contains 30 and 114 slots. Untestable slots contribute p=1 to the BH denominator and remain NA in displayed results. All planned matched-sham, direct late/acute and sham-adjusted temporal contrasts are retained.

Effects describe whole-kidney transcript-associated fluorescence. Acute arrays were acquired in 2012 and recovery arrays in 2014; temporal contrasts therefore require assumptions about acquisition era. Selected recovery comparisons lack sufficient shared scan dates for adjustment. DCT magnesium scores are lower in selected acute contrasts, whereas Trpm6 has no FDR-supported matched-sham effect. Human lineage-level RNA and mouse whole-kidney fluorescence are interpreted on their own scales.

See the root `DATA_LICENSE.md` for attribution and reuse terms. Project code is MIT; third-party data retain their source terms.
