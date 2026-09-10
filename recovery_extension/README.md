# GSE96102 public mouse obstruction and recovery

The primary analysis unit is one animal. This extension retains all 75 source samples in the inventory and uses 71 matched-comparison arrays for normalization and inference. Four normal animals are source/QC records only. The five frozen human transport programs retain their original membership; 17 genes have unambiguous one-to-one mouse orthology. CLCNKA and CLCNKB map many-to-many, leaving TAL salt and DCT NaCl programs untestable in this exact-membership analysis.

## Reproduction

From the repository root, use the versions in requirements.txt and an R installation with limma 3.68.5:

```sh
python reproduce_recovery.py --figures --check-normalization
```

RSCRIPT selects the Rscript executable and KIDNEY_R_LIB can select an installed R package library. No dependencies are installed automatically. The default command (with no options) requires only Python, checks the input and frozen-plan hashes, reconstructs the gene/program summaries from the included normalized fluorescence and target quality flags, and reproduces 13 result files byte-for-byte in the recorded environment. The expected hashes are validators, never effect inputs. It writes regenerated results under reproduced_results and a receipt under validation.

The optional normalization check runs limma in a temporary directory from the supplied foreground/background matrices, compares all 43,379 x 71 normalized values within absolute tolerance 1e-10, and leaves the frozen inputs intact. This verifies the normalization stage; it is not a fresh acquisition of the original Agilent files. The optional figure run uses the regenerated tables. PDF/SVG metadata and font differences can change figure hashes on another system.

## Source acquisition and upstream reconstruction

raw_data_sources.tsv lists the 75 official source URLs, byte sizes and SHA-256 values. The source compressed arrays total 625,754,441 bytes and are excluded from Git. The preserved GPL4134 annotation, public sample metadata and dated Ensembl responses supply the mapping provenance. Metadata describe individual biological replicates with unique filename animal IDs; no samples are repeatedly measured animals.

In a separate working copy of the repository, the following commands download/extract the source arrays and rebuild the pre-normalization inputs:

```sh
python recovery_extension/04_scripts/02_acquire_feature_tables.py
python recovery_extension/04_scripts/03_extract_cy5_inputs.py
python reproduce_recovery.py --check-normalization
```

Check every downloaded source file against raw_data_sources.tsv before interpreting a reconstruction. The code never treats absent optional flags as measured zeroes. The original 43-column and 112-column schemas, shared feature/probe identities, and single-channel decision are documented in 02_audit and the frozen protocol. GPL annotation and Ensembl mapping snapshots are fixed inputs rather than a request to silently accept future annotation changes. The released default acceptance run does not claim a fresh download/re-extraction; the original source extraction is documented in 09_qa/raw_acquisition_receipt.json and cy5_input_receipt.json.

For the independent statistical calculation, run Rscript 04_scripts/07_crosscheck_statistics.R with this extension directory as the working directory. It compares separate Welch t.test, four-group Satterthwaite and date-adjusted HC3 calculations against the supplied animal values/results, including BH families, and reproduces 427 tested contrasts within absolute tolerance 1e-8.

## Interpretation and complete reporting

The raw rMedianSignal and rBGMedianSignal kidney channel receives normexp correction (offset 50) and log2 quantile normalization. The universal-reference channel is not used. This is a new single-channel reanalysis, not a reproduction of the source two-colour ratio pipeline. Gene estimates first aggregate replicate spots within ProbeName and then distinct unambiguous gene probes by median. Complete program membership is required for each animal. The detectability sensitivity can change the probes contributing to a gene estimate, in addition to removing unmeasured samples.

Families contain 50 program and 190 gene slots; date sensitivity has 30 and 114 slots. Untestable slots use p=1 only in the BH denominator and remain NA in displayed results. Same-time sham comparisons, direct late/acute differences, and sham-adjusted changes are all retained with separate status/caveat fields. The inference is about whole-kidney transcript-associated fluorescence, not transporter activity or within-lineage regulation. Acute acquisitions are from 2012 and recovery acquisitions from 2014. Direct temporal changes cannot separate era from biology; sham-adjusted change requires an untestable group-by-era assumption. Several later comparisons lack enough shared scan dates for adjustment. Uncertain late contrasts do not establish equivalence to sham.

The DCT magnesium aggregate supports selected acute contrasts, but Trpm6 has no FDR-supported matched-sham effect and several individual transporter estimates depend on probe detectability. All positive, negative, uncertain and untestable comparisons are supplied. Human atlas and mouse fluorescence effect scales are not pooled.

PXD039314 sample mapping and GSE145053 expression-scale/recovery-label clarification were requested from the source authors on 10 September 2026. They remain resource dependencies and are not counted as completed protein or additional RNA validation in this release. Private correspondence is not part of this public package.

See the root DATA_LICENSE.md for source attribution and reuse terms. Code is MIT; third-party data retain their source terms.
