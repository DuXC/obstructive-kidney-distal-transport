# scCv3 representation review before secondary comparisons

The original frozen external plan required source representation and donor identity review before the optional second assay. No scCv3 program comparison has yet been calculated.

The GSE183276 RDS contains nonnegative count-scale values with fractional entries. Its GEO sample records state SoupX v1.5.0 ambient-RNA preprocessing. The SoupX author implementation returns fractional corrected counts by default (roundToInt=FALSE): https://github.com/constantAmateur/SoupX/blob/master/R/adjustCounts.R .

There are 149,561,540 fractional stored values among 174,049,344 stored entries; no negative or nonfinite values were found. Cell library totals vary from 616.8793 to 50,744.8031 and correlate with metadata nCount_RNA at r=1. In target donor-lineage aggregates, relative differences from the metadata count sums range from -1.81e-7 to 3.58e-8, compatible with finite-precision export. These checks support an author-processed count scale, rather than an integrated or log-normalized assay. The exact processing history of each numeric entry is not independently reconstructed.

For the optional secondary analysis, aggregate these supplied values without rounding, use the already planned TMM/log2 CPM program summaries and Welch tests, and retain limma-voom gene results as exploratory. The primary snCv3 matrix remains verified integer RNA counts. This representation clarification does not change cohort membership, fixed programs, thresholds, contrasts or multiplicity families. It is not evidence of independent obstruction replication.

The GEO source GSE169285 and atlas metadata identify reference biopsy libraries. Three PRE019 libraries share a source patient ID and are merged. EO-suffixed IDs are distinct patient entries in the atlas; common number prefixes are not collapsed (some carry different sex annotations). Fourteen patient IDs shared with the full snCv3 metadata are excluded. Describe the result as a second assay from the same atlas with shared reported IDs excluded, not as a separate independent study; unrecorded aliases cannot be ruled out from public identifiers alone.

