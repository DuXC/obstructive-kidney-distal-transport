# External validation figure contract

Core claim: the fixed distal transport estimates have context-dependent uncertainty across public injury datasets. Figures report the entire frozen comparison family and do not select a favorable sensitivity result.

Authoritative generator: Python matplotlib, continuing the existing project stack. Source data: sn_program_comparisons.tsv, sc_program_comparisons.tsv (only after count-representation acceptance), sn_donor_scores.tsv, and sc_donor_scores.tsv.

Fig 4: quantitative grid with AKI and CKD primary snRNA estimates. Every fixed program is shown with donor n, 95% Welch interval and the ten-test FDR. Include corresponding donor score panels to show small-reference-group variation.

Fig S5: secondary scRNA assay estimates and donor scores, subject to representation and eligibility gates. Describe source preprocessing and shared-patient exclusion.

Fig S6: snRNA sensitivity effects across thresholds, anatomy, reference-state restriction, COVID exclusion and sex adjustment. Fixed rows include uncertain outcomes. Panel claims focus on interpretation stability; the primary family remains authoritative.

Export: 183 mm width, white background, editable SVG/PDF text, 300 dpi PNG review and 600 dpi compressed TIFF. Sans serif labels, minimum 7 pt, blue/orange/neutral palette consistent with the existing manuscript. All axes, intervals, n, denominators and multiplicity definitions supplied in legends. Inspect every standalone image and final embedded page.
