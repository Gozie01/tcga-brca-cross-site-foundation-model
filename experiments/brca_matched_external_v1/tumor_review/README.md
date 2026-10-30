# Next step: outcome-blinded tumor review

No prediction experiment or feature extraction is required to prepare this review.

1. Locate the previously downloaded brca_stage2_full_qc_20261021_041052_426646.zip. Place a copy beside prepare_review.py.
2. Open PowerShell in this tumor_review folder and run:

   py -3.11 .\prepare_review.py --qc-zip .\brca_stage2_full_qc_20261021_041052_426646.zip

   If the file is elsewhere, provide its full quoted path after --qc-zip. Python standard library only; no installation or GPU is required.
3. Give ONLY the generated reviewer_pack folder and these instructions to a qualified histology reviewer. Do not share the parent result package, predictions, labels, or review_linkage_PRIVATE.csv with the reviewer.
4. The reviewer opens review_index.html, inspects each contact sheet and records findings in review_decisions.csv. Save reviewer name, qualification, date, review scope and reasons. The study custodian retains the private linkage for locating the original WSI when needed.
5. Return the completed review_decisions.csv for reconciliation before any eligibility change or recomputation.

## What to record

sampled_tumor_present: YES / NO / UNCERTAIN. This applies only to the displayed sampled patches.
review_scope: CONTACT_SHEET / WHOLE_SLIDE / CONTACT_SHEET_AND_WHOLE_SLIDE.
slide_decision: RETAIN / EXCLUDE / DEFER, with a reason. Start every row as PENDING.

A contact sheet samples only some retained tissue patches. Tumor not seen on that sheet does not prove tumor absence on the full slide. Use DEFER and request WSI review when uncertain; do not automatically exclude negative previews. Retaining a slide does not establish that every analysed patch is tumor or that it is the molecularly assayed block. Record these distinctions in the manuscript.

Review IDs are randomly ordered independently of molecular labels. The generated pack contains no patient IDs, molecular labels or model predictions. Existing review_decisions.csv is never overwritten when preparation is repeated.
