# Stage 12 Results Synthesis v1

This frozen package converts the verified Stage 10–11 outputs into manuscript-
ready numerical evidence.

## Contents

- `RESULTS_SECTION_DRAFT.md`: evidence-grounded Results narrative, working title,
  central claim, and interpretation boundary.
- `figures/`: four publication figures in vector PDF and 400-dpi PNG formats.
- `tables/`: five publication tables covering cohorts, internal performance,
  site-held-out performance, representation ablations, and statistically
  supported paired effects.
- `stage12_results_synthesis_metadata_v1.json`: synthesis status and scope.
- `build_stage12_results_synthesis.py`: reproducible figure/table generator.

## Required interpretation boundary

The A2, C8, E2, and BH experiments are institution-held-out evaluations within
TCGA-BRCA. They are not independent external-cohort validation. Claims should
therefore concern cross-site robustness within TCGA-BRCA, not clinical
generalization to unrelated cohorts or prospective deployment.
