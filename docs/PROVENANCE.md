# Public-release provenance

The numerical outputs in this release are copied from the locked Stage 10–11 analysis artifacts. Output-file SHA-256 values embedded in the metadata were verified after packaging.

To avoid publishing author-specific Windows directory information, absolute paths inside four JSON metadata files were replaced with `<PROJECT_ROOT>`. No numerical value, validation gate, model specification, sample count, seed, output hash, or statistical result was modified. The SHA-256 values of the original locked metadata files before path sanitization were:

| Metadata file | Original SHA-256 |
|---|---|
| `stage10a_cohort_split_freeze_metadata_v1.json` | `ec0f61903bf84c88c8edb256e2b886ce01c605817215fdb91b6dc6686e733d9d` |
| `stage10b_internal_prediction_calibration_metadata_v1.json` | `4ab217522a5ee22485fe0069efc0c7ad7b4a206dd23764dae1ea4774c98d7747` |
| `stage10c_site_holdout_metadata_v1.json` | `9f523a65a5474295429da859517ebcdae792cec95dba514e213ad9c511b2663c` |
| `stage11_statistical_analysis_ablations_metadata_v1.json` | `4f754a44d16961e1480856ea858855bcc27fef317a6c81c06095b9ebfb141ac1` |

The public Stage 9–11 scripts likewise differ from the production copies only in path portability and terminology clarifying that A2, C8, E2, and BH are site holdouts within TCGA-BRCA rather than independent external cohorts. The locked algorithms, thresholds, seeds, model configurations, and expected cohort counts are unchanged.
