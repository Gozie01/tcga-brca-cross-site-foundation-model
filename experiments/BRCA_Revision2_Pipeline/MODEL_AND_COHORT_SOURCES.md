# Encoder and external-cohort sources

Checked 19 October 2026 against primary publisher/maintainer sources.

| Item | Source | Consequence |
| --- | --- | --- |
| Hibou-B | https://huggingface.co/histai/hibou-b | Publisher model card identifies an 85.7M-parameter pathology encoder and Apache-2.0 licence. Model files require a Hugging Face access agreement. |
| Official model loading and output | https://github.com/HistAI/hibou/blob/main/example.ipynb | Use the publisher image processor and `pooler_output`; record exact downloaded revision and file hashes. |
| Hibou paper | https://arxiv.org/abs/2406.05074 | Pathology pretraining on a proprietary collection; do not assert complete pretraining/evaluation independence without an overlap audit. |
| TCIA CPTAC-BRCA collection | https://www.cancerimagingarchive.net/collection/cptac-brca/ | Registry lists 134 subjects and about 113 GB of pathology images. These are registry counts, not our eligible validation cohort. |
| cBioPortal CPTAC 2020 study | https://www.cbioportal.org/study/summary?id=brca_cptac_2020 | A molecular/clinical cohort candidate; imaging overlap must be established. |
| Curated data provenance | https://github.com/cBioPortal/datahub/blob/master/public/brca_cptac_2020/README.md | Clinical and molecular records were curated from Krug et al. 2020 supplements. Review assay fields and mutation profiling before assigning labels. |
| Historical CPTAC-TCGA distinction | https://www.cancerimagingarchive.net/research/cptac-imaging-proteomics/ | Some earlier CPTAC resources include TCGA-derived cases. A CPTAC name alone is not proof of independent validation. |

No model access conditions were accepted on the user's behalf. The public clinical metadata API probe timed out in this workspace; no endpoint labels or external-validation results were inferred from search snippets.
