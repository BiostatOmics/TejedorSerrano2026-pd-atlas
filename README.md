# TejedorSerrano2026-pd-atlas

[![Preprint](https://img.shields.io/badge/Preprint-10.21203%2Frs.3.rs--9220144%2Fv1-blue)](https://doi.org/10.21203/rs.3.rs-9220144/v1)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

Code and data accompanying the manuscript:

> **Decoding the multiregional atlas of Parkinson's disease at single-cell resolution**
> Tejedor-Serrano JA, Arzalluz-Luque Á, Tarazona S\*, Martorell-Marugán J\*, 2026

## Repository structure

```
00-envs/            Conda environment files for reproducibility
01-preprocessing/   Raw data download and cell QC/clustering
02-analysis/        Differential abundance (tascCODA) and DEA/GSEA
03-data/            Metadata and gene-level reference files (AMP PD access required)
04-analysis/        Pre-computed analysis results (AMP PD access required)
05-figures/         R scripts that generate all manuscript figures
utils/              Shared R helper functions
```

## Environments

All environment files are in `00-envs/`. Recreate with `conda env create -f <file>`.

| File | Language | Used for |
|------|----------|----------|
| `r-sc-full.yml` | R | All figure scripts (`05-figures/`) and DEA/GSEA analyses |
| `py-sc-full.yml` | Python | Preprocessing and general single-cell analysis |
| `pertpy-full.yml` | Python | Compositional differential abundance (`02-analysis/02-*`) |

The `-history.yml` variants list only explicitly installed packages (more readable; less strictly pinned).

## Data access

> For full details on data processing, analysis design, and interpretation please refer to the manuscript (preprint linked above).

The data used in this study are from the **Accelerating Medicines Partnership® Parkinson's Disease (AMP® PD)** controlled-access dataset and cannot be redistributed. To reproduce the analyses, register and apply for access at [amp-pdrd.org](https://www.amp-pdrd.org), then re-run the preprocessing and analysis scripts in `01-preprocessing/` and `02-analysis/`.

## Acknowledgements

Data used in the preparation of this article were obtained from the Accelerating Medicine Partnership® (AMP®) Parkinson's Disease (AMP PD) and Parkinson's Disease & Related Disorders (AMP PDRD) Knowledge Platform. For up-to-date information on the study, visit https://www.amp-pdrd.org.

The AMP® PD program is a public-private partnership managed by the Foundation for the National Institutes of Health and funded by the National Institute of Neurological Disorders and Stroke (NINDS) in partnership with the Food and Drug Administration (FDA), National Institute on Aging (NIA), Aligning Science Across Parkinson's (ASAP) initiative; Celgene Corporation, a subsidiary of Bristol-Myers Squibb Company; GlaxoSmithKline plc (GSK); The Michael J. Fox Foundation for Parkinson's Research (MJFF); AbbVie Inc.; Pfizer Inc.; Sanofi US Services Inc.; and Verily Life Sciences LLC.

The AMP® PDRD program is a public-private partnership managed by the Foundation for the National Institutes of Health and funded by the National Institute of Neurological Disorders and Stroke (NINDS) in partnership with the Food and Drug Administration (FDA), National Institute on Aging (NIA), AbbVie Inc.; Aligning Science Across Parkinson's (ASAP) initiative; C2N Diagnostics, LLC; CurePSP; GlaxoSmithKline plc (GSK); Denali Therapeutics Inc.; Laboratory Corporation of America Holdings (Labcorp); The Michael J. Fox Foundation for Parkinson's Research (MJFF); Sanofi US Services Inc.; and Verily Life Sciences LLC.

ACCELERATING MEDICINES PARTNERSHIP and AMP are registered service marks of the U.S. Department of Health and Human Services.

Clinical data and biosamples used in preparation of this article were obtained from the (i) Michael J. Fox Foundation for Parkinson's Research (MJFF) and National Institutes of Neurological Disorders and Stroke (NINDS) BioFIND study, (ii) Harvard Biomarkers Study (HBS) and the Stephen & Denise Adams Center for Parkinson's Disease Research of Yale School of Medicine (CPDR-Y), (iii) National Institute on Aging (NIA) International Lewy Body Dementia Genetics Consortium Genome Sequencing in Lewy Body Dementia Case-control Cohort (LBD), (iv) MJFF LRRK2 Cohort Consortium (LCC), (v) NINDS Parkinson's Disease Biomarkers Program (PDBP), (vi) MJFF Parkinson's Progression Markers Initiative (PPMI), (vii) NINDS Study of Isradipine as a Disease-modifying Agent in Subjects With Early Parkinson Disease, Phase 3 (STEADY-PD3), and (viii) the NINDS Study of Urate Elevation in Parkinson's Disease, Phase 3 (SURE-PD3). Full cohort acknowledgement text is provided in the manuscript.

Data used in the preparation of this article were also obtained from the Global Parkinson's Genetics Program (GP2). GP2 is funded by the Aligning Science Against Parkinson's (ASAP) Initiative and implemented by The Michael J. Fox Foundation for Parkinson's Research (www.gp2.org).
