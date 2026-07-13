# Multi-tissue molecular QTL analyses in major depressive disorder

This repository contains the analysis code accompanying the manuscript **“Multi-tissue molecular QTL analyses prioritize candidate genes and constrain causal claims for circulating inflammation in major depressive disorder.”**

The workflow separates three evidence layers:

1. whole-blood and cortical cis-eQTL analyses using SMR, HEIDI and colocalization;
2. cis-pQTL Mendelian-randomization analyses of 91 circulating inflammatory proteins in SCALLOP and UKB-PPP;
3. peripheral-blood and exploratory acupuncture-related expression analyses.

The repository contains code and metadata, not redistributed third-party molecular-QTL or GWAS summary statistics. Obtain each dataset under its source terms and place it under `data/` as described below.

## Main findings reproduced by the workflow

- 80 blood and 68 cortical tier-1 genes, yielding 142 unique candidates after tissue integration.
- 15 candidates in the prespecified tier-A evidence group.
- 62 of 72 testable candidates meeting the independent-eQTL-resource robustness criteria.
- 74 of 91 circulating inflammatory proteins testable across SCALLOP and UKB-PPP; none meeting both panel-corrected MR and PP.H4 >0.80.
- Model-dependent and heterogeneous peripheral-blood expression findings.
- Exploratory manual-acupuncture expression comparisons that are not treated as causal or clinical evidence.

## Repository map

| Stage | Script | Purpose |
|---|---|---|
| 01 | `01_acupuncture_reversal.R` | Within-study and cross-species direction comparisons for GSE86392 |
| 02 | `02_module1_blood_eQTL_MR.sh`, `02_module1_blood_eQTL_MR.R` | Whole-blood eQTLGen SMR/HEIDI and colocalization |
| 03 | `03_module1_brain_eQTL_MR.sh`, `03_module1_brain_eQTL_MR.R` | BrainMeta cortical SMR/HEIDI, colocalization and tissue integration |
| 04 | `04_drive_pQTL.sh`, `04_pQTL_inflammation_MR.R` | SCALLOP cis-pQTL instrument extraction, MR and colocalization |
| 05 | `05_module2_transcriptome_meta_WGCNA_immune.R`, `05b_immune_deconvolution.py` | Four-cohort blood meta-analysis, WGCNA and marker-score analysis |
| 06 | `06_module2_GO_KEGG_enrichment.R` | Exploratory GO and KEGG enrichment |
| 07 | `07_module5_network_STRING_OpenTargets_tiering.py` | STRING/Open Targets annotation and transparent tier scoring |
| 08 | `08_build_workflow_figure.py` | Vector workflow figure |
| 09 | `09_download_ukbpp_user.sh`, `09b_download_ukbpp_missing.sh` | UKB-PPP archive retrieval through Synapse |
| 10 | `10_gtex_replication_SMR.sh`, `10_gtex_download_manifest.md` | GTEx/PsychENCODE eQTL-resource robustness stage |

`run_smr_module1.sh` is a compact whole-blood SMR launcher retained for provenance. The incomplete early prescan is isolated under `legacy/` and is not part of the reported workflow.

## Expected directory structure

```text
data/
  gwas/mdd2025.ma
  ld_ref/1000G.EUR.QC.{bed,bim,fam}
  eqtl/
  raw/
    scallop/
    ukbpp/
    gtex_v8/
metadata/
out/
```

Large inputs and generated outputs are ignored by Git. Stable accession maps and small non-restricted metadata should be placed in `metadata/`.

## Data sources

- MDD outcome: Psychiatric Genomics Consortium MDD2025 European summary statistics excluding 23andMe and UK Biobank; 357,636 cases and 1,281,936 controls.
- Whole-blood eQTL: eQTLGen cis-eQTL summary data.
- Cortical eQTL: BrainMeta v2 cis-eQTL summary data.
- Replication eQTL: GTEx v8 and PsychENCODE SMR BESD resources from the Yang Lab portal.
- SCALLOP inflammatory-protein pQTL: NHGRI-EBI GWAS Catalog GCST accessions listed in the manuscript Supplementary Table S17.
- UKB-PPP: European discovery pGWAS summary statistics in Synapse collection `syn51365303`; per-protein identifiers are listed in Supplementary Table S3.
- Expression datasets: GEO `GSE86392`, `GSE98793`, `GSE76826`, `GSE38206`, `GSE39653`, `GSE102556` and `GSE124387`.
- LD reference: 1000 Genomes phase 3 European panel.

Users are responsible for following the access and redistribution conditions of each source dataset.

## Configuration

Run scripts from the repository root. Defaults are repository relative. Override locations with environment variables when needed:

```bash
export PROJECT_ROOT="$PWD"
export MDD_MA="data/gwas/mdd2025.ma"
export LDREF_PREFIX="data/ld_ref/1000G.EUR.QC"
export PQTL_DIR="data/raw/scallop"
export PQTL_CIS_DIR="data/pqtl/cis"
export PQTL_OUT="out/pqtl_inflammation_MR.csv"
export UKBPPP_DEST="data/raw/ukbpp"
export GTEX_DIR="data/raw/gtex_v8"
```

## Software

The analysis environment used R 4.5.3. Principal R packages included `data.table`, `TwoSampleMR`, `coloc`, `ieugwasr`, `GEOquery`, `limma`, `metafor`, `clusterProfiler`, `WGCNA` and `biomaRt`. Principal Python packages included `pandas`, `numpy`, `scipy`, `matplotlib` and `networkx`. See `00_sessionInfo.txt` for the recorded R session.

SMR 1.03 x86_64 was used through Rosetta on Apple Silicon for the reported run. PLINK was used for LD clumping.

## Minimal validation

Syntax-only checks do not require the restricted datasets:

```bash
for f in *.R; do Rscript -e 'parse(file=commandArgs(TRUE)[1])' "$f"; done
for f in *.sh; do bash -n "$f"; done
python3 -m py_compile *.py
```

Full numerical reproduction requires the source datasets, access credentials where applicable and the exact accession/probe maps supplied with the submission.

## Interpretation boundaries

“No robust evidence” refers only to proteins that were instrumentable under the prespecified strong cis-instrument, multiple-testing and colocalization criteria. It does not exclude smaller effects, uninstrumented proteins, brain-local or stimulus-dependent inflammation. The acupuncture analysis is exploratory and does not establish treatment mechanism or efficacy.

## Citation, permanent archive and licence

The manuscript citation, repository URL, Zenodo DOI and code licence will be added before public release.
