# Independent eQTL replication: data manifest

This stage re-runs SMR and HEIDI for the 142 tier-1 genes, with emphasis on the 15 tier-A genes, using eQTL resources distinct from the discovery datasets. The MDD outcome GWAS is unchanged, so this is an exposure-resource robustness analysis rather than independent outcome replication.

## Required: GTEx v8 cis-eQTL SMR BESD

- Source: Yang Lab SMR data portal.
- Resource: GTEx v8 cis-eQTL summary, lite release (all tissues).
- URL: `https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary_lite.tar`
- Approximate size: 506 MB.
- Default repository location after extraction: `data/raw/gtex_v8/`.
- Tissues used: `Whole_Blood`, `Brain_Cortex`, and `Brain_Frontal_Cortex_BA9`.

Example:

```bash
mkdir -p data/raw/gtex_v8
cd data/raw/gtex_v8
curl -L -O "https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary_lite.tar"
tar xf GTEx_V8_cis_eqtl_summary_lite.tar
```

Set `GTEX_DIR` if the resource is stored elsewhere.

## Optional: PsychENCODE prefrontal cortex cis-eQTL

- URL: `https://yanglab.westlake.edu.cn/data/SMR/PsychENCODE_cis_eqtl_PEER50_summary.tar.gz`
- Suggested location: `data/raw/psychencode/`.

## Coordinate convention

The listed SMR resources, the MDD `.ma` file and the 1000 Genomes European LD panel used in this project are in hg19/b37 coordinates. No liftover is needed for this replication stage.

## Output and interpretation

Run `10_gtex_replication_SMR.sh` after configuring the paths. A gene is counted as replicated when it is testable in the matched tissue, has BH-FDR <0.05, has HEIDI P >0.01 and agrees in direction with discovery. Missing probes are reported as not testable, not as negative results.
