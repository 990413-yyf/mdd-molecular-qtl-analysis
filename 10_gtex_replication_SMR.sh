#!/usr/bin/env bash
# ============================================================
# tier-1 / tier-A 因果基因在独立 eQTL(GTEx v8)中的 SMR 复现
# 坐标 hg19/b37,与我们 MDD .ma 和 1000G 面板一致,无需 liftover。
# 目的:对发现集(血 eQTLGen / 脑 BrainMeta)锁定的 142 tier-1(重点 15 tier-A),
#       在 GTEx 全血 / 脑皮层重跑 SMR+HEIDI,核对显著性与方向一致性→复现率。
# ============================================================
set -euo pipefail

# ---- Configure with environment variables when paths differ ----
SMR="${SMR_BIN:-smr}"
MDD_MA="${MDD_MA:-data/gwas/mdd2025.ma}"
BFILE="${LDREF_PREFIX:-data/ld_ref/1000G.EUR.QC}"
GTEX="${GTEX_DIR:-data/raw/gtex_v8}"          # directory containing unpacked GTEx BESD resources
OUT="${GTEX_OUT:-out/replication_gtex}"; mkdir -p "$OUT"

# ---- 各组织跑 SMR(cis,±1Mb,peqtl 5e-8,HEIDI 开) ----
# GTEx besd 组织名(解压后):
for TISS in Whole_Blood Brain_Cortex Brain_Frontal_Cortex_BA9; do
  echo ">>> SMR replication in GTEx $TISS"
  "$SMR" --bfile "$BFILE" \
         --gwas-summary "$MDD_MA" \
         --beqtl-summary "$GTEX/$TISS" \
         --peqtl-smr 5e-8 --cis-wind 1000 --heidi-mtd 1 \
         --thread-num 8 \
         --out "$OUT/GTEx_${TISS}_SMR"
done

# ---- 汇总:把 tier-1/tier-A 对到 GTEx SMR 结果,判复现 ----
# 复现判据(逐基因):在该组织 GTEx 有该 probe、SMR 达显著(先报名义 p<0.05 与 BH-FDR<0.05 两档)、
#   通过 HEIDI(p_HEIDI>0.01)、且 b_SMR 方向与发现集一致。
# Merge $OUT/*.smr (probe/Gene/b_SMR/p_SMR/p_HEIDI) with the discovery table and
# stable symbol-to-Ensembl mapping to produce:
#   out/replication_gtex.csv(gene,tier,tissue,in_gtex,b_SMR,p_SMR,p_HEIDI,dir_concordant,replicated_nominal,replicated_fdr)
#   并汇总:15 个 tier-A 中在 GTEx 血/脑至少一处复现(方向一致+显著)的个数、142 tier-1 的复现率。
# 诚实框定:GTEx 样本量小于 eQTLGen/BrainMeta,部分基因可能无强 cis 工具或不复现——
#   如实报"N/15 复现",不硬凑;未复现者列出并说明(多因 GTEx 该组织无显著 cis-eQTL)。
echo ">>> SMR stage complete. Summarize against the locked discovery-gene list and report missing probes explicitly."
