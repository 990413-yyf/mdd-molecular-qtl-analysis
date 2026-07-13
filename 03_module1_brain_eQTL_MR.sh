#!/bin/bash
## ============================================================================
## 03_module1_brain_eQTL_MR.sh  — 脑 eQTL-MR SMR+HEIDI (逐染色体)
## 与 03_module1_brain_eQTL_MR.R (annotation+coloc+脑血比较) 配套
## ----------------------------------------------------------------------------
## 输入:
##   data/brain_eqtl/chr{1..22}.{besd,esi,epi}  (BrainMeta v2 皮层 cis-eQTL,
##       解压自 BrainMeta_cis_eqtl_summary.tar.gz; 16744 探针, 11630972 SNP)
##   data/gwas/mdd2025.ma        (复用模块一, 不重下)
##   data/ld_ref/1000G.EUR.QC    (复用模块一 EUR 1000G)
## 输出:
##   out/smr_brain_mdd2025.smr   (22 条染色体 concat, 15579 探针)
## 对应正文: 模块一 脑 tier-1 因果基因 + 脑血比较
## ----------------------------------------------------------------------------
## 工具: SMR 1.03 x86_64 经 Rosetta 2 (同血, arm64 build bug)。
## 注: --besd-flist --make-besd 合并在 CJK 路径下报 "Column number not correct" 失败,
##     故改为逐染色体跑 SMR+HEIDI 再 concat .smr (BESD 本已按染色体拆分)。
## 样本重叠自查: BrainMeta(非UKB) × MDD2025(无UKB) ≈ 0。
## 注: BrainMeta v2 皮层文献(Qi et al. 2022)与样本量 n≈2865 为记忆, 未联网核对,
##     .summary 文件不含样本量字段 (投稿前须核实)。
## ============================================================================
set -e
SMR="tools/smr_x86/smr_Mac"
BDIR="data/brain_eqtl"
GWAS="data/gwas/mdd2025.ma"
LDREF="data/ld_ref/1000G.EUR.QC"
mkdir -p out/brain_perchr
for chr in $(seq 1 22); do
  arch -x86_64 "$SMR" \
    --bfile "$LDREF" \
    --gwas-summary "$GWAS" \
    --beqtl-summary "$BDIR/chr${chr}" \
    --peqtl-smr 5e-8 --heidi-mtd 1 --peqtl-heidi 1.57e-3 \
    --cis-wind 1000 --diff-freq-prop 0.2 --thread-num 6 \
    --out "out/brain_perchr/chr${chr}"
  echo "chr${chr} done"
done
## concat: 保留一份表头 + 各染色体数据行
head -1 out/brain_perchr/chr1.smr > out/smr_brain_mdd2025.smr
for chr in $(seq 1 22); do tail -n +2 "out/brain_perchr/chr${chr}.smr"; done >> out/smr_brain_mdd2025.smr
echo "brain SMR concat -> out/smr_brain_mdd2025.smr ($(($(wc -l < out/smr_brain_mdd2025.smr)-1)) probes)"
## 结果: 15579 探针; p_SMR<0.05=2639, FDR<0.05=569, HEIDI-pass=449。
