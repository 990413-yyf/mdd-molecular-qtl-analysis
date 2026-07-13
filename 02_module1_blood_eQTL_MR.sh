#!/bin/bash
## ============================================================================
## 02_module1_blood_eQTL_MR.sh  — 血 eQTL-MR SMR+HEIDI 全基因组扫描
## 与 02_module1_blood_eQTL_MR.R (annotation+coloc+tier) 配套
## ----------------------------------------------------------------------------
## 输入:
##   data/eqtl/cis-eQTL-SMR_20191212  (eQTLGen 全血 cis-eQTL BESD 三件套 .besd/.esi/.epi;
##                                      19250 gene probes, 8932843 SNP, N≈31684)
##   data/gwas/mdd2025.ma             (MDD2025 no23andMe-noUKBB EUR daner→.ma, hg19)
##   data/ld_ref/1000G.EUR.QC         (EUR 1000G PLINK bfile, 489 EUR, GRCh37)
## 输出:
##   out/smr_eqtlgen_mdd2025.smr      (15339 探针 SMR+HEIDI 结果)
## 对应正文: 模块一 血液因果基因 (SMR+HEIDI 主 + coloc 交叉验证)
## ----------------------------------------------------------------------------
## 工具: SMR 1.03 macOS-x86_64 经 Rosetta 2 (arch -x86_64)。
##   注: SMR 1.4.1 macOS-arm64 官方二进制读 .esi 段错误 (arm64 build bug, 自建
##       2-SNP BESD 亦崩), 故降级 1.03 x86_64 经 Rosetta 跑。
## 样本重叠自查: eQTLGen × MDD2025(no-UKB) ≈ 0。
## 参数: cis 窗口 ±1Mb; --peqtl-smr 5e-8; HEIDI 保留 p_HEIDI>1.57e-3;
##   --diff-freq-prop 0.2。
## 注: CJK 路径会使 R system()/plink 报 "cannot open the connection", 故所有输入
##   先复制到 workspace ASCII 路径 (data/, out/) 再跑。
## ============================================================================
set -e
SMR="tools/smr_x86/smr_Mac"                 # SMR 1.03 x86_64
BESD="data/eqtl/cis-eQTL-SMR_20191212"      # BESD prefix (.besd/.esi/.epi)
GWAS="data/gwas/mdd2025.ma"
LDREF="data/ld_ref/1000G.EUR.QC"            # PLINK bfile prefix
OUT="out/smr_eqtlgen_mdd2025"

arch -x86_64 "$SMR" \
    --bfile "$LDREF" \
    --gwas-summary "$GWAS" \
    --beqtl-summary "$BESD" \
    --peqtl-smr 5e-8 \
    --heidi-mtd 1 \
    --peqtl-heidi 1.57e-3 \
    --cis-wind 1000 \
    --diff-freq-prop 0.2 \
    --thread-num 6 \
    --out "$OUT"
echo "SMR done -> $OUT.smr"
## 结果: 4409879 SNP 过等位检查, 15339 探针; p_SMR<0.05=2586, FDR<0.05=721,
##        FDR<0.05 且 HEIDI-pass=513。
