#!/bin/bash
# 模块一 血 eQTL-MR: SMR+HEIDI 全基因组扫描
# exposure=eQTLGen cis-eQTL (BESD), outcome=MDD2025 daner, LD ref=EUR 1000G
set -e
PROJ="${PROJECT_ROOT:-$PWD}"
BESD="$PROJ/data/eqtl/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
GWAS="$PROJ/data/gwas/mdd2025.ma"          # SMR .ma format (built from daner)
LDREF="${LDREF_PREFIX:-$PROJ/data/ld_ref/EUR}"
OUT="$PROJ/out/smr_eqtlgen_mdd2025"
smr --bfile "$LDREF" \
    --gwas-summary "$GWAS" \
    --beqtl-summary "$BESD" \
    --peqtl-smr 5e-8 \
    --heidi-mtd 1 \
    --peqtl-heidi 1.57e-3 \
    --cis-wind 1000 \
    --thread-num 6 \
    --diff-freq-prop 0.2 \
    --out "$OUT"
echo "SMR done -> $OUT.smr"
