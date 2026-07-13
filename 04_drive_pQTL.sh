#!/bin/bash
## ============================================================================
## 04_drive_pQTL.sh  — 逐蛋白驱动: awk 提取 cis 区 ±1Mb → 调 04_pQTL_inflammation_MR.R
## 用法: bash 04_drive_pQTL.sh <gene> <GCST> <chr> <gene_start> <gene_end>
## 对 91 个蛋白循环调用 (GCST→蛋白映射见 out/zhao2023_gcst_protein_map.csv,
##   基因 GRCh38 坐标见 out/zhao2023_gene_coords_grch38.json)。
## ============================================================================
# drive_pqtl.sh <gene> <GCST> <chr> <gene_start> <gene_end>
# extracts cis region (±1Mb) via awk then runs R MR
gene="$1"; gcst="$2"; chrom="$3"; gstart="$4"; gend="$5"
DL="${PQTL_DIR:-data/raw/scallop}"
CIS_DIR="${PQTL_CIS_DIR:-data/pqtl/cis}"
CIS="$CIS_DIR/${gene}.txt"
win=1000000
lo=$(( ${gstart%.*} - win )); [ $lo -lt 0 ] && lo=0
hi=$(( ${gend%.*} + win ))
mkdir -p "$CIS_DIR"
# stream-extract cis window (deCODE/Zhao cols: chr=1 pos=2 ... p=8)
gzip -dc "$DL/${gcst}.tsv.gz" | awk -F'\t' -v c="$chrom" -v lo="$lo" -v hi="$hi" \
  'NR==1{print;next} $1==c && $2>=lo && $2<=hi' > "$CIS"
nrow=$(($(wc -l < "$CIS")-1))
echo "[$gene] cis SNPs extracted: $nrow (chr$chrom:$lo-$hi)"
Rscript 04_pQTL_inflammation_MR.R "$gene" "$gcst" "$chrom" "$gstart" "$gend" 2>&1 | \
  grep -vE "Loading|Attaching|masked|The following|^\s*$|Warning message|was built under|^Note:|^ *NOTE" | tail -30
