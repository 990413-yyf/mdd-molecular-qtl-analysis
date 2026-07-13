## ============================================================================
## 02_module1_blood_eQTL_MR.R  — 血 SMR 结果注释 + coloc 交叉验证 + tier 表
## 前置: 先跑 02_module1_blood_eQTL_MR.sh 得 out/smr_eqtlgen_mdd2025.smr
## ----------------------------------------------------------------------------
## 输入:
##   out/smr_eqtlgen_mdd2025.smr      (SMR+HEIDI 结果)
##   out/mdd2025.ma                   (GWAS, coloc 用)
##   out/tier1_cis_eqtl.txt           (513 个 tier1 探针的 cis-eQTL 记录; 由 SMR
##                                      --query 或 besd 提取)
##   ensg_symbol_map.tsv              (ENSG→HGNC, 由 eQTLGen sig 文件构建)
## 输出:
##   out/blood_causal_gene_tier_table.csv   (全 15339 探针分层)
##   out/blood_causal_tier1_dual.csv        (80 个 SMR+HEIDI+coloc 三重通过)
##   out/coloc_tier1.rds                    (513 探针 coloc 结果, checkpoint)
##   out/fig_smr_volcano.png / fig_forest_top20.png (见 python 段)
## 对应正文: 模块一 血液 tier-1 因果基因表 + 火山/森林图
## ----------------------------------------------------------------------------
## 双方法: SMR FDR<0.05 & HEIDI-pass(p_HEIDI>0.01) & coloc PP.H4>0.8 = tier1_dual。
## 阳性对照 IL6R/CRP/TNF/P2RX7 全阴性 → 无经典炎症基因进入因果层。
## ============================================================================
# 运行前将工作目录设为 workspace 根 (含 out/、data/ 的 ASCII 路径), 例如:
# setwd("/path/to/workspace")   # 切勿用含中文的路径 (R read/system 会 "cannot open the connection")
suppressMessages({library(data.table); library(coloc)})

## ---------------- 1. 注释 SMR 结果 (ENSG→symbol, FDR, tier) ----------------
smr <- fread("out/smr_eqtlgen_mdd2025.smr")
map <- fread("ensg_symbol_map.tsv", header = TRUE)               # ENSG→GeneSymbol
smr[, Gene := map$GeneSymbol[match(probeID, map$ENSG)]]
smr[, FDR := p.adjust(p_SMR, "BH")]
smr[, heidi_pass := (!is.na(p_HEIDI) & p_HEIDI > 0.01)]
smr[, tier := fifelse(FDR < 0.05 & heidi_pass, "tier1",
              fifelse(FDR < 0.05 & !is.na(p_HEIDI) & !heidi_pass, "tier3_HEIDI_fail",
              fifelse(FDR < 0.05 & is.na(p_HEIDI), "tier4_HEIDI_NA", "ns")))]
saveRDS(smr, "out/smr_annotated.rds")
## 结果: tier1(HEIDI-pass)=513, tier3_HEIDI_fail=207, tier4_HEIDI_NA=1(LRP1)。

## ---------------- 2. coloc.abf 交叉验证 (513 tier1 探针) ----------------
gw <- fread("out/mdd2025.ma")
setnames(gw, c("SNP", "A1_g", "A2_g", "freq_g", "b_g", "se_g", "p_g", "n_g"))
gw <- gw[!is.na(b_g) & !is.na(se_g) & se_g > 0]; setkey(gw, SNP)
eq <- fread("out/tier1_cis_eqtl.txt")     # cols: SNP Chr BP A1 A2 Freq Probe ... b SE p
s_frac <- 357636 / (357636 + 1281936)     # MDD2025 case 比例
N_eqtl <- 31684
probes <- unique(eq$Probe); res <- vector("list", length(probes))
for (i in seq_along(probes)) {
  pr <- probes[i]; e <- eq[Probe == pr]; e <- e[!is.na(b) & !is.na(SE) & SE > 0]
  m <- merge(e, gw, by = "SNP")
  if (nrow(m) < 10) { res[[i]] <- data.table(probe = pr, nsnp = nrow(m), PP.H4 = NA_real_, PP.H3 = NA_real_); next }
  D_eqtl <- list(beta = m$b,   varbeta = m$SE^2,   snp = m$SNP, type = "quant", N = N_eqtl,        MAF = pmin(m$Freq, 1 - m$Freq))
  D_gwas <- list(beta = m$b_g, varbeta = m$se_g^2, snp = m$SNP, type = "cc", s = s_frac, N = round(m$n_g), MAF = pmin(m$freq_g, 1 - m$freq_g))
  ct <- tryCatch(coloc.abf(D_eqtl, D_gwas), error = function(e) NULL)
  if (is.null(ct)) { res[[i]] <- data.table(probe = pr, nsnp = nrow(m), PP.H4 = NA_real_, PP.H3 = NA_real_); next }
  pp <- ct$summary
  res[[i]] <- data.table(probe = pr, nsnp = nrow(m), PP.H4 = pp["PP.H4.abf"], PP.H3 = pp["PP.H3.abf"])
}
colocdt <- rbindlist(res); saveRDS(colocdt, "out/coloc_tier1.rds")
## 结果: 513/513 有结果; PP.H4>0.8 = 80, >0.5 = 163。

## ---------------- 3. tier1_dual 表 (三重通过) ----------------
smr <- merge(smr, colocdt[, .(probeID = probe, coloc_nsnp = nsnp, coloc_PPH4 = PP.H4)],
             by = "probeID", all.x = TRUE)
smr[, dual_pass := (tier == "tier1" & !is.na(coloc_PPH4) & coloc_PPH4 > 0.8)]
smr[, final_tier := fifelse(dual_pass, "tier1_dual",
                    fifelse(tier == "tier1", "tier2_SMRonly",
                    fifelse(tier == "tier3_HEIDI_fail", "tier3_HEIDIfail",
                    fifelse(tier == "tier4_HEIDI_NA", "tier4_HEIDI_NA", "ns"))))]
smr[, direction := fifelse(b_SMR > 0, "risk_up", "protective_down")]
fwrite(smr, "out/blood_causal_gene_tier_table.csv")
fwrite(smr[final_tier == "tier1_dual"], "out/blood_causal_tier1_dual.csv")
## 结果: tier1_dual=80 (risk↑50/protective↓30); 无经典炎症基因。
## 阳性对照: IL6R b=-0.0009 p=0.978, TNF b=0.008 p=0.481, P2RX7 b=0.014 p=0.395,
##           CRP 缺失 → 全阴性。

## ---------------- 火山/森林图 (python 段, 示意) ----------------
## fig_smr_volcano.png: x=b_SMR, y=-log10(p_SMR), tier1_dual 高亮, FDR<0.05 阈值线。
## fig_forest_top20.png: top20 tier1_dual 森林图 (b_SMR ± 1.96*se_SMR)。
