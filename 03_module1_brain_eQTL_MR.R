## ============================================================================
## 03_module1_brain_eQTL_MR.R  — 脑 SMR 注释 + coloc + 脑血比较 + 统一因果表
## 前置: 先跑 03_module1_brain_eQTL_MR.sh 得 out/smr_brain_mdd2025.smr
## ----------------------------------------------------------------------------
## 输入:
##   out/smr_brain_mdd2025.smr, out/mdd2025.ma, out/brain_tier1_cis_eqtl.txt
##   out/blood_causal_tier1_dual.csv  (来自 02, 脑血比较用)
## 输出:
##   out/brain_causal_gene_tier_table.csv, out/brain_causal_tier1_dual.csv (68)
##   out/coloc_brain.rds, out/brain_blood_comparison.csv
##   out/unified_causal_tier1.csv     (血80+脑68去重=142; 供模块五)
##   out/fig_brain_smr_volcano.png / fig_brain_forest_top20.png /
##     fig_brain_blood_comparison.png (见 python 段)
## 对应正文: 模块一 脑因果基因 + 脑血异同 + 统一 142 因果集
## ----------------------------------------------------------------------------
## coloc 与血同法 (coloc.abf, type cc/quant, s=0.218)。BrainMeta N 见 .sh 注。
## ============================================================================
# 运行前将工作目录设为 workspace 根 (含 out/、data/ 的 ASCII 路径), 例如:
# setwd("/path/to/workspace")   # 切勿用含中文的路径 (R read/system 会 "cannot open the connection")
suppressMessages({library(data.table); library(coloc)})

## ---------------- 1. 注释 + tier (同血脚本逻辑) ----------------
smr <- fread("out/smr_brain_mdd2025.smr")
map <- fread("ensg_symbol_map.tsv", header = TRUE)
smr[, Gene := map$GeneSymbol[match(sub("\\..*", "", probeID), map$ENSG)]]  # 去 ENSG 版本号
smr[, FDR := p.adjust(p_SMR, "BH")]
smr[, heidi_pass := (!is.na(p_HEIDI) & p_HEIDI > 0.01)]
smr[, tier := fifelse(FDR < 0.05 & heidi_pass, "tier1",
              fifelse(FDR < 0.05 & !is.na(p_HEIDI) & !heidi_pass, "tier3_HEIDI_fail",
              fifelse(FDR < 0.05 & is.na(p_HEIDI), "tier4_HEIDI_NA", "ns")))]

## ---------------- 2. coloc.abf (tier1 探针; 同血, s=0.218) ----------------
gw <- fread("out/mdd2025.ma"); setnames(gw, c("SNP","A1_g","A2_g","freq_g","b_g","se_g","p_g","n_g"))
gw <- gw[!is.na(b_g) & se_g > 0]; setkey(gw, SNP)
eq <- fread("out/brain_tier1_cis_eqtl.txt")
s_frac <- 357636 / (357636 + 1281936)
probes <- unique(eq$Probe); res <- vector("list", length(probes))
for (i in seq_along(probes)) {
  pr <- probes[i]; e <- eq[Probe == pr][!is.na(b) & SE > 0]; m <- merge(e, gw, by = "SNP")
  if (nrow(m) < 10) { res[[i]] <- data.table(probe = pr, nsnp = nrow(m), PP.H4 = NA_real_); next }
  D_e <- list(beta = m$b,   varbeta = m$SE^2,   snp = m$SNP, type = "quant", N = 2865,           MAF = pmin(m$Freq, 1 - m$Freq))
  D_g <- list(beta = m$b_g, varbeta = m$se_g^2, snp = m$SNP, type = "cc", s = s_frac, N = round(m$n_g), MAF = pmin(m$freq_g, 1 - m$freq_g))
  ct <- tryCatch(coloc.abf(D_e, D_g), error = function(e) NULL)
  res[[i]] <- data.table(probe = pr, nsnp = nrow(m), PP.H4 = if (is.null(ct)) NA_real_ else ct$summary["PP.H4.abf"])
}
colocdt <- rbindlist(res); saveRDS(colocdt, "out/coloc_brain.rds")
smr <- merge(smr, colocdt[, .(probeID = probe, coloc_PPH4 = PP.H4)], by = "probeID", all.x = TRUE)
smr[, dual_pass := (tier == "tier1" & !is.na(coloc_PPH4) & coloc_PPH4 > 0.8)]
smr[, direction := fifelse(b_SMR > 0, "risk_up", "protective_down")]
fwrite(smr, "out/brain_causal_gene_tier_table.csv")
fwrite(smr[dual_pass == TRUE], "out/brain_causal_tier1_dual.csv")
## 结果: 脑 tier1_dual=68 (risk↑43/protective↓25)。

## ---------------- 3. 脑血比较 + 统一 142 因果集 ----------------
blood <- fread("out/blood_causal_tier1_dual.csv"); brain <- smr[dual_pass == TRUE]
bl <- unique(blood$Gene); br <- unique(brain$Gene)
shared <- intersect(bl, br)                     # 6 共有
uni <- data.table(gene = union(bl, br))
uni[, tissue := fifelse(gene %in% shared, "both", fifelse(gene %in% bl, "blood", "brain"))]
## 方向/coloc 从各自表填入 (略); 已知 MDD 基因保守标注 (DCC/KLF7/... 拿得准才标)。
fwrite(uni, "out/unified_causal_tier1.csv")
## 结果: 142 基因 (血特异74 + 脑特异62 + 共有6); 共有 6 方向一致 3/6。
## 脑/血 tier-1 ∩ 手针89-DEG = 空; ∩ causal_34 = 空(脑)/CNNM2(血) → 因果与针刺两互补层。
