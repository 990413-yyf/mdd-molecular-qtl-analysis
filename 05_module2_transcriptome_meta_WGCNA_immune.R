## ============================================================================
## 05_module2_transcriptome_meta_WGCNA_immune.R
## 模块二: 四队列外周血转录组 meta + WGCNA + 免疫浸润
## ----------------------------------------------------------------------------
## 输入 (GEOquery 直下, 存 data/geo/<GSE>_raw.rds):
##   GSE98793 (GPL570, 128 MDD/64 CTRL, 全血)
##   GSE76826 (GPL17077 Agilent, 20/12)
##   GSE38206 (GPL13607 Agilent, 18 MDE/18 CTRL, PBMC)
##   GSE39653 (GPL10558 Illumina, 21 MDD/24 CTRL, 剔除 8 双相)
##   out/unified_causal_tier1.csv  (来自模块一, 142 因果基因, 用于两层对照)
## 输出:
##   out/module2_mdd_blood_meta_DEG.csv                (14857 基因 meta 签名)
##   out/module2_inflammatory_cluster_in_disease_sig.csv (针刺反转簇17基因方向)
##   out/module2_immune_infiltration.csv               (7 细胞类型 MDD vs CTRL)
##   out/module2_wgcna_module_membership.csv
##   out/fig_module2_volcano.png / fig_module2_wgcna.png / fig_module2_immune_GO.png
## 对应正文: 模块二 (下游炎症 vs 上游因果对照, 疾病签名, WGCNA, 免疫浸润)
## ----------------------------------------------------------------------------
## 铁律: 每队列先打印分组核实; 多重检验 BH; 换方法自我复核(RE/FE/Stouffer 三法)。
## R 包: GEOquery/limma/data.table/metafor (conda), WGCNA (overlay CRAN source)。
## 免疫浸润在 python kernel 用标志基因打分 (避 CIBERSORT-LM22 组织错配); 见文末。
## ============================================================================
.libPaths(c("./.r-libs/geo", .libPaths()))
suppressMessages({library(GEOquery); library(limma); library(data.table); library(metafor)})
options(timeout = 900)

## ---------------- Step 1. 下载 + 核实分组 (每队列打印病例/对照) ----------------
gses <- c("GSE98793", "GSE76826", "GSE38206", "GSE39653")
for (g in gses) {
  gse <- getGEO(g, destdir = "data/geo", GSEMatrix = TRUE, getGPL = TRUE)
  e <- gse[[1]]
  saveRDS(list(ex = exprs(e), pd = pData(e), fd = fData(e), annot = annotation(e)),
          paste0("data/geo/", g, "_raw.rds"))
}
## 分组定义 (先打印 pData 特征列核实后再用)
get_groups <- function(g, pd) {
  if (g == "GSE98793") {
    ifelse(grepl("CASE", pd$characteristics_ch1), "MDD",
    ifelse(grepl("CNTL", pd$characteristics_ch1), "CTRL", NA))
  } else if (g == "GSE76826") {
    dx <- pd[["diagnosis:ch1"]]; if (is.null(dx)) dx <- pd$characteristics_ch1.1
    ifelse(grepl("depress", dx, ignore.case = TRUE), "MDD",
    ifelse(grepl("Healthy", dx, ignore.case = TRUE), "CTRL", NA))
  } else if (g == "GSE38206") {
    ifelse(grepl("MDE|major depress", pd$characteristics_ch1, ignore.case = TRUE), "MDD",
    ifelse(grepl("control", pd$characteristics_ch1, ignore.case = TRUE), "CTRL", NA))
  } else if (g == "GSE39653") {
    dx <- pd[["disease:ch1"]]
    ifelse(grepl("major depress", dx, ignore.case = TRUE), "MDD",
    ifelse(grepl("healthy", dx, ignore.case = TRUE), "CTRL", NA))   # 剔除双相→NA
  }
}
## 核实结果: 98793=128/64, 76826=20/12, 38206=18/18, 39653=21/24 (8双相 NA)

## ---------------- Step 2. 各队列 limma DE (probe→symbol max-mean 折叠) ----------------
sym_col <- list(GSE98793 = "Gene Symbol", GSE76826 = "GENE_SYMBOL",
                GSE38206 = "GeneName", GSE39653 = "ILMN_Gene")
de_one <- function(g) {
  d <- readRDS(paste0("data/geo/", g, "_raw.rds"))
  ex <- d$ex; pd <- d$pd; fd <- d$fd
  grp <- get_groups(g, pd); keep <- !is.na(grp); ex <- ex[, keep]; grp <- grp[keep]
  if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)          # GSE38206 需 log2
  sym <- as.character(fd[[sym_col[[g]]]]); names(sym) <- rownames(fd); sym <- sym[rownames(ex)]
  sym <- trimws(sub("[ ]*///.*", "", sym))
  ok <- !is.na(sym) & sym != "" & sym != "---"; ex <- ex[ok, ]; sym <- sym[ok]
  o <- order(rowMeans(ex, na.rm = TRUE), decreasing = TRUE)   # max-mean 探针
  ex <- ex[o, ]; sym <- sym[o]; dup <- duplicated(sym); ex <- ex[!dup, ]; rownames(ex) <- sym[!dup]
  design <- model.matrix(~factor(grp, levels = c("CTRL", "MDD")))
  fit <- eBayes(lmFit(ex, design)); tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(ex), logFC = tt$logFC, t = tt$t, P = tt$P.Value,
             adjP = tt$adj.P.Val, se = tt$logFC / tt$t,
             n_mdd = sum(grp == "MDD"), n_ctrl = sum(grp == "CTRL"), stringsAsFactors = FALSE)
}
res <- lapply(gses, de_one); names(res) <- gses
saveRDS(res, "data/geo/de_per_cohort.rds")

## ---------------- Step 2b. Meta 分析 (三法: RE 主保守 / FE / Stouffer 主报) ----------------
common <- Reduce(intersect, lapply(res, function(x) x$gene))   # 14857 共有基因
FC <- sapply(res, function(x) x[match(common, x$gene), "logFC"])
SE <- sapply(res, function(x) x[match(common, x$gene), "se"])
rownames(FC) <- common; rownames(SE) <- common
SE[!is.finite(SE) | SE <= 0] <- NA
## 随机效应 REML (metafor) — 最保守, 惩罚跨平台异质性 → 0 显著
meta_one <- function(i) {
  yi <- FC[i, ]; sei <- SE[i, ]; ok <- is.finite(yi) & is.finite(sei)
  if (sum(ok) < 3) return(c(NA, NA, NA, NA))
  mm <- tryCatch(rma(yi = yi[ok], sei = sei[ok], method = "REML"),
                 error = function(e) tryCatch(rma(yi = yi[ok], sei = sei[ok], method = "DL"),
                                              error = function(e2) NULL))
  if (is.null(mm)) return(c(NA, NA, NA, NA))
  c(mm$beta, mm$se, mm$pval, mm$I2)
}
M <- t(sapply(seq_along(common), meta_one))
colnames(M) <- c("meta_logFC", "meta_se", "meta_P", "I2")
meta <- data.frame(gene = common, M); meta <- meta[is.finite(meta$meta_P), ]
meta$meta_adjP <- p.adjust(meta$meta_P, "BH")               # RE: 0 显著

## 固定效应 IVW (第二法)
meta_fe <- function(i) {
  yi <- FC[i, ]; sei <- SE[i, ]; ok <- is.finite(yi) & is.finite(sei)
  if (sum(ok) < 3) return(c(NA, NA, NA))
  w <- 1 / sei[ok]^2; b <- sum(w * yi[ok]) / sum(w); se <- sqrt(1 / sum(w))
  c(b, se, 2 * pnorm(-abs(b / se)))
}
FE <- t(sapply(seq_along(common), meta_fe)); colnames(FE) <- c("fe_logFC", "fe_se", "fe_P")
fe <- data.frame(gene = common, FE); fe <- fe[is.finite(fe$fe_P), ]; fe$fe_adjP <- p.adjust(fe$fe_P, "BH")

## Stouffer N-加权 Z (主报; 跨平台 GEO meta 标准)
Zmat <- sapply(res, function(x) { r <- x[match(common, x$gene), ]
  sign(r$logFC) * qnorm(pmax(r$P / 2, 1e-300), lower.tail = FALSE) })
Ns <- sapply(res, function(x) x$n_mdd[1] + x$n_ctrl[1]); rownames(Zmat) <- common
stouffer <- function(i) { z <- Zmat[i, ]; ok <- is.finite(z); if (sum(ok) < 3) return(c(NA, NA))
  w <- sqrt(Ns[ok]); Z <- sum(w * z[ok]) / sqrt(sum(w^2)); c(Z, 2 * pnorm(-abs(Z))) }
ST <- t(sapply(seq_along(common), stouffer)); colnames(ST) <- c("stouffer_Z", "stouffer_P")
st <- data.frame(gene = common, ST); st <- st[is.finite(st$stouffer_P), ]; st$stouffer_adjP <- p.adjust(st$stouffer_P, "BH")
## 结果: RE=0, FE=380, Stouffer=121 显著 (BH<0.05)。

## 统一签名: primary_logFC=固定效应, primary_adjP=Stouffer
mm <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE),
   list(meta[, c("gene", "meta_logFC", "meta_P", "meta_adjP", "I2")],
        fe[, c("gene", "fe_logFC", "fe_P", "fe_adjP")],
        st[, c("gene", "stouffer_Z", "stouffer_P", "stouffer_adjP")]))
mm$primary_logFC <- mm$fe_logFC; mm$primary_adjP <- mm$stouffer_adjP
mm <- mm[order(mm$stouffer_P), ]
mm$n_cohorts <- sapply(mm$gene, function(g) sum(sapply(res, function(x) g %in% x$gene)))
write.csv(mm, "out/module2_mdd_blood_meta_DEG.csv", row.names = FALSE)

## ---------------- Step 3. 针刺反转炎症簇在疾病签名中的方向 ----------------
infl <- c("CXCL9","CXCL10","CCL2","C3","A2M","SERPINA3","PTX3","LGALS3","S100A9",
          "S100A4","UBD","CXCL11","CCL4","IL15RA","IL1A","CD14","IRG1","ACOD1","TTR")
ic <- mm[mm$gene %in% infl, ]
ic$direction_in_MDD <- ifelse(ic$fe_logFC > 0, "UP", "DOWN")
ic <- ic[order(-ic$stouffer_Z), ]
write.csv(ic, "out/module2_inflammatory_cluster_in_disease_sig.csv", row.names = FALSE)
## 结果: 12/17 上调 (binom p=0.072); PTX3 显著上调 (Stouffer adjP=0.005)。
## 对照: 142 因果 tier-1 基因在疾病签名 0/76 达 Stouffer BH<0.05 (下游 vs 上游)。

## ---------------- Step 4. WGCNA (GSE98793, 最大队列 192 样本) ----------------
suppressMessages(library(WGCNA)); enableWGCNAThreads(4)
d <- readRDS("data/geo/GSE98793_raw.rds"); ex <- d$ex; pd <- d$pd; fd <- d$fd
grp <- ifelse(grepl("CASE", pd$characteristics_ch1), 1, ifelse(grepl("CNTL", pd$characteristics_ch1), 0, NA))
anx <- ifelse(grepl("anxiety: yes", pd$characteristics_ch1.1), 1, 0)
sym <- trimws(sub("[ ]*///.*", "", as.character(fd[["Gene Symbol"]]))); names(sym) <- rownames(fd); sym <- sym[rownames(ex)]
ok <- !is.na(sym) & sym != "" & sym != "---"; ex <- ex[ok, ]; sym <- sym[ok]
o <- order(rowMeans(ex), decreasing = TRUE); ex <- ex[o, ]; sym <- sym[o]
dup <- duplicated(sym); ex <- ex[!dup, ]; rownames(ex) <- sym[!dup]
top <- names(sort(apply(ex, 1, var), decreasing = TRUE))[1:8000]  # top-8000 变异基因
datExpr <- t(ex[top, ])
sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)), networkType = "signed", verbose = 0)
## soft-power=7 (首个 signed R2>0.8)
net <- blockwiseModules(datExpr, power = 7, networkType = "signed", TOMType = "signed",
        minModuleSize = 40, mergeCutHeight = 0.20, numericLabels = TRUE,
        pamRespectsDendro = FALSE, maxBlockSize = 8000, verbose = 0, saveTOMs = FALSE)
mc <- labels2colors(net$colors)
MEs <- orderMEs(moduleEigengenes(datExpr, mc)$eigengenes)
mt <- cor(MEs, data.frame(MDD = grp, anx = anx), use = "p")
mp <- corPvalueStudent(mt, nrow(datExpr))
write.csv(data.frame(gene = colnames(datExpr), module = mc),
          "out/module2_wgcna_module_membership.csv", row.names = FALSE)
## 结果: 6 模块; 唯一显著关联 MDD = MEgrey (未分配基因, r=-0.26, p=3e-4);
##       无相干共表达模块追踪 MDD → 信号弥散分布式。模块-性状热图见 fig_module2_wgcna.png
##       (labeledHeatmap, blueWhiteRed)。

## ---------------- Step 6. 免疫浸润 (Python 段, 示意 — 在 python kernel 运行) ----------------
## 用血液适配的细胞类型标志基因集打分 (per-sample marker z-score 均值), 避 LM22 组织错配。
## 见附 05b_immune_deconvolution.py:
##   markers = {Neutrophils, Monocytes, T_cells, CD8_T, B_cells, NK_cells, Dendritic}
##   z = (logx - mean)/sd per gene; score[ct] = mean z of markers; MannWhitney MDD vs CTRL; BH。
## 结果: 单核↑(FDR=0.034) 中性粒↑; CD8 T↓(0.006) NK↓ DC↓ (0.034) → 髓系偏移炎症浸润景观。
