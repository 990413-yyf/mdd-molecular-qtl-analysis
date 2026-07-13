## ============================================================================
## 06_module2_GO_KEGG_enrichment.R
## 模块二 疾病签名功能富集: GO-BP + KEGG (clusterProfiler, BH-FDR)
## ----------------------------------------------------------------------------
## 输入:
##   - out/module2_mdd_blood_meta_DEG.csv   (来自 05; 14857 基因 meta 签名)
## 输出:
##   - out/module2_GO_up_enrichment.csv     (GO-BP, 上调基因, 411 项 p.adjust<0.1)
##   - out/module2_KEGG_up_enrichment.csv   (KEGG hsa, 316 通路, 42 项 BH<0.05)
##   - out/fig_module2_immune_GO.png        (免疫浸润 + GO, 见 05/本脚本)
##   - out/fig_module2_KEGG.png             (KEGG 条形图; 见附 Python 段)
## 对应正文: 模块二 第五节功能富集 (炎症在表达层失调=下游状态)
## ----------------------------------------------------------------------------
## 输入基因集: 疾病签名中 Stouffer 名义 p<0.01 的上调基因 (623 个), GO/KEGG 一致。
## 多重检验: 一律 BH。
## R 包: clusterProfiler(conda), org.Hs.eg.db + GO.db (overlay ./.r-libs/geo,
##        由 bioconductor tarball 经 R CMD INSTALL 装入)。
## ============================================================================
.libPaths(c("./.r-libs/geo", .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db)})

m <- read.csv("out/module2_mdd_blood_meta_DEG.csv")
## 复现与 GO 完全一致的输入: Stouffer 名义 p<0.01 且固定效应 logFC>0
sig_up <- m$gene[m$stouffer_P < 0.01 & m$fe_logFC > 0]
stopifnot(length(sig_up) == 623)                       # 核实: 623 个上调基因
eg_up <- bitr(sig_up, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID   # → 604 映射
eg_bg <- bitr(m$gene,  "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID  # 背景 14013

## ---------------- GO-BP 富集 (上调基因) ----------------
sig_dn <- m$gene[m$stouffer_P < 0.01 & m$fe_logFC < 0]  # 下调 200 (备用)
go_up <- enrichGO(eg_up, org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                  universe = eg_bg, pvalueCutoff = 0.1, qvalueCutoff = 0.25,
                  readable = TRUE)
write.csv(go_up@result[go_up@result$p.adjust < 0.1, ],
          "out/module2_GO_up_enrichment.csv", row.names = FALSE)
## 结果: 411 项 p.adjust<0.1; top: leukocyte mediated immunity (p.adj=5.7e-6),
##       cytokine-mediated signaling (1.1e-4), chemotaxis, cellular response to LPS。

## ---------------- KEGG 富集 (上调基因) ----------------
## 技术注 (重要, 可重复性): clusterProfiler 的通路分类资源 kegg_category.rda
## 默认从 https://yulab-smu.top/clusterProfiler 下载, 本运行环境该域名 403;
## 其 GitHub gh-pages 镜像可达。解决: 从镜像取得 .rda, 注入包缓存环境后跑通。
## KEGG 通路本体经 rest.kegg.jp 在线获取 (enrichKEGG 默认)。
##
## (a) 取 kegg_category.rda (bash, 已在会话中执行 — 示意命令):
##   curl -L -o ./.cache/R/yulab.utils/kegg_category.rda \
##     https://raw.githubusercontent.com/YuLab-SMU/clusterProfiler/refs/heads/gh-pages/kegg_category.rda
##
## (b) 注入 yulab.utils 缓存环境, 并 patch get_cached_kegg_data 直接读缓存:
load("./.cache/R/yulab.utils/kegg_category.rda")        # 载入 kegg_category (586x4)
nsy <- getNamespace("yulab.utils")
cache_env <- get(".yulabCache", envir = nsy)
assign("kegg_category", kegg_category, envir = cache_env)
nscp <- getNamespace("clusterProfiler")
patched <- function(type = "category") {
  type <- match.arg(type, c("category", "species"))
  basefile <- sprintf("kegg_%s", type)
  env <- yulab.utils:::get_cache()
  if (!exists(basefile, envir = env)) stop("no cached ", basefile)
  get(basefile, envir = env, inherits = FALSE)
}
environment(patched) <- nscp
unlockBinding("get_cached_kegg_data", nscp)
assign("get_cached_kegg_data", patched, envir = nscp)
lockBinding("get_cached_kegg_data", nscp)

## (c) 跑 enrichKEGG (organism='hsa', BH 校正):
options(clusterProfiler.download.method = "auto", timeout = 300)
kk <- enrichKEGG(gene = eg_up, organism = "hsa", keyType = "kegg",
                 pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                 universe = eg_bg, minGSSize = 5, maxGSSize = 500)
kk <- setReadable(kk, org.Hs.eg.db, "ENTREZID")
write.csv(kk@result, "out/module2_KEGG_up_enrichment.csv", row.names = FALSE)
## 结果: 316 通路受检, 42 项 BH<0.05, 64 项 BH<0.10。
## 免疫/炎症命中 (BH<0.05): Phagocytosis(0.0013), MAPK(0.0047), TNF signaling(0.0077),
##   NF-kappa B signaling(0.0081), leukocyte transendothelial migration(0.0099),
##   chemokine signaling(0.012), Toll-like receptor signaling(0.013)。
## → 精确落在 TLR→NF-κB→TNF 轴, 与 GO 一致, 坐实炎症在表达层失调=下游状态。

## ---------------- KEGG 条形图 (附 Python 段, 示意 — 在 python kernel 绘) ----------
## import pandas as pd, numpy as np, matplotlib.pyplot as plt
## r = pd.read_csv("out/module2_KEGG_up_enrichment.csv")
## kw = "cytokine|chemokine|Toll|NF-kappa|TNF|MAPK|phagocy|leukocyte|IL-17|Neutrophil extracellular"
## r["immune"] = r.Description.str.contains(kw, case=False)
## top = r.sort_values("p.adjust").head(18).iloc[::-1]
## ... barh(-log10(p.adjust)), 免疫红/其他灰, 存 out/fig_module2_KEGG.png
