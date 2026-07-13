## ================================================================
## 预扫脚本骨架：针刺信号 vs 抑郁"因果基因/疾病签名"的交集与方向反转
## 目的：判定 module 4（针刺差异化主线）是否成立的命门测试
## 用法：在 Claude Science 中逐步补全 TODO 并运行；每步产物落盘、留存代码与环境
## 环境：R >= 4.3 + Bioconductor
## 重要：本文件是"骨架"，所有数据集编号、分组列、方向定义都必须先核实（见协议文档"AI 验证清单"）
## ================================================================

## ---------- Step 0. 依赖 ----------
## install.packages("BiocManager")
## BiocManager::install(c("GEOquery","limma","GeneOverlap","fgsea",
##                        "clusterProfiler","org.Hs.eg.db","babelgene"))
## remotes::install_github("RRHO2/RRHO2")   # 方向性签名比较
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(GeneOverlap)
  library(fgsea); library(clusterProfiler); library(org.Hs.eg.db)
  library(babelgene)
  ## library(RRHO2)
})
set.seed(1)
OUT <- "prescan_out"; dir.create(OUT, showWarnings = FALSE)

## 约定：所有签名统一为 data.frame(gene, logFC, adjP)；gene = 人类 HGNC symbol
## 方向约定：
##   针刺签名 logFC_acu  > 0 → 针刺使该基因升高（EA vs 模型组）
##   疾病签名 logFC_dis  > 0 → 抑郁使该基因升高（MDD vs 对照）
##   "反转" = sign(logFC_acu) != sign(logFC_dis)  （针刺把疾病方向扳回）

## ---------- helper：大鼠→人 直系同源映射 ----------
map_rat_to_human <- function(genes_rat) {
  ## babelgene::orthologs 返回 rat->human 映射；TODO 核实保留率（应报告丢失多少基因）
  m <- babelgene::orthologs(genes = genes_rat, species = "rat", human = TRUE)
  ## 期望列：symbol(人), rat_symbol；按需去重（多对一取一）
  message(sprintf("[ortholog] 输入 %d 个大鼠基因，映射到 %d 个人类 symbol",
                  length(unique(genes_rat)), length(unique(m$symbol))))
  m
}

## ================================================================
## Step 1. 构建 AcuDEG（针刺响应差异基因，带方向）
## ================================================================
## 候选数据集（抑郁相关优先；须逐一核实 accession/分组/平台）：
##   [首选] 大鼠 CUMS + 电针(百会/印堂) 海马  —— 抑郁范式，穴位对口
##   [次选] 大鼠 CUMS + 电针 下丘脑 (Neurochem Res 2024, PMID 38522048)
##   [次选] 大鼠 CUMS + 电针 (Brain & Behavior 2024)
##   [人但非抑郁] ST36 电针健康人全血 (BMC Complement Med 2017) —— 仅证"针刺可调控"，不用于反转
## 数据来源两条路：
##   A) GEO 原始数据：getGEO() + limma（下方骨架）
##   B) 若未存 GEO：从论文附表提取 DEG（gene + logFC + adjP），直接读入
##   —— 现实中不少针刺转录组只在附表给 DEG，B 路常用

build_acu_deg_from_GEO <- function(gse_id, group_col, level_EA, level_model) {
  gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
  ex <- exprs(gset); ph <- pData(gset)
  ## TODO 核实：分组列名 group_col、EA 组/模型组的取值；必要时 log2 变换与归一化
  grp <- factor(ph[[group_col]], levels = c(level_model, level_EA))
  design <- model.matrix(~ grp)                     # 第2列 = EA vs 模型
  fit <- eBayes(lmFit(ex, design))
  tt <- topTable(fit, coef = 2, number = Inf)
  ## TODO：把探针注释到 rat symbol（用 fData(gset) 或平台注释），再映射到人
  tt$rat_symbol <- NA                               # TODO fill
  map <- map_rat_to_human(tt$rat_symbol)
  ## 合并 → 人类 symbol 级 logFC（多探针取 |logFC| 最大或均值，TODO 决定）
  acu <- data.frame(gene = map$symbol,
                    logFC = tt$logFC[match(map$rat_symbol, tt$rat_symbol)],
                    adjP  = tt$adj.P.Val[match(map$rat_symbol, tt$rat_symbol)])
  acu <- aggregate(cbind(logFC) ~ gene, data = acu, FUN = function(x) x[which.max(abs(x))])
  acu
}

## 或：从附表读入（人类或大鼠 symbol）
## acu_deg <- read.csv("acu_supp_DEG.csv")  # 列: gene/logFC/adjP；若大鼠则先 map_rat_to_human

acu_deg <- NULL      # TODO：赋值为最终 data.frame(gene, logFC, adjP)
## 显著子集（阈值可调；针刺数据样本小，可放宽到 P<0.05 未校正并在协议中说明）
## acu_sig <- subset(acu_deg, adjP < 0.05)

## ================================================================
## Step 2. 构建疾病签名 disease_deg（MDD vs 对照，带方向）
## ================================================================
## 候选：
##   [血] GSE98793（MDD 全血，常用）；GSE76826 / GSE38206 作复现
##   [脑] GSE102556（MDD 多脑区尸检，Labonté 2017）—— 脑更贴近致病组织
## 血/脑分别做，后续比较脑血一致性

build_disease_deg_GEO <- function(gse_id, group_col, level_case, level_ctrl) {
  gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
  ex <- exprs(gset); ph <- pData(gset)
  ## TODO 核实分组列与取值；批次/协变量（年龄性别）如可得则纳入 design
  grp <- factor(ph[[group_col]], levels = c(level_ctrl, level_case))
  design <- model.matrix(~ grp)
  fit <- eBayes(lmFit(ex, design))
  tt <- topTable(fit, coef = 2, number = Inf)
  ## TODO：探针 → 人类 symbol（fData/平台注释），多探针聚合
  tt$gene <- NA                                     # TODO fill (human symbol)
  dis <- data.frame(gene = tt$gene, logFC = tt$logFC, adjP = tt$adj.P.Val)
  dis <- subset(dis, !is.na(gene) & gene != "")
  aggregate(cbind(logFC) ~ gene, data = dis, FUN = function(x) x[which.max(abs(x))])
}

disease_deg <- NULL  # TODO：赋值（血液版）
## disease_deg_brain <- build_disease_deg_GEO("GSE102556", ...)  # 可选脑版

## ================================================================
## Step 3. 抑郁"因果基因"集 causal_genes（带方向 if 可得）
## ================================================================
## 路线 A（快，零算力）：从已发表抑郁 MR/SMR/TWAS 论文编制基因清单
##   参考：Transl Psychiatry 2021 (eQTL-MR)、Molecular Psychiatry 2022 (脑+血 蛋白+转录组)、
##         J Affect Disord 2025 (脑单细胞 eQTL-MR)、Neuropsychopharmacology 2022 (全基因组 MR 药靶)
##   产物：causal_genes = character 向量（若能提取效应方向则另存 sign）
## 路线 B（可选，需算力/数据）：自跑 SMR
##   前置：MDD GWAS 用 rmUKBB 版（避免与 eQTLGen/UKB-PPP 样本重叠）；LD 参考(1000G EUR)
##   命令示例（shell，非 R）：
##     smr --bfile g1000_eur --gwas-summary MDD_rmUKBB.ma \
##         --beqtl-summary eQTLGen --peqtl-smr 5e-8 --out smr_blood
##     smr --bfile g1000_eur --gwas-summary MDD_rmUKBB.ma \
##         --beqtl-summary MetaBrain_cortex --out smr_brain
##   读取 *.smr：保留 p_SMR < 0.05/FDR 且 p_HEIDI > 0.05
read_smr <- function(path) {
  x <- read.table(path, header = TRUE, sep = "\t")
  subset(x, p_HEIDI > 0.05 & p.adjust(p_SMR, "fdr") < 0.05)$Gene
}
causal_genes <- NULL # TODO：character 向量（人类 symbol）

## ================================================================
## Step 4. 交集显著性（针刺 ∩ 因果、针刺 ∩ 疾病）
## ================================================================
## 背景基因数：取三方共同可测基因的并/交，作 genome.size（TODO 决定，建议用共同背景）
BG <- 20000  # TODO：改为实际共同背景基因数

overlap_test <- function(a, b, bg = BG) {
  go <- newGeneOverlap(unique(a), unique(b), genome.size = bg)
  go <- testGeneOverlap(go)
  list(n_overlap = length(intersect(a, b)),
       jaccard   = getJaccard(go),
       p         = getPval(go),
       genes     = intersect(a, b))
}

## 需要显著基因子集（阈值见协议）
acu_sig     <- unique(subset(acu_deg,     adjP < 0.05)$gene)   # TODO 阈值
dis_sig     <- unique(subset(disease_deg, adjP < 0.05)$gene)   # TODO 阈值
res_acu_causal  <- overlap_test(acu_sig, causal_genes)
res_acu_disease <- overlap_test(acu_sig, dis_sig)
print(res_acu_causal[c("n_overlap","jaccard","p")])
print(res_acu_disease[c("n_overlap","jaccard","p")])

## ================================================================
## Step 5. 方向反转一致性（核心：针刺是否把疾病方向扳回）
## ================================================================
## 5a. 交集基因的符号反转检验
m <- merge(acu_deg[,c("gene","logFC")], disease_deg[,c("gene","logFC")],
           by = "gene", suffixes = c("_acu","_dis"))
m$reversal <- sign(m$logFC_acu) != sign(m$logFC_dis)
cat(sprintf("交集基因 %d 个，其中方向反转 %d 个\n", nrow(m), sum(m$reversal)))
bt <- binom.test(sum(m$reversal), nrow(m), p = 0.5, alternative = "greater")
print(bt)     # p 小 → 反转显著多于随机

## 5b.（推荐）全签名 RRHO2：不依赖阈值，看两个排序签名的方向一致/反转象限
## 排序度量：signed score = sign(logFC) * -log10(adjP)
## la <- data.frame(gene=acu_deg$gene,     score = sign(acu_deg$logFC)*-log10(acu_deg$adjP))
## ld <- data.frame(gene=disease_deg$gene, score = sign(disease_deg$logFC)*-log10(disease_deg$adjP))
## common <- intersect(la$gene, ld$gene)
## obj <- RRHO2::RRHO2_initialize(la[la$gene %in% common,], ld[ld$gene %in% common,],
##                                labels=c("Acupuncture","MDD"), log10.ind=TRUE)
## RRHO2::RRHO2_heatmap(obj)   # 强"discordant"象限 = 反转证据

## ================================================================
## Step 6. 通路层面一致性（比基因级更稳，跨物种更可比）
## ================================================================
## 用疾病签名做排序，fgsea 检验"针刺升高基因集""针刺降低基因集"的富集方向
rank_dis <- with(disease_deg, setNames(sign(logFC)*-log10(adjP), gene))
rank_dis <- sort(rank_dis[is.finite(rank_dis)], decreasing = TRUE)
acu_sets <- list(
  Acu_UP = unique(subset(acu_deg, logFC > 0 & adjP < 0.05)$gene),
  Acu_DN = unique(subset(acu_deg, logFC < 0 & adjP < 0.05)$gene)
)
fg <- fgsea(pathways = acu_sets, stats = rank_dis, minSize = 5)
print(fg[, c("pathway","NES","padj")])
## 期望"反转"信号：Acu_UP 在疾病 DOWN 端富集(NES<0)、Acu_DN 在疾病 UP 端富集(NES>0)

## （可选）通路富集交叉：分别对 acu_sig / dis_sig / causal_genes 跑 KEGG/GO，看共享通路
## ek_acu <- enrichKEGG(bitr(acu_sig, "SYMBOL","ENTREZID", org.Hs.eg.db)$ENTREZID)
## ek_dis <- enrichKEGG(...); 比较 top 通路重叠

## ================================================================
## Step 7. 汇总为一张决策表
## ================================================================
decision <- data.frame(
  metric = c("针刺∩因果 富集p", "针刺∩疾病 富集p",
             "交集方向反转 binom p", "通路反转(fgsea)是否成立"),
  value  = c(signif(res_acu_causal$p,3), signif(res_acu_disease$p,3),
             signif(bt$p.value,3), "TODO：看 fg 方向")
)
write.csv(decision, file.path(OUT, "prescan_decision.csv"), row.names = FALSE)
print(decision)

## —— 决策标准见协议文档"go/no-go"——
## 绿灯：交集富集显著(p<0.05) + 反转显著(binom p<0.05) + 通路反转成立 → 按修改稿铺开
## 黄灯：部分成立 → 限定更贴合的针刺数据/穴位后重试，或降低模块四权重
## 红灯：均不显著 → 放弃"针刺反转"主线，改走纯遗传学(脑+血多组学 MR 找可成药靶)路线
