## ============================================================================
## 01_reversal_稳版反转.R
## 手针针刺 GSE86392 反转分析：同研究内反转 + 跨物种反转
## ----------------------------------------------------------------------------
## 输入:
##   - 1676808.f1.pdf                论文补充表原件 (Tab.S1a/S2a=疾病轴 C-vs-M;
##                                    Tab.S1c/S2c=针刺轴 A-vs-M)
##   - prescan_output/MDD_GSE98793_blood_human.csv   人类 MDD 血签名 (limma)
##   - prescan_output/MDD_GSE102556_*_human.csv      人类 MDD 脑 6 区签名
## 输出:
##   - out/within_study_reversal_intersection.csv
##   - out/cross_species_reversal.csv
##   - out/inflammatory_cluster_directions.csv
##   - out/reversal_summary_stats.csv
##   - out/fig1_within_study_scatter.png / fig2_reversal_heatmap.png /
##     fig3_reversal_summary.png
## 对应正文: 图1(同研究内反转散点)、图2(反转热图)、表(反转率+二项检验)
## ----------------------------------------------------------------------------
## 方向约定: logFC>0 = 该处理使基因升高;
##   针刺轴 A-vs-M 直接用 (PDF log2FC=log2(A/M));
##   疾病轴 C-vs-M 翻转符号 (PDF log2FC=log2(C/M) → 疾病=模型升高, 故取负);
##   "反转" = sign(logFC_针刺) != sign(logFC_疾病)。
## 注: PDF 解析与within-study反转在 Python 跑; biomaRt 直系同源映射在 R(geo env)跑。
##     本文件把两段合并记录; 实际运行时 Python 段见文末 reticulate/示意标注。
## 多重检验: 同研究内为单一主检验 FDR=p; 跨物种 7 检验族(6脑+血) BH 校正。
## ============================================================================

## ---------- [Python 段, 示意 — 实际在 python kernel 运行] ----------
## 见 01b_reversal_pdf_parse.py (下方以注释保留 Python 源, 交互探索所得)
## 关键结果: counts S1a=20 S2a=114 S1c=17 S2c=72 | A-vs-M=89DEG;
##           同研究内反转 26/26 (100%), binom p=1.490116e-08 (单一主检验)。

## ---------- R 段: 大鼠→人 直系同源映射 (biomaRt, 主; Ensembl REST 交叉验证) ----------
suppressMessages(library(biomaRt))
## acu_axis_rat.csv: 89 个大鼠针刺轴基因 (ensembl, symbol, logFC)
acu <- read.csv("handoff/acu_axis_rat.csv", stringsAsFactors = FALSE)
rat <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl")
res <- getBM(
  attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type"),
  filters = "ensembl_gene_id", values = acu$ensembl, mart = rat)
res <- res[res$hsapiens_homolog_ensembl_gene != "", ]
write.csv(res, "handoff/biomart_orthologs.csv", row.names = FALSE)
## 结果: 64/89 大鼠基因有≥1 人类直系同源 (65 对); 保留率 71.9%。
## 交叉验证: Ensembl REST /homology/id/rattus_norvegicus/ 给 63 基因/64 对,
##           两法一致 64 对, biomaRt 独有补回 CD14 → 采用 biomaRt 为主。

## ---------- 跨物种反转 (Python 段, 示意) ----------
## 把 65 个人类针刺轴基因分别与 GSE98793 血签名 / GSE102556 脑 6 区签名取交集,
## 算 sign(logFC_acu)!=sign(logFC_disease) 的反转率 + binom.test, BH 校正(7检验族)。
## 结果: 血 GSE98793 反转 35/54=64.8%, binom p=0.0201, BH-FDR=0.1408 (黄灯);
##       脑 6 区 50-69%, 全部 BH>0.05 不显著。

## ============================================================================
## 附: Python 源 (01b_reversal_pdf_parse.py) — PDF 解析 + 同研究内反转
## ============================================================================
## import pypdfium2 as pdf, re, pandas as pd, numpy as np
## from scipy import stats
## Run from the repository root; all paths below are repository relative.
## doc = pdf.PdfDocument(f'{PROJ}/1676808.f1.pdf')
## full = ''.join(doc[i].get_textpage().get_text_range() for i in range(8)).replace('\r','')
## tags = ['Tab.S1a','Tab.S1b','Tab.S1c','Tab.S2a','Tab.S2b','Tab.S2c']
## pos = {t: full.find(t) for t in tags}
## ordered = sorted(pos.items(), key=lambda x:x[1])
## slices = {t: full[p:(ordered[i+1][1] if i+1<len(ordered) else len(full))]
##           for i,(t,p) in enumerate(ordered)}
## row_re3 = re.compile(r'(ENSRNOG\d+)\s+(\S+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+E[+-]\d{2})\s+(-?[\d.]+E[+-]\d{2})')
## def parse3(tag):
##     out=[]
##     for m in row_re3.finditer(slices[tag].replace('\n',' ')):
##         ens,sym,r1,r2,pv,lfc = m.groups()
##         out.append((ens,sym,float(r1),float(r2),float(pv),float(lfc)))
##     return out
## cols = ['ensembl','symbol','rpkm_num','rpkm_M','pval','log2fc_raw']
## s1a,s2a,s1c,s2c = [pd.DataFrame(parse3(t),columns=cols) for t in ['Tab.S1a','Tab.S2a','Tab.S1c','Tab.S2c']]
## acu = pd.concat([s1c,s2c],ignore_index=True); acu['logFC']=acu['log2fc_raw']       # A-vs-M as-is
## dis = pd.concat([s1a,s2a],ignore_index=True); dis['logFC']=-dis['log2fc_raw']      # C-vs-M sign-flip
## def dedup(df):
##     df=df.assign(a=df.logFC.abs()).sort_values('a',ascending=False)
##     return df.drop_duplicates('symbol').drop(columns='a')
## acu_d,dis_d = dedup(acu),dedup(dis)
## merge = pd.merge(dis_d[['symbol','logFC']].rename(columns={'logFC':'logFC_disease'}),
##                  acu_d[['symbol','logFC']].rename(columns={'logFC':'logFC_acu'}), on='symbol')
## merge['reversed'] = np.sign(merge.logFC_disease)!=np.sign(merge.logFC_acu)
## n_int=len(merge); n_rev=int(merge.reversed.sum())
## binom = stats.binomtest(n_rev, n_int, 0.5, alternative='greater')  # 26/26, p=1.490116e-08
