#!/usr/bin/env python3
## ============================================================================
## 05b_immune_deconvolution.py  (模块二 Step 6; 由 05 调用/示意)
## 免疫浸润: 血液适配的细胞类型标志基因集打分 (避 CIBERSORT-LM22 组织错配)
## ----------------------------------------------------------------------------
## 输入:  data/geo/GSE98793_expr_linear.csv, data/geo/GSE98793_groups.csv
##         (由 05 从 GSE98793_raw.rds 导出的 gene-level 线性表达 + 分组)
## 输出:  out/module2_immune_infiltration.csv  (7 细胞类型 MDD vs CTRL, MannWhitney+BH)
## 对应正文: 模块二 免疫浸润 (髓系↑淋巴↓炎症景观)
## 方法: 标志基因 z-score (per gene) → per-sample marker 均值 = 富集分; MDD vs CTRL BH。
## ============================================================================
import pandas as pd, numpy as np
from scipy import stats

WS = "."
expr = pd.read_csv(f"{WS}/data/geo/GSE98793_expr_linear.csv", index_col=0)
grp  = pd.read_csv(f"{WS}/data/geo/GSE98793_groups.csv")

markers = {
 "Neutrophils": ["FCGR3B","CSF3R","FUT4","CEACAM3","S100A8","S100A9","S100A12","MMP9","CXCR1","CXCR2"],
 "Monocytes":   ["CD14","CSF1R","LYZ","FCN1","VCAN","S100A12","CLEC7A","CD68"],
 "T_cells":     ["CD3D","CD3E","CD3G","CD2","CD28","LCK","IL7R","CD5"],
 "CD8_T":       ["CD8A","CD8B","GZMK","GZMA","NKG7","CCL5"],
 "B_cells":     ["CD19","MS4A1","CD79A","CD79B","IGHM","BANK1"],
 "NK_cells":    ["NCAM1","KLRD1","KLRF1","NKG7","GNLY","KLRB1","NCR1"],
 "Dendritic":   ["FCER1A","CLEC10A","CD1C","IRF8","ITGAX"],
}
logx = np.log2(expr + 1)
z = logx.sub(logx.mean(axis=1), axis=0).div(logx.std(axis=1) + 1e-9, axis=0)
scores = {ct: z.loc[[g for g in gs if g in z.index]].mean(axis=0) for ct, gs in markers.items()}
S = pd.DataFrame(scores)
S["group"] = grp.set_index("sample").loc[S.index, "group"].values
S.to_csv(f"{WS}/data/geo/immune_scores.csv")

rows = []
for ct in markers:
    a = S.loc[S.group == "MDD", ct]; b = S.loc[S.group == "CTRL", ct]
    u, p = stats.mannwhitneyu(a, b)
    rows.append((ct, a.mean(), b.mean(), a.mean() - b.mean(), p))
imm = pd.DataFrame(rows, columns=["cell_type", "MDD_mean", "CTRL_mean", "diff", "p"])
imm["adjP_BH"] = stats.false_discovery_control(imm.p)
imm.to_csv(f"{WS}/out/module2_immune_infiltration.csv", index=False)
## 结果: 单核↑(FDR=0.034) 中性粒↑(趋势); CD8 T↓(0.006) NK↓(0.034) DC↓(0.034) T↓(趋势)
##   → 髓系偏移+淋巴减少的炎症浸润景观, 与上调炎症签名互印。
