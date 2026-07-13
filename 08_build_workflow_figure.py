#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
图1 研究流程（中/英双版，黑白、直角线条）——仅流程图，不含证据总览表。
- 纯黑白:白底黑框,连接线为直角折线(orthogonal),黑色实心箭头。
- 措辞已按因果语言校准（候选因果基因 / 未见支持性证据 / 手针描述性反转）。
用法:python3 08_build_workflow_figure.py
"""
import os, math

W, H = 1300, 800
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Figures_Main")

FONT = {
  "cn": "'Noto Sans CJK SC','PingFang SC','Microsoft YaHei',sans-serif",
  "en": "'Helvetica Neue','Helvetica','Arial',sans-serif",
}

NODES = {
  "gwas":  dict(x=350, y=88,  w=600, h=58),
  "blood": dict(x=40,  y=196, w=386, h=86),
  "brain": dict(x=457, y=196, w=386, h=86),
  "pqtl":  dict(x=874, y=196, w=386, h=86),
  "mr":    dict(x=40,  y=330, w=1220,h=58),
  "causal":dict(x=40,  y=452, w=590, h=126),
  "null":  dict(x=670, y=452, w=590, h=126),
  "expr":  dict(x=40,  y=648, w=590, h=104),
  "acu":   dict(x=670, y=648, w=590, h=104),
}

TXT = {
 "cn": {
  "title": "抑郁症：遗传因果信号与疾病相关炎症状态——研究流程",
  "subtitle": "多组织多组学孟德尔随机化 · 数字对应正文方法与结果",
  "band": "疾病相关状态层（描述性 / 假设生成）",
  "gwas": ["结局 GWAS：PGC MDD2025（Adams 2025）",
           "357,636 病例 / 1,281,936 对照（去 UKB 欧洲人群）"],
  "blood": ["暴露①  全血 cis-eQTL", "eQTLGen（n=31,684）"],
  "brain": ["暴露②  脑皮层 cis-eQTL", "BrainMeta v2（2,865 样本/2,443 个体）"],
  "pqtl":  ["暴露③  91 循环炎症蛋白 cis-pQTL",
            "SCALLOP（n≈14,736）+ UKB-PPP（n≈34,557）"],
  "mr":    ["孟德尔随机化：SMR + HEIDI + 共定位（PP.H4>0.8）",
            "＋ Steiger 方向性检验 ＋ Brion 效力/最小可检测效应"],
  "causal":["候选基因优先排序",
            "142 个 tier-1 候选基因 → 15 个 tier-A 候选靶点",
            "替代 eQTL 资源稳健性：可测者 86% 达复现标准、方向一致 99%"],
  "null":  ["当前可工具化循环炎症蛋白未获稳健支持",
            "两平台可测 74 种（含 TNF/CCL2/IFNG）：0 通过 MR＋共定位",
            "SCALLOP 可测蛋白在假设模型下对 OR≥1.10 有≥80% 效力；17 种无工具"],
  "expr":  ["表达状态层：四队列外周血 meta（178/109）",
            "炎症相关状态信号（PTX3 BH=0.038）；髓系浸润升高",
            "142 因果基因中 76 个可测者 0/76 显著（效应小、效力有限）"],
  "acu":   ["手针签名（GSE86392）——假设生成、不与主 MR 并列",
            "同研究内方向一致反转 26/26（描述性）",
            "跨物种血 64.8%（FDR=0.14，未达显著）"],
 },
 "en": {
  "title": "Depression: genetic causal signal vs disease-related inflammatory state — study workflow",
  "subtitle": "Multi-omics Mendelian randomization; figures correspond to Methods and Results",
  "band": "Disease-associated state layer (descriptive / hypothesis-generating)",
  "gwas": ["Outcome GWAS: PGC MDD2025 (Adams 2025)",
           "357,636 cases / 1,281,936 controls (non-UKB European)"],
  "blood": ["Exposure 1  whole-blood cis-eQTL", "eQTLGen (n=31,684)"],
  "brain": ["Exposure 2  brain cortex cis-eQTL", "BrainMeta v2 (2,865 samples/2,443 indiv.)"],
  "pqtl":  ["Exposure 3  91 inflammatory-protein pQTL",
            "SCALLOP (n=14,736) + UKB-PPP (n=34,557)"],
  "mr":    ["Mendelian randomization: SMR + HEIDI + colocalization (PP.H4>0.8)",
            "+ Steiger directionality + Brion power / minimum detectable effect"],
  "causal":["Candidate-gene prioritization",
            "142 tier-1 candidate genes -> 15 tier-A candidates",
            "Alternative-eQTL robustness: 86% of testable met criteria; 99% direction concordance"],
  "null":  ["No robust support among currently instrumentable circulating proteins",
            "74 testable across two platforms (incl. TNF/CCL2/IFNG): 0 pass MR+coloc",
            "SCALLOP testable proteins: >=80% power for OR>=1.10 under assumptions; 17 lack instruments"],
  "expr":  ["State layer: 4-cohort blood meta (178/109)",
            "Inflammation-related state signal (PTX3 BH=0.038); myeloid skew",
            "76 testable causal genes: 0/76 significant (small effects, underpowered)"],
  "acu":   ["Hand-acupuncture (GSE86392): hypothesis-generating",
            "Within-study directional reversal 26/26 (descriptive)",
            "Cross-species blood 64.8% (FDR=0.14, ns)"],
 },
}

def esc(s): return s.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

def box(n, lines, bold_first=True):
    x,y,w,h=n["x"],n["y"],n["w"],n["h"]; cx=x+w/2
    s=f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="#ffffff" stroke="#000000" stroke-width="1.6"/>'
    nl=len(lines); hf,bf,lh=13.5,11.5,17.5
    total=hf+(nl-1)*lh if nl>1 else hf
    ty=y+h/2-total/2+hf*0.78
    for i,ln in enumerate(lines):
        fs=hf if i==0 else bf; fw="700" if (i==0 and bold_first) else "400"
        s+=(f'<text x="{cx:.0f}" y="{ty:.0f}" text-anchor="middle" font-size="{fs}" '
            f'font-weight="{fw}" fill="#000000">{esc(ln)}</text>')
        ty+=lh
    return s

def head(px,py,ang,size=8):
    a1=ang+math.radians(148); a2=ang-math.radians(148)
    return (f'<path d="M {px:.1f} {py:.1f} L {px+size*math.cos(a1):.1f} {py+size*math.sin(a1):.1f} '
            f'L {px+size*math.cos(a2):.1f} {py+size*math.sin(a2):.1f} Z" fill="#000000"/>')

def elbow(pts, arrow=True):
    d="M "+" L ".join(f"{px:.0f} {py:.0f}" for px,py in pts)
    s=f'<path d="{d}" fill="none" stroke="#000000" stroke-width="1.4"/>'
    if arrow:
        (x1,y1),(x2,y2)=pts[-2],pts[-1]
        s+=head(x2,y2,math.atan2(y2-y1,x2-x1))
    return s

def build(lang):
    t=TXT[lang]; N=NODES
    p=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="{FONT[lang]}">']
    p.append(f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>')
    p.append(f'<text x="{W/2:.0f}" y="40" text-anchor="middle" font-size="22" font-weight="700" fill="#000000">{esc(t["title"])}</text>')
    p.append(f'<text x="{W/2:.0f}" y="64" text-anchor="middle" font-size="13" fill="#000000">{esc(t["subtitle"])}</text>')
    for cx in (230,650,1067):
        p.append(elbow([(650,146),(650,172),(cx,172),(cx,196)]))
    for cx in (233,650,1067):
        p.append(elbow([(cx,282),(cx,330)]))
    for cx in (335,965):
        p.append(elbow([(650,388),(650,414),(cx,414),(cx,452)]))
    p.append(elbow([(335,578),(335,648)]))
    p.append(elbow([(965,578),(965,648)]))
    p.append(f'<line x1="40" y1="620" x2="1260" y2="620" stroke="#000000" stroke-width="1" stroke-dasharray="3 3"/>')
    p.append(f'<text x="{W/2:.0f}" y="638" text-anchor="middle" font-size="13" font-weight="700" fill="#000000">{esc(t["band"])}</text>')
    for k,n in N.items(): p.append(box(n,t[k]))
    p.append("</svg>")
    return "\n".join(p)

os.makedirs(OUT_DIR, exist_ok=True)
for lang,fn in (("cn","Fig_workflow_CN.svg"),("en","Fig_workflow_EN.svg")):
    with open(os.path.join(OUT_DIR,fn),"w",encoding="utf-8") as f: f.write(build(lang))
    print("wrote",fn)
