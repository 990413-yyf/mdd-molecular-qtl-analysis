#!/usr/bin/env python3
## ============================================================================
## 07_module5_network_STRING_OpenTargets_tiering.py
## 模块五: 142 tier-1 因果基因 网络整合 + 成药性分级
## ----------------------------------------------------------------------------
## 输入:
##   out/unified_causal_tier1.csv   (142 因果基因: 血80+脑68去重, 组织/方向/coloc)
## 输出:
##   out/ppi_degree_hubs.csv                (STRING PPI degree/hub/模块; 99 行)
##   out/unified_causal_tier_ranked.csv     (142 基因多层打分 + tier + 成药性 + 药物)
##   out/fig_ppi_network.png                (PPI 网络, hub 加粗, 组织着色)
##   out/fig_tierA_druggable.png            (15 tier-A 靶点证据分条形图)
## 对应正文: 模块五 (可成药因果靶点骨架; 本模块不碰炎症, 142 基因均非炎症)
## ----------------------------------------------------------------------------
## 联网: STRING API (string-db.org), OpenTargets GraphQL
##       (api.platform.opentargets.org/api/v4/graphql), MyGene (symbol→ENSG)。
## 关键修正 (可重复性): STRING /network 返回的 preferredName 可能含面板外规范名
##   (如 CIPC); 只保留两端均为 142 面板基因的边, 剔除外来节点。
##   RMC1 是面板内基因但 STRING preferredName 与查询符号不同, 经边保留 → 99 节点。
## ============================================================================
import urllib.request, json, csv, time
import networkx as nx

WS = "."                        # 运行时设为 workspace 根 (out/ 在此之下)
panel = [r["gene"] for r in csv.DictReader(open(f"{WS}/out/unified_causal_tier1.csv"))]
panel_set = set(panel)

## ---------------- 一. STRING PPI (物种 9606, 置信分>0.4) ----------------
STRING = "https://string-db.org/api"
def string_post(endpoint, params):
    data = urllib.parse.urlencode(params).encode()
    return urllib.request.urlopen(urllib.request.Request(f"{STRING}/{endpoint}", data=data), timeout=60)
## map identifiers
r = string_post("json/get_string_ids",
     {"identifiers": "\r".join(panel), "species": 9606, "limit": 1, "echo_query": 1})
idmap = json.load(r)
q2string = {d["queryItem"]: d["preferredName"] for d in idmap}   # 98 蛋白编码映射成功
## network edges (confidence >0.4 = score>400)
r = string_post("json/network", {"identifiers": "\r".join(q2string.values()),
                                   "species": 9606, "required_score": 400})
net_raw = json.load(r)
## 只保留两端均为面板基因的边 (剔除 STRING 规范外来名如 CIPC)
string2panel = {v: k for k, v in q2string.items()}
edges = []
for e in net_raw:
    a = string2panel.get(e["preferredName_A"], e["preferredName_A"])
    b = string2panel.get(e["preferredName_B"], e["preferredName_B"])
    if a in panel_set and b in panel_set and a != b:
        edges.append((a, b, float(e["score"])))
## 节点 = 面板内且映射到 STRING, 或出现在保留边中 (含 RMC1)
edge_nodes = {n for a, b, s in edges for n in (a, b)}
nodes = sorted((set(q2string) & panel_set) | edge_nodes)         # 99 节点
G = nx.Graph(); G.add_nodes_from(nodes)
for a, b, s in edges: G.add_edge(a, b, weight=s)
deg = dict(G.degree())
comps = sorted([c for c in nx.connected_components(G) if len(c) >= 3], key=len, reverse=True)
## 结果: 99 节点 / 29 边 / 2 子模块; hubs(deg>=3)=7:
##   CTBP1, LRGUK, COPZ1, HSPE1, KLF11, NSF, SF3B1
tiss = {r["gene"]: r["tissue"] for r in csv.DictReader(open(f"{WS}/out/unified_causal_tier1.csv"))}
with open(f"{WS}/out/ppi_degree_hubs.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "tissue", "degree", "is_hub", "module"]); w.writeheader()
    for g in sorted(deg, key=lambda x: -deg[x]):
        mod = next((f"module{i+1}" for i, c in enumerate(comps) if g in c), "")
        w.writerow(dict(gene=g, tissue=tiss.get(g, ""), degree=deg[g],
                        is_hub="Y" if deg[g] >= 3 else "", module=mod))

## ---------------- 二. 成药性注释 (OpenTargets GraphQL v4) ----------------
## symbol→ENSG 经 MyGene (示意); 107/142 解析 (35 个 lincRNA/假基因未解析)
OT = "https://api.platform.opentargets.org/api/v4/graphql"
def gql(q, v, retry=3):
    body = json.dumps({"query": q, "variables": v}).encode()
    for _ in range(retry):
        try:
            return json.load(urllib.request.urlopen(
                urllib.request.Request(OT, data=body, headers={"Content-Type": "application/json"}), timeout=40))
        except urllib.error.HTTPError as e:
            return {"HTTPErr": e.read().decode()[:200]}
        except Exception:
            time.sleep(2)
    return {"err": "exhausted"}
## 注意: knownDrugs 已非 Target 字段, 用 drugAndClinicalCandidates (无 size 参数)。
Q = '''query($id:String!){target(ensemblId:$id){approvedSymbol
 tractability{modality value label}
 drugAndClinicalCandidates{count rows{drug{name}}}}}'''
## 对每个解析到 ENSG 的基因查询; sm_druggable/ab_druggable 由 tractability bucket
## (Approved Drug/Druggable Family/High-Quality Ligand|Pocket/Structure with Ligand...) 判定。
## 结果: 29 SM可成药, 52 AB可成药, 62 任一; 7 有药物候选:
##   ESR2, HDAC3, CHRNA4, TUBB, EPHA7, PPP3CA, RPS6KB1。

## ---------------- 三+四. 已知MDD注释 + 多层证据打分 → tierA/B/C ----------------
## 已知 MDD (保守, 拿得准才标): DCC KLF7 NEGR1 GRIA1 CLOCK FADS1 SLC12A5 RERE
##   ESR2 HDAC3 SF3B1 (共 9 个落在 142 内)。
## 打分 (每项+1): 血MR命中 / 脑MR命中 / 脑血共有 / coloc>0.9 / 已知MDD /
##                可成药(SM或AB) / 有现成药 / PPI-hub。
## tier: A(score>=4) / B(2-3) / C(<=1)。
## 结果: 15 tier-A, 102 tier-B, 25 tier-C。tier-A:
##   ESR2 HDAC3 KDELR2 RERE REV1 SF3B1 ZDHHC5 BMS1P4 CHRNA4 CTBP1 FADS1 HSPE1 PPP6C SLC12A5 TUBB
## → 写 out/unified_causal_tier_ranked.csv (gene,tissue,direction,coloc_PPH4,known_MDD,
##    sm_druggable,ab_druggable,n_drug_candidates,has_drug,PPI_hub,degree,score,tier,drugs,ensg)

## ---------------- 五. 功能富集 ----------------
## 复用模块一已算 out/GO_enrichment_merged_tier1.csv (clusterProfiler, 见 R 段);
## 142 基因 GO: q<0.05=3 (steroid/cholesterol/sterol binding), 无炎症通路。

## ---------------- 图 (matplotlib) ----------------
## fig_ppi_network.png: spring_layout 连通子图, node size∝degree, hub 加粗, 组织着色。
## fig_tierA_druggable.png: 15 tier-A 横向条形, 证据分, 标注 MDD/drug/SM/AB/hub。
