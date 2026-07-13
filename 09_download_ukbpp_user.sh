#!/usr/bin/env bash
# ============================================================
# Download the 91 UKB-PPP inflammatory-protein pGWAS archives on the analysis host.
#
# 前置(只需一次):
#   brew install pipx && pipx install synapseclient
#   凭证用本机 ~/.synapseConfig(已配好,无需再登录)
#
# 体量:每个 .tar 约 300–550MB,全 91 个约 30GB
#   可随时 Ctrl-C 中断;重跑会自动跳过已下好的(synapse 缓存)
#   新覆盖的 40 个(含 IL6/TNF/IFNG/CCL2)排在最前,先下到最有价值的
# ============================================================
set -u
MAP="${UKBPPP_MAP:-metadata/ukbpp_coverage_map.csv}"
DEST="${UKBPPP_DEST:-data/raw/ukbpp}"

mkdir -p "$DEST"
command -v synapse >/dev/null 2>&1 || {
  echo "✗ 未找到 synapse 命令。请先安装:brew install pipx && pipx install synapseclient"
  exit 1
}
[ -f "$MAP" ] || { echo "✗ 找不到清单 $MAP"; exit 1; }

total=$(($(wc -l < "$MAP") - 1)); i=0
echo ">> 共 $total 个蛋白,下载到 $DEST"
# 第3列=ukbpp_synid;清单已按"新覆盖优先"排序
tail -n +2 "$MAP" | while IFS=, read -r gene in_ukbpp synid rest; do
  [ -z "${synid:-}" ] && continue
  i=$((i+1))
  echo ""; echo "=== [$i/$total] $gene ($synid) ==="
  synapse get "$synid" --downloadLocation "$DEST" || echo "  ! $gene 下载失败,跳过(可稍后重跑补下)"
done
echo ""; echo ">> 全部处理完。文件在:$DEST"
echo ">> Download stage complete. Continue with the analysis steps in README.md."
