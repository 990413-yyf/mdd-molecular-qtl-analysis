#!/usr/bin/env bash
# 只补下 UKBPPP_tmp 里还缺的 3 个蛋白(在你自己 Mac 上跑)
# 凭证用本机 ~/.synapseConfig;失败会跳过,可重复跑
set -u
DEST="${UKBPPP_DEST:-data/raw/ukbpp}"
mkdir -p "$DEST"
command -v synapse >/dev/null 2>&1 || { echo "先装: brew install pipx && pipx install synapseclient"; exit 1; }
echo ">> Downloading three previously missing proteins"
echo "=== [1/3] IL17C (syn51468579) ==="; synapse get syn51468579 --downloadLocation "$DEST" || echo "  ! IL17C 失败,跳过"
echo "=== [2/3] IL18R1 (syn51469422) ==="; synapse get syn51469422 --downloadLocation "$DEST" || echo "  ! IL18R1 失败,跳过"
echo "=== [3/3] KITLG (syn51468698) ==="; synapse get syn51468698 --downloadLocation "$DEST" || echo "  ! KITLG 失败,跳过"
echo ">> 完成。再看还缺不缺,或回对话说\"继续\"。"
