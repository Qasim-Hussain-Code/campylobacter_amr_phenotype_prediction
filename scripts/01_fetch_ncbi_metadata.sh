#!/usr/bin/env bash
set -euo pipefail

PDG="PDG000000003.2859"
BASE="https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Campylobacter/${PDG}"
RAW="$(dirname "$0")/../data/raw"

mkdir -p "${RAW}"

fetch() {
    local url="$1" dest="$2"
    if [[ -s "${dest}" ]]; then
        echo "[skip] $(basename "${dest}")"
        return
    fi
    echo "[get ] $(basename "${dest}")"
    curl -fSL --retry 3 -o "${dest}.part" "${url}"
    mv "${dest}.part" "${dest}"
}

fetch "${BASE}/Metadata/${PDG}.metadata.tsv" \
      "${RAW}/${PDG}.metadata.tsv"

fetch "${BASE}/Clusters/${PDG}.reference_target.cluster_list.tsv" \
      "${RAW}/${PDG}.cluster_list.tsv"

echo
echo "[done] release ${PDG}"
wc -l "${RAW}/${PDG}.metadata.tsv" "${RAW}/${PDG}.cluster_list.tsv"
