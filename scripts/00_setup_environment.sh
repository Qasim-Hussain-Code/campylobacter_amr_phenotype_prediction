#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="campy-amr"
ENV_FILE="$(dirname "$0")/../envs/environment.yml"

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "[info] ${ENV_NAME} already exists, updating from spec"
    conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
    echo "[info] creating ${ENV_NAME}"
    conda env create -f "${ENV_FILE}"
fi

echo
echo "[done] activate with: conda activate ${ENV_NAME}"
