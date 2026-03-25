#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

: "${BIOWULF_HOST:?Set BIOWULF_HOST, for example obriena2@biowulf.nih.gov}"
: "${BIOWULF_REPO_DIR:?Set BIOWULF_REPO_DIR, for example ~/scScalpChromatin}"

if [[ "${BIOWULF_REPO_DIR}" == /Users/* ]]; then
  echo "BIOWULF_REPO_DIR looks like a local macOS path: ${BIOWULF_REPO_DIR}" >&2
  echo "Use a remote path instead, for example '~/scScalpChromatin' or '/data/\$USER/scScalpChromatin'." >&2
  exit 1
fi

printf -v REMOTE_REPO_DIR_Q '%q' "${BIOWULF_REPO_DIR}"
ssh "${BIOWULF_HOST}" "mkdir -p ${REMOTE_REPO_DIR_Q}"

rsync -av \
  --exclude ".git/" \
  --exclude ".DS_Store" \
  --exclude "._.DS_Store" \
  "${ROOT_DIR}/" \
  "${BIOWULF_HOST}:${BIOWULF_REPO_DIR}/"

echo "Synced ${ROOT_DIR} to ${BIOWULF_HOST}:${BIOWULF_REPO_DIR}"
