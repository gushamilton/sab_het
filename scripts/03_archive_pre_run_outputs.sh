#!/bin/bash
set -euo pipefail

# Archive outputs/logs older than a cutoff timestamp to avoid run provenance confusion.
# Usage:
#   bash scripts/03_archive_pre_run_outputs.sh 2026-02-11T10:35:08 pre_highres_15612319

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <cutoff-iso8601> [archive_label]"
  exit 1
fi

CUTOFF="$1"
LABEL="${2:-archive_$(date +%Y%m%d_%H%M%S)}"
ARCHIVE_ROOT="archive/${LABEL}"
MANIFEST="${ARCHIVE_ROOT}/manifest_moved_files.txt"

mkdir -p "${ARCHIVE_ROOT}/results" "${ARCHIVE_ROOT}/logs"

echo "Archive label: ${LABEL}"
echo "Cutoff: ${CUTOFF}"
echo "Archive root: ${ARCHIVE_ROOT}"

move_tree() {
  local src_root="$1"
  local dst_root="$2"

  if [[ ! -d "${src_root}" ]]; then
    echo "Skip missing directory: ${src_root}"
    return 0
  fi

  local tmp_list
  tmp_list="$(mktemp)"
  find "${src_root}" -type f ! -newermt "${CUTOFF}" | sort > "${tmp_list}"
  local n_files
  n_files="$(wc -l < "${tmp_list}")"
  echo "Files to move from ${src_root}: ${n_files}"

  if [[ "${n_files}" -eq 0 ]]; then
    rm -f "${tmp_list}"
    return 0
  fi

  while IFS= read -r src_file; do
    rel_path="${src_file#${src_root}/}"
    dst_file="${dst_root}/${rel_path}"
    mkdir -p "$(dirname "${dst_file}")"
    mv "${src_file}" "${dst_file}"
    echo "${src_file} -> ${dst_file}" >> "${MANIFEST}"
  done < "${tmp_list}"

  rm -f "${tmp_list}"
}

: > "${MANIFEST}"
move_tree "results" "${ARCHIVE_ROOT}/results"
move_tree "logs" "${ARCHIVE_ROOT}/logs"

echo "Archive complete. Manifest: ${MANIFEST}"
