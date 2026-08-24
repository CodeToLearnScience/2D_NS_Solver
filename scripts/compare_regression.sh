#!/usr/bin/env bash
#
# Byte-compare two regression directories (e.g. baseline vs post-refactor run).
#
#   scripts/compare_regression.sh regression/baseline regression/new
#
# Compares every file listed in SHA256SUMS of DIR_A against DIR_B.
# run.log is skipped (banner output may legitimately change during refactor).
# Exits non-zero on any mismatch or missing file.
set -euo pipefail

A="${1:?usage: compare_regression.sh <dirA> <dirB> [case ...]}"
B="${2:?usage: compare_regression.sh <dirA> <dirB> [case ...]}"
shift 2 || true
CASES=("$@")

fail=0
compare_case() {
    local ca="$1/$2" cb="$3/$2"
    if [[ ! -f "$ca/SHA256SUMS" ]]; then echo "MISSING manifest: $ca/SHA256SUMS"; fail=1; return; fi
    if [[ ! -d "$cb" ]]; then echo "MISSING result dir: $cb"; fail=1; return; fi
    local status=0
    while read -r sum path; do
        local fname="${path#./}"
        [[ "$fname" == "run.log" || "$fname" == "SHA256SUMS" ]] && continue
        if [[ ! -f "$cb/$fname" ]]; then
            echo "[$2] MISSING file: $fname"; status=1
        elif ! echo "$sum  $cb/$fname" | sha256sum -c --quiet >/dev/null 2>&1; then
            echo "[$2] MISMATCH:     $fname"; status=1
        fi
    done < "$ca/SHA256SUMS"
    if (( status == 0 )); then echo "[$2] OK ($(wc -l < "$ca/SHA256SUMS") files match)"; fi
    fail=$((fail | status))
}

if (( ${#CASES[@]} == 0 )); then
    mapfile -t CASES < <(find "$A" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
fi

for c in "${CASES[@]}"; do compare_case "$A" "$c" "$B"; done

if (( fail )); then
    echo "REGRESSION CHECK FAILED"
    exit 1
fi
echo "All compared cases match."
