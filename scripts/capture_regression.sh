#!/usr/bin/env bash
#
# Capture golden outputs from the LEGACY solver for refactor-parity testing.
#
#   scripts/capture_regression.sh [outdir]        (default: regression/baseline)
#
# For each case a shortened config is generated (few iterations), the legacy
# binary runs it, and every produced file is archived with a SHA256 manifest.
# The legacy code resolves ../input and ../output relative to its CWD, so the
# run happens in build/ (repo/build/../input = repo/input). repo/output/* is
# cleared between cases -- do not keep precious results there when running.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PRESET="${PRESET:-release}"
BUILD="$ROOT/build/$PRESET"
RUNDIR="$ROOT/build"
OUTDIR="${1:-$ROOT/regression/baseline}"
ITER="${ITER:-200}"
OUTFREQ="${OUTFREQ:-100}"

mkdir -p "$OUTDIR"

# --- locate cmake -----------------------------------------------------------
if ! command -v cmake >/dev/null 2>&1; then
    for c in "$HOME"/.local/opt/cmake-*/bin/cmake; do
        [[ -x "$c" ]] && export PATH="$(dirname "$c"):$PATH" && break
    done
fi
command -v cmake >/dev/null 2>&1 || { echo "ERROR: cmake not found" >&2; exit 1; }

# --- configure & build ------------------------------------------------------
cmake -S "$ROOT" -B "$BUILD" -D CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"
cmake --build "$BUILD" -j "$(nproc)"
BIN="$(find "$BUILD" -maxdepth 1 -type f -name '2D_NS_Solver' -perm -u+x | head -1)"
[[ -n "$BIN" ]] || { echo "ERROR: legacy binary not built" >&2; exit 1; }

# --- config rewriter: override selected fields of a legacy .cfg -------------
rewrite_cfg() { # <src> <dst> <test_case> <mesh> <topo> <dims>
    python3 - "$1" "$2" "$3" "$4" "$5" "$6" "$ITER" "$OUTFREQ" <<'PY'
import sys
src, dst, tc, mesh, topo, dims, maxiter, outfreq = sys.argv[1:9]
KEYMAP = {
    "#Test case name": ("scalar", tc),
    "#mesh file": ("scalar", mesh),
    "#Mesh topology file containing boundary informations": ("scalar", topo),
    "#grid size (Nx*Ny)": ("dims", dims),
    "#Maximum number of iterations": ("scalar", maxiter),
    "#Output frequency. Output will be written": ("scalar", outfreq),
}
lines = open(src).read().splitlines()
out, i = [], 0
while i < len(lines):
    ln = lines[i]
    out.append(ln)
    hit = next((v for k, v in KEYMAP.items() if ln.startswith(k)), None)
    if hit and i + 1 < len(lines):
        i += 1
        out.append(" ".join(hit[1].split()) if hit[0] == "dims" else str(hit[1]))
    i += 1
open(dst, "w").write("\n".join(out) + "\n")
PY
}

# --- case table: name | source cfg | test_case | mesh | topo | dims ---------
run_case() {
    local name="$1" srccfg="$2" tc="$3" mesh="$4" topo="$5" dims="$6"
    local cfg="input/generated_${name}.cfg"
    rewrite_cfg "$ROOT/input/$srccfg" "$ROOT/$cfg" "$tc" "$mesh" "$topo" "$dims"
    rm -rf "$ROOT/output"; mkdir -p "$ROOT/output"
    echo "--- running case: $name ($dims cells/iters=$ITER)"
    ( cd "$RUNDIR" && "$BIN" "../$cfg" > "$OUTDIR/${name}.run.log.tmp" 2>&1 ) \
        || { cat "$OUTDIR/${name}.run.log.tmp"; echo "ERROR: solver failed for $name" >&2; exit 1; }
    local dest="$OUTDIR/$name"; mkdir -p "$dest"
    mv "$OUTDIR/${name}.run.log.tmp" "$dest/run.log"
    cp "$ROOT/$cfg" "$dest/generated.cfg"
    find "$ROOT/output" -type f -exec mv {} "$dest/" \;
    ( cd "$dest" && sha256sum ./*.dat ./*.vtk ./*.cfg 2>/dev/null > SHA256SUMS || true )
    echo "    archived $(ls "$dest" | wc -l) files -> $dest"
}

run_case ramp_euler_24080      InputRampEuler.cfg            RegRampEuler           RampGrid24080.dat          GridTopRampEuler24080.dat           "240 80"
run_case blasius_ns_16064      InputBlasius.cfg              RegBlasiusNS           BlasiusGrid16064.dat       GridTopBlasius16064.dat             "160 64"
run_case hypercyl_euler_80160  InputHypersonicFlowEuler.cfg  RegHyperCylEuler80160  HalfCylinderGrid80160.dat  GridTopHypersonicFlowEuler80160.dat "80 160"

rm -rf "$ROOT/output"
echo "Done. Goldens in $OUTDIR"
