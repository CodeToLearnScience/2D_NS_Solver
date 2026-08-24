#!/usr/bin/env bash
#
# MPI end-to-end check: builds the MPI flavour (if needed) and runs the ramp
# case at 1/2/4 ranks, reporting residual agreement vs the LEGACY solver.
#
#   scripts/check_mpi_parity.sh [iterations]     (default 200)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ITERS="${1:-200}"
BUILD="$ROOT/build-mpi"

if [[ ! -x "$BUILD/ns_solver" ]]; then
    cmake -S "$ROOT" -B "$BUILD" -D CMAKE_BUILD_TYPE=Release -D NS_ENABLE_MPI=ON \
          -D NS_BUILD_TESTS=OFF >/dev/null
    cmake --build "$BUILD" -j "$(nproc)" >/dev/null
fi

# legacy reference (fixed sources)
python3 - "$ROOT/input/InputRampEuler.cfg" "$ROOT/input/generated_mpi.cfg" RegRampEuler \
          RampGrid24080.dat GridTopRampEuler24080.dat "240 80" "$ITERS" $((ITERS)) <<'PY'
import sys
src,dst,tc,mesh,topo,dims,mi,of = sys.argv[1:9]
K={"#Test case name":("s",tc),"#mesh file":("s",mesh),
   "#Mesh topology file containing boundary informations":("s",topo),
   "#grid size (Nx*Ny)":("d",dims),"#Maximum number of iterations":("s",mi),
   "#Output frequency. Output will be written":("s",of)}
lines=open(src).read().splitlines(); out=[]; i=0
while i<len(lines):
    ln=lines[i]; out.append(ln)
    hit=next((v for k,v in K.items() if ln.startswith(k)),None)
    if hit and i+1<len(lines):
        i+=1
        out.append(" ".join(hit[1].split()) if hit[0]=="d" else str(hit[1]))
    i+=1
open(dst,"w").write("\n".join(out)+"\n")
PY
rm -rf "$ROOT/output" && mkdir -p "$ROOT/output"
LEG_BIN=""
for cand in "$ROOT/build/release/2D_NS_Solver" "$BUILD/ns_legacy"             "$ROOT/bin/2D_NS_Solver"; do
    [[ -x "$cand" ]] && LEG_BIN="$cand" && break
done
[[ -n "$LEG_BIN" ]] || { echo "legacy binary not found"; exit 1; }
(cd "$BUILD" && "$LEG_BIN" ../input/generated_mpi.cfg >/dev/null 2>&1)
LEG=$(ls "$ROOT/output"/Residue_*.dat 2>/dev/null | head -1)
[[ -n "$LEG" ]] || { echo "legacy reference run failed"; exit 1; }

for NP in 1 2 4; do
    OUT=$(mktemp -d)
    mpirun --allow-run-as-root -np $NP "$BUILD/ns_solver" \
        configs/RampEuler.toml input/RampGrid24080.dat "$OUT" "$ITERS" >/dev/null 2>&1
    NEW="$OUT/Residue_Ramp_KFDS_maxU_24080_2nd240_80.dat"
    python3 - "$LEG" "$NEW" "$NP" <<'PY'
import sys
def load(p):
    return [[float(x) for x in ln.split()] for ln in open(p)
            if ln.split() and ln.split()[0].isdigit()]
L,N=load(sys.argv[1]),load(sys.argv[2])
wdt=max(abs(a[1]-b[1]) for a,b in zip(L,N))
print(f"np={sys.argv[3]}: max |dt-legacy| = {wdt:.2e}, "
      f"final res leg={L[-1][3]:.5e} new={N[-1][3]:.5e}")
PY
    rm -rf "$OUT"
done
rm -rf "$ROOT/output" "$ROOT/input/generated_mpi.cfg"
echo "done."
