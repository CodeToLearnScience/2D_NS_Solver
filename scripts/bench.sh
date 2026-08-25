#!/usr/bin/env bash
#
# Benchmark: wall-clock per-iteration timings for the ramp-Euler case.
#   - legacy binary (frozen reference)
#   - new stack, serial (ns_solver without mpirun)
#   - new stack, MPI at 2 and 4 ranks
#
#   scripts/bench.sh [iterations]     (default 300)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ITERS="${1:-300}"

LEG_BIN=""
for cand in "$ROOT/build/release/2D_NS_Solver" "$ROOT/bin/2D_NS_Solver"; do
    [[ -x "$cand" ]] && LEG_BIN="$cand" && break
done
NEW=""
for cand in "$ROOT/build/release/ns_solver" "$ROOT/build/ns_solver" \
            "$ROOT/build-mpi/ns_solver"; do
    [[ -x "$cand" ]] && NEW="$cand" && break
done
[[ -n "$NEW" ]] || { echo "no ns_solver built"; exit 1; }
MPISOLVER="$([[ -x "$ROOT/build-mpi/ns_solver" ]] && echo "$ROOT/build-mpi/ns_solver")"

bench() { # <label> <cmd...>
    local label="$1"; shift
    local t0 t1
    t0=$(date +%s.%N)
    "$@" >/dev/null 2>&1
    t1=$(date +%s.%N)
    python3 -c "print(f'  {\"$label\":28s} {$t1-$t0:7.2f} s  ({$ITERS/($t1-$t0):8.1f} iter/s)')"
}

echo "benchmark: ramp-Euler, $ITERS iterations"
rm -rf "$ROOT/output" && mkdir -p "$ROOT/output"

python3 - "$ROOT/input/InputRampEuler.cfg" "$ROOT/input/generated_bench.cfg" Bench \
          RampGrid24080.dat GridTopRampEuler24080.dat "240 80" "$ITERS" $((ITERS + 1)) <<'PY'
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

if [[ -n "$LEG_BIN" ]]; then
    bench "legacy" bash -c "cd '$ROOT/build' && '$LEG_BIN' ../input/generated_bench.cfg"
fi
OUT=$(mktemp -d)
bench "new serial" "$NEW" configs/RampEuler.toml input/RampGrid24080.dat "$OUT" "$ITERS"
if [[ -n "$MPISOLVER" ]]; then
    for NP in 2 4; do
        bench "new mpi -np $NP" mpirun --oversubscribe --allow-run-as-root -np $NP \
              "$MPISOLVER" configs/RampEuler.toml input/RampGrid24080.dat "$OUT" "$ITERS"
    done
fi
rm -rf "$OUT" "$ROOT/output" "$ROOT/input/generated_bench.cfg"
